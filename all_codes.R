#package required for simulation
library(SimDesign)
library(R.matlab)
library(Matrix)
##### ALGORITHM FUNTIONS #####

#function 1: Cool inverse that guarantees symmetricity
#input: matrix want to invert(mat)
#output: inverted matrix that is guaranteed to be symmetric
cool_inverse=function(mat){
  inv_mat=solve(mat)
  inv_mat=(inv_mat+t(inv_mat))/2
  return(inv_mat)
}

#function 2: Make matrix positive definite by adding quantity on diagonal, thanks Dr. Rothman for the code!!
#inpute: a symmetric matrix(matrix), extra add-on number for all diagonal elements(addon)
#output: positive definite matrix
pd_matrix=function(matrix,addon=.2){
  p=dim(matrix)[1]
  eigens=eigen(matrix)
  ev=eigens$values
  dd=(p*min(ev)-max(ev))/(1-p)
  matrix=matrix+diag(dd+addon,p)
  return(matrix)
}

#function 3: soft-threshold function
#input: variable(x), threshold(t)
#output: shresholded x
soft_thresh=function(x,t){
  sign(x)*ifelse(abs(x)-t>0,abs(x)-t,0)
}

#function 4: pathwise coordinate descent
#input: symmetric p*p matrix(W), p dim vector(s), p dim penalty vector or scalar(rho), whether print information for debugging(print_iter)
#output: optimization dual parameter which is a middle step value(beta0)
coord_descent=function(W,s,rho,print_iter=FALSE){
  if(!all(W == t(W))){
    print('Error: W is not symmetrix!')
    return(NULL)
  }
  if(length(rho)==1 && length(s)>1){
    rho=rep(rho,length(s))
  }
  threshold=1e-6
  dis=1
  beta0=1/diag(W)*s/3
  beta1=rep(0,length(s))
  iter=0
  max_iter=5e3
  while (dis>threshold && iter<max_iter) {
    iter=iter+1
    beta00=beta0
    for (j in (1:length(s))) {
      x=s[j]-sum(as.vector(W[j,-j])*beta0[-j])
      beta1[j]=soft_thresh(x,rho[j])/W[j,j]
      beta0[j]=beta1[j]
    }
    dis=sqrt(sum((beta00-beta1)^2))
    beta0=beta1
    if(dis==0){
      break
    }
    if(print_iter==T){
      cat('iterate(',iter,'), update distance:',dis,'\n')
    }
  }
  return(beta0)
}

#function 5: glasso_costm: our replication of glasso in glasso package
#input: sample covariance matrix(S), penalty matrix or scalar(rho),whether print information for debugging(print_iter)
#output: estimated covariance matrix that is sparse(W), estimated precision matrix that is sparse(WI)
glasso_costm=function(S,rho,print_iter=F){
  p=dim(S)[1]
  if(length(rho)==1&&p>1){
    rho=matrix(rep(rho,p*p),ncol = p)
  }
  threshold=mean(as.vector(abs(S-diag(diag(S)))))*1e-4
  W0=S+diag(diag(rho))
  W1=W0
  dis=1
  iter=0
  max_iter=10000
  disturb_on=F
  while(dis>threshold && iter<max_iter){
    iter=iter+1
    W1=W0
    for (j in 1:p) {
      W11=W1[-j,-j]
      s12=as.vector(S[j,-j])
      rho0=as.vector(rho[j,-j])
      beta=coord_descent(W11,s12,rho0)
      w12=as.vector(W11%*%beta)
      W1[j,-j]=w12
      W1[-j,j]=w12
    }
    dis=mean(as.vector(abs(W1-W0)))
    W0=W1
    if(print_iter==T){
      cat('iterate(',iter,'), update distance:',dis,'\n')
    }
    if(iter>1500 && disturb_on==F){
        W0=W1=W0+diag(diag(rho*0.1))
        disturb_on=T
    }
  }
  WI=cool_inverse(W1)
  WI=ifelse(abs(WI)<1e-6,0,WI)
  W=W1
  return(list(W=W,WI=WI))
}


#function 6: middle step for joint Glasso: calculate the entrywise penalty parameter given estimation of precision matrices list in a single iteration
#input: estimated precision matrices list in an iteration(Omega_list)
#output: corresponding penalty matrix(tao)
Cal_weight=function(Omega_list){
  K=length(Omega_list)
  lis=lapply(Omega_list,function(x){x-diag(diag(x))+diag(1,dim(Omega_list[[1]])[1])})
  tao=abs(lis[[1]])
  for (k in 2:K) {
    tao=tao+abs(lis[[k]])
  }
  tao=ifelse(tao<10e-10,10e-10,tao)
  tao=tao^(-1/2)
  return(tao)
}



#function 7: Cool Glasso: joint Glasso function
#input: covariance matrices list of all categories(S_list), penalty term parameter(lambda)
#output: estimated precision matrices list that are sparse
cool_glasso=function(S_list,lambda,nu_tune=0.5){
  K=length(S_list)
  threshold=1e-4
  dis=1
  iter=0
  max_iter=1e4
  p=dim(S_list[[1]])[1]
  Omega0=lapply(S_list, function(x){
    evl=eigen(x)$values
    dd=(p*min(evl)-max(evl))/(1-p)
    ci=cool_inverse(x + dd*diag(1,p))
    return(ci)
  })
  Omega1=Omega0
  tao0=Cal_weight(Omega0)
  tao1=Cal_weight(Omega0)+100
  #print(tao)
  while(dis>threshold && iter<max_iter){
    cat('iter=',iter,'\n')
    iter=iter+1
    diag_nonsgl_set=0
    if(iter<20){
      diag_nonsgl_set=0.2
    }
    dis=0
    for (k in 1:K) {
      #print(k)
      Omega1[[k]]=glasso_costm(cool_inverse(Omega0[[k]]),lambda*tao0)$WI
      dis=dis+sum(abs((Omega1[[k]]-diag(diag(Omega1[[k]])))-(Omega0[[k]]-diag(diag(Omega0[[k]])))))
      
    }
    #dis=dis/(K*p*(p-1))
    tao1=Cal_weight(Omega1)
    dis=sum(abs(tao1-tao0))/sum(abs(tao0))
    tao0=tao1
    #calculate off diagonal distance
    # if(iter>0){
    #   cat('dis=',dis,'\n')
    # }
    Omega0=Omega1
  }
  return(Omega1)
}





##### PARAMETER TUNING FUNCTIONS #####

#function 8: Middle step of parameter tuning for lambda in cool_glasso using cross-validation, calculating predictive negative log likelihood of a single fold
#input: covariance matrices calculated from one fold(S_list), estimated precision matrices from rest folds
#output: predictive negative log likelihood for single fold
cv_function=function(S_list,est_omega){
  #S_list=sigma list, est_omega=glasso list
  K=length(S_list)
  fest=0
  for(i in (1:K)){
    ff=sum(diag(S_list[[i]]%*%est_omega[[i]]))-log(det(est_omega[[i]]))
    fest=fest+ff
  }
  return(fest)
}

#function 9: Middle step of parameter tuning for lambda in cool_glasso using cross-validation, calculating average predictive negative log likelihood of all folds
#input: sample matrices list in all categories(X_list), how many folds for CV(K_fold), penalty parameter value(lambda), whether we do it for joint estimation result(joint)
#output: average predictive negative log likelihood of all folds(the larger, the better)
cv_estimation=function(X_list,k_fold,lambda,joint=T){
  #X=n*p matrix
  K=length(X_list)
  KK=k_fold
  n=nrow(X_list[[1]])
  idx=sample(1:n,n,replace=F)
  idx2=rep(1:KK,length=n)
  fcv=0
  for(i in (1:KK)){
  siglist=lapply(X_list,function(x){
    tridx=idx[idx2==i]
    testx=x[tridx,]
    sigcov=cov(testx)
    return(sigcov)#scov used for sigma, wcov used for glasso
   }
  )
  wlist=lapply(X_list,function(x){
    tridx=idx[idx2==i]
    trainx=x[-tridx,]
    wcov=cov(trainx)
    return(wcov)#scov used for sigma, wcov used for glasso
    }
   )
  if(joint==T){
    cg=cool_glasso(wlist,lambda)
  }else{
    cg=list()
    for(k in 1:K){
      cg[[k]]=glasso_costm(wlist[[k]],lambda)$WI  
    }
    
  }
  tfcv=cv_function(siglist,cg)
  fcv=fcv+tfcv
  }
  return(fcv)
}

# function 10: calculate BIC value
# input: sample covariance matrices list(S_list), estimated precision matrices list(est_omega), sample size vector(ns)
# output: bic value
BIC_custom=function(S_list,est_omega,ns){
  df=numeric(0)
  K=length(S_list)
  p=dim(S_list[[1]])[1]
  bics=0
  df_term=0
  for (k in 1:K) {
    df[k]=(sum(est_omega[[k]]!=0)-p)/2
    bics=bics+sum(diag(S_list[[k]]%*%est_omega[[k]]))-log(det(est_omega[[k]]))+df[k]*log(ns[k])/ns[k]
    df_term=df_term+df[k]*log(ns[k])/ns[k]
  }
  print(df_term)
  return(bics)
}


##### METRICS FUNCTIONS #####
#below are 5 metric functions evaluating how good an estimation is

#false negative rate
#if high, sparsity is over-estimated
cal_FNR=function(res,true){
  K=length(res)
  p=dim(true[[1]])[1]
  numerator=c()
  denominator=c()
  for (k in 1:K) {
    numerator[k]=sum(res[[k]]==0)-sum(res[[k]]==true[[k]])
    denominator[k]=sum(true[[k]]!=0)-p
  }
  sum(numerator)/sum(denominator)
}

#false positive rate
#if high, sparsity is under estimated
cal_FPR=function(res,true){
  p=dim(true[[1]])[1]
  K=length(true)
  numerator=c()
  denominator=c()
  for (k in 1:K) {
    numerator[k]=sum( (res[[k]]!=0)*(true[[k]]==0) )
    denominator[k]=sum(true[[k]]==0)
  }
  return(sum(numerator/denominator)/K)
}

#entropy loss
#comes from KL divergence, the lower the better
cal_EL=function(res,true){
  #true:true omega, res:estimated omega
  p=dim(true[[1]])[1]
  K=length(true)
  el=0
  for(k in (1:K)){
    tel=sum(diag(cool_inverse(true[[k]]) %*% res[[k]])) -
      log(det(cool_inverse(true[[k]]) %*% res[[k]]))
    el=el+tel
  }
  return((el-p)/K)
}

#Frobenius loss
#the lower the better
cal_FL=function(res,true){
  #true:true omega, res:estimated omega
  K=length(true)
  fl=0
  for(k in(1:K)){
    num=sum((true[[k]]-res[[k]])^2)
    denom=sum((true[[k]])^2)
    tfl=num/denom
    fl=fl+tfl
  }
  return(fl/K)
}


#common zero error rate
#the proportion of common zeros that are not detected, the lower the better.
cal_CZ=function(res,true){
  p=dim(res[[1]])[1]
  K=length(res)
  res_ary=array(0,dim=c(p,p,K))
  true_ary=array(0,dim=c(p,p,K))
  for(k in (1:K)){
    res_ary[,,k]=abs(res[[k]])
    true_ary[,,k]=abs(true[[k]])
  }
  rres=apply(res_ary,MARGIN=c(1,2),sum)
  ttrue=apply(true_ary,MARGIN=c(1,2),sum)
  num=sum((ttrue==0)*(rres!=0))
  denom=sum(ttrue==0)
  return(num/denom)
}


##### SIMULATION GENERATING FUNCTION #####

#function 11: middle step for generating chain case, generates common structure as tridiagonal matrix.
#input: dimension(p), random seed for generating common structure(seed)
#output: generated true covariance matrix(covmat), true precision matrix(premat)
generate_case1=function(p,seed=1234){
  set.seed(seed)
  s=numeric(0)
  s[1]=runif(1,0,1/2)
  for (i in 2:p) {
    s[i]=s[i-1]+runif(1,1/2,1)
  }
  covmat=matrix(0,nrow = p,ncol = p)
  for (i in 1:p) {
    for (j in i:p) {
      covmat[i,j]=covmat[j,i]=exp(-abs(s[i]-s[j])/2)
    }
  }
  precmat=cool_inverse(covmat)
  precmat=ifelse(abs(precmat)<1e-6,0,precmat)
  set.seed(NULL)
  return(list(covmat=covmat,precmat=precmat))
}

#function 12: generating chain case for a single category, adding category-specific structure to common structure
#input: number of samples in each category(n), proportion of category-specific to common structure(rho), random seed for generating common structure(common.seed)
#output: generated true precision matrix
generate_chain=function(n=100,p=100,rho=0.25,common.seed=111){
  M=p-1
  num=as.integer(M*rho)
  print(num)
  mat1=generate_case1(p,common.seed)$precmat
  if(rho!=0){
    for (k in 1:num) {
      i=sample(1:p,1)
      j=sample(1:p,1)
      while (mat1[i,j]!=0) {
        j=sample(1:p,1)
      }
      non0val=runif(1,-1,1)
      while(non0val<0.5&&non0val>-0.5){
        non0val=runif(1,-1,1)
      }
      mat1[i,j]=mat1[j,i]=non0val
    }
  }
  #mat1 is the maybe non pd version of true omega for a single category
  return(mat1)
}


#function 13: generating chain case for a all categories, adding category-specific structure to common structure
#input: number of samples in each category(n), proportion of category-specific to common structure(rho), random seed for generating common structure(common.seed), number of categories(K)
#output: generated true precision matrices list(true_omega_list), true covariance matrices list(S_list), sample matrices list in all categories(X_list)
generate_joint_chain=function(n=100,p=100,rho=0.25,common.seed=111,K=3){
  unmdfd_omega_list=list()
  for (k in 1:K) {
    unmdfd_omega_list[[k]]=generate_chain(n,p,rho,common.seed)
  }
  dd_list=lapply(unmdfd_omega_list, function(x){
    eigens=eigen(x)
    ev=eigens$values
    dd=(p*min(ev)-max(ev))/(1-p)
    return(dd+0.2)
  })
  dd=max(unlist(dd_list))
  true_omega_list=lapply(unmdfd_omega_list, function(x){
    x+diag(dd,p)
  })
  S_list=list()
  X_list=list()
  for (k in 1:K) {
    X_list[[k]]=rmvnorm(n,sigma = cool_inverse(true_omega_list[[k]]))
    S_list[[k]]=cov(X_list[[k]])
  }
  return(list(true_omega_list=true_omega_list,S_list=S_list,X_list=X_list))
}


#function 14: generate common structure with knn
#input: number of samples in each category(n), proportion of category-specific to common structure(rho), number of most near points to connect(m), random seed for generating common structure(common.seed)
#output: common structure precision matrix
generate_knn=function(n=100,p=100,rho=0.25,m=5,common.seed=111){
  set.seed(common.seed)
  M=p-1
  num=as.integer(M*rho)
  x=runif(p,0,1)
  y=runif(p,0,1)
  xy=cbind(x,y)
  idx=seq(1,p)
  dmat=matrix(0,p,p)
  for(i in (1:p)){
    for(j in (1:p)){
      dmat[i,j]=dist_function(xy[i,],xy[j,])
    }
  }
  diag(dmat)=rep(max(dmat)+1,p)
  sp_list=list()
  for(i in (1:p)){
    hd=head(dmat[order(dmat[,i]),i],m)
    sp_list[[i]]=idx[dmat[,i] %in% hd]
  }
  knn_mat=matrix(0,p,p)
  for(j in (1:p)){
    knn_mat[j,sp_list[[j]]]=1
  }
  knn_mat=(knn_mat+t(knn_mat))/2
  diag(knn_mat)=1
  set.seed(NULL)
  if(rho!=0){
     for(k in (1:num)){
       i=sample(1:p,1)
       j=sample(1:p,1)
     while (knn_mat[i,j]!=0) {
         j=sample(1:p,1)
       }
       non0val=runif(1,-1,1)
     while(non0val<0.5&&non0val>-0.5){
         non0val=runif(1,-1,1)
      }
       knn_mat[i,j]=knn_mat[j,i]=non0val
     }
  }
  return(knn_mat)
}

#function 15: generating knn case for a all categories, adding category-specific structure to common structure
#input: number of samples in each category(n), proportion of category-specific to common structure(rho), number of most near points to connect(m), random seed for generating common structure(common.seed), number of categories(K)
#output: generated true precision matrices list(true_omega_list), true covariance matrices list(S_list), sample matrices list in all categories(X_list)
generate_joint_knn=function(n=100,p=100,rho=0.25,m=5,common.seed=111,K=3){
  unmdfd_omega_list=list()
  for(k in (1:K)){
    unmdfd_omega_list[[k]]=generate_knn(n,p,rho,
                                        m,common.seed)
  }
  # print('here')
  levl=lapply(unmdfd_omega_list,function(x){
    eigens=eigen(x)
    ev=eigens$values
    dd=(p*min(ev)-max(ev))/(1-p)
    return(dd+0.2)
  })
  mevl=max(unlist(levl))
  true_omega_list=lapply(unmdfd_omega_list,function(x){
    x+diag(mevl,p)
  })
  S_list=list()
  X_list=list()
  for(k in (1:K)){
    X_list[[k]]=rmvnorm(n,sigma=cool_inverse(true_omega_list[[k]]))
    S_list[[k]]=cov(X_list[[k]])
  }
  return(list(true_omega_list=true_omega_list,
              S_list=S_list,
              X_list=X_list))
}


##### RUN SIMULATION-----CHAIN CASE #####

#generate
rho=0
sample_num=50
p=30
generate_result=generate_joint_chain(n=sample_num,p=p,rho=rho,common.seed=111,K=3)
true_omega_list=generate_result$true_omega_list
S_list=generate_result$S_list
X_list=generate_result$X_list
ns=rep(sample_num,length(S_list))


#parameter tuning for separate estimation case
print('start tuning for SEPARATE ESTIMATION')
lambda_seq1=10^(seq(-9,1,1))
cvs=numeric(0)
FNRs=numeric(0)
BICs=numeric(0)
ELs=numeric(0)
FLs=numeric(0)
FPRs=numeric(0)
CZs=numeric(0)
for (i in 1:length(lambda_seq1)) {
  result=list()
    for(k in (1:length(S_list))){
      result[[k]]=glasso_costm(S_list[[k]],lambda_seq1[i])$WI
     }
  cvs[i]=cv_estimation(X_list,5,lambda_seq1[i],joint=F)
  FNRs[i]=cal_FNR(result,true_omega_list)
  BICs[i]=BIC_custom(S_list,result,ns)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}
cat('lambdas:',lambda_seq1,'\n')
cat('cv:',cvs,'\n')
cat('bic:',BICs,'\n')
cat('fnr:',FNRs,'\n')
cat('fpr:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')

#finer grid search with CV
lambda_seq2=c((5:9)*lambda_seq1[which.min(cvs)]/10,(1:6)*lambda_seq1[which.min(cvs)])
cv1s=numeric(0)
for (i in 1:length(lambda_seq2)) {
  result=list()
    for(k in (1:length(S_list))){
      result[[k]]=glasso_costm(S_list[[k]],lambda_seq2[i])$WI
     }
  cv1s[i]=cv_estimation(X_list,5,lambda_seq2[i],joint=F)
  FNRs[i]=cal_FNR(result,true_omega_list)
  #BICs[i]=BIC_custom(S_list,result,ns)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}
cat('lambdas:',lambda_seq2,'\n')
cat('cv:',cv1s,'\n')
#cat('bic:',BICs,'\n')
cat('fnr:',FNRs,'\n')
cat('fpr:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')

#finer grid search with bic
lambda_seq3=c((5:9)*lambda_seq1[which.min(BICs)]/10,(1:6)*lambda_seq1[which.min(BICs)])
BIC1s=numeric(0)
for (i in 1:length(lambda_seq2)) {
  result=list()
    for(k in (1:length(S_list))){
      result[[k]]=glasso_costm(S_list[[k]],lambda_seq3[i])$WI
     }
  #cv1s[i]=cv_estimation(X_list,5,lambda_seq2[i])
  FNRs[i]=cal_FNR(result,true_omega_list)
  BIC1s[i]=BIC_custom(S_list,result,ns)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}
cat('lambdas:',lambda_seq3,'\n')
#cat('cv:',cv1s,'\n')
cat('bic:',BIC1s,'\n')
cat('fnr:',FNRs,'\n')
cat('fpr:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')
print('end tuning for SINGLE')




#parameter tuning for joint estimation case
print('start tuning for joint')
FPRs=numeric(0)
FNRs=numeric(0)
lambda_seq1=10^(seq(-9,1,1))
cvs=numeric(0)
BICs=numeric(0)
ELs=numeric(0)
FLs=numeric(0)
CZs=numeric(0)
for (i in 1:length(lambda_seq1)) {
  print(i)
  result=cool_glasso(S_list,lambda_seq1[i])
  cvs[i]=cv_estimation(X_list,5,lambda_seq1[i])
  BICs[i]=BIC_custom(S_list,result,ns)
  FNRs[i]=cal_FNR(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}
cat('lambdas:',lambda_seq1,'\n')
cat('cv:',cvs,'\n')
cat('bic:',BICs,'\n')
cat('FNR:',FNRs,'\n')
cat('FPR:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')

#finer grid search with CV
print('start tuning with cv')
lambda_seq2=c((5:9)*lambda_seq1[which.min(cvs)]/10,(1:6)*lambda_seq1[which.min(cvs)])
cv1s=numeric(0)
for (i in 1:length(lambda_seq2)) {
  print(i)
  result=cool_glasso(S_list,lambda_seq2[i])
  cv1s[i]=cv_estimation(X_list,5,lambda_seq2[i])
  FNRs[i]=cal_FNR(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}
cat('lambdas:',lambda_seq2,'\n')
cat('cv:',cv1s,'\n')
cat('lambda selected with cv:',lambda_seq2[which.min(cv1s)],'\n')
cat('FNR:',FNRs,'\n')
cat('FPR:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')

#finer grid search with bic
print('start tuning with bic')
lambda_seq3=c((5:9)*lambda_seq1[which.min(BICs)]/10,(1:6)*lambda_seq1[which.min(BICs)])
BIC1s=numeric(0)
for (i in 1:length(lambda_seq3)) {
  print(i)
  result=cool_glasso(S_list,lambda_seq3[i])
  BIC1s[i]=BIC_custom(S_list,result,ns)
  FNRs[i]=cal_FNR(result,true_omega_list)
  FPRs[i]=cal_FPR(result,true_omega_list)
  ELs[i]=cal_EL(result,true_omega_list)
  FLs[i]=cal_FL(result,true_omega_list)
  CZs[i]=cal_CZ(result,true_omega_list)
}

cat('lambdas:',lambda_seq3,'\n')
cat('bic:',BIC1s,'\n')
cat('lambda selected with bic:',lambda_seq3[which.min(BIC1s)],'\n')
cat('FNR:',FNRs,'\n')
cat('FPR:',FPRs,'\n')
cat('el:',ELs,'\n')
cat('fl:',FLs,'\n')
cat('cz:',CZs,'\n')

##### RUN SIMULATION-----KNN CASE #####

#generate
generate_result=generate_joint_knn(n=sample_num,p=p,rho=rho,m=5,common.seed=111,K=3)
true_omega_list=generate_result$true_omega_list
S_list=generate_result$S_list
X_list=generate_result$X_list
ns=rep(sample_num,length(S_list))

#parameter tuning is exactly the same as last part, simply copy paste would work.





##### REAL DATA ANALYSIS #####


#set the number of variables in p
p=100


data <- readMat("~/proj_gl/train_converted.mat")[[1]]
sparse_matrix <- sparseMatrix(i = data[,1], j = data[,2], x = data[,3])
label=read.table('~/proj_gl/train.label')
maps=read.table('~/proj_gl/train.map')
#select group 12:15
print(maps[12:15,])
index=which(label==12|label==13|label==14|label==15)
index1=which(label==12)
index2=which(label==13)
index3=which(label==14)
index4=which(label==15)
mat=sparse_matrix[index,]
counts=colSums(mat)
tmat=apply((mat),MARGIN=2,function(x){x/sum(x)})
tmat[is.nan(tmat)]=0
cal_entropy=function(x){
  y=ifelse(x==0,0,log(x))
  sum(y*x)/length(x)+1
}
entropy=apply((tmat),MARGIN=2,cal_entropy)
entropy=ifelse(entropy==1,-Inf,entropy)
index_words=order(entropy,decreasing = T)[1:p]
selected_entropy=entropy[index_words]
selected_mat=t(sparse_matrix[,index_words])
X=t(apply(as.matrix(selected_mat), 2, function(x){selected_entropy*log(x+1)}))
X=apply(X, 2, scale)
index=list(index1,index2,index3,index4)

S_list=list()
ns=numeric(0)
X_list=list()
for(i in 1:4){
    X_mat=X[index[[i]],]
    S_list[[i]]=cov(X_mat)
    ns[i]=length(index[[i]])
    X_list[[i]]=X_mat
}


lambda_seq1=10^seq(-8,-2,1)
#cvs=numeric(0)
BICs=numeric(0)

print('start tuning')
for (i in 1:length(lambda_seq1)) {
  print(i)
  result=cool_glasso(S_list,lambda_seq1[i])
#  cvs[i]=cv_estimation(X_list,5,lambda_seq1[i])
  BICs[i]=BIC_custom(S_list,result,ns)
  print(BICs[i])
}
#cat('cv:',cvs,'\n')
cat('bic:',BICs,'\n')

#for bics
print('start tuning with bic')
lambda_seq3=c((5:9)*lambda_seq1[which.min(BICs)]/10,(1:6)*lambda_seq1[which.min(BICs)])
BIC1s=numeric(0)
for (i in 1:length(lambda_seq3)) {
  print(i)
  result=cool_glasso(S_list,lambda_seq3[i])
  BIC1s[i]=BIC_custom(S_list,result,ns)
}

cat('bic:',BIC1s,'\n')
cat('lambda selected with bic:',lambda_seq3[which.min(BIC1s)])

result=cool_glasso(S_list,lambda_seq3[which.min(BIC1s)])

y=lapply(result,function(x){x-diag(diag(x))})
v1=which(y[[1]]!=0)
v2=which(y[[2]]!=0)
v3=which(y[[3]]!=0)
v4=which(y[[4]]!=0)
common_values <- Reduce(intersect, list(v1, v2, v3, v4))
#this gives the index of common nonzeros values detected in matrices
common_values

save(result,file='result.RData')
save(index_words,file='index_words.RData')
