## documents intro
all_codes.R  Contains all codes including GLASSO, Joint GLASSO, simulation, metrics functions and real world data analysis;  
train.label  Contains the category information of each document;  
vocabulary.txt  Contains the words appears in at least one document;  
train_converted.mat  Is a n*p sparse matrix, ij-th element denotes the number of times that j-th word appears in i-th document;  
train.map  Contains a map from the matrix column to the vocabulary.  

## references:
### main references:
Friedman, J., Hastie, T., Tibshirani, R.: Sparse inverse covariance estimation with the graphical lasso. Biostatistics 9(3), 432–441 (12 2007)  
Guo, J., Levina, E., Michailidis, G., Zhu, J.: Joint estimation of multiple graphical models. Biometrika 98(1), 1–15 (02 2011)  
Friedman, J., Hastie, T., H ̈ofling, H., Tibshirani, R.: Pathwise coordinate optimization. The Annals of Applied Atatistics (2007)  


### other references:
Banerjee, O., El Ghaoui, L., d’Aspremont, A.: Model selection through sparse maximum likelihood estimation for multivariate gaussian or binary data. The Journal of Machine Learning Research 9, 485–516 (2008)  
Dumais, S.T.: Improving the retrieval of information from external sources. Behavior research methods,instruments, & computers 23(2), 229–236 (1991)  
Fan, J., Feng, Y., Wu, Y.: Network exploration via the adaptive lasso and scad penalties. The Annals of Applied Statistics 3(2), 521 (2009)  
Lang, K.: Newsweeder: Learning to filter netnews. In: Machine learning proceedings 1995, pp. 331–339. Elsevier (1995)  
Li, H., Gui, J.: Gradient directed regularization for sparse gaussian concentration graphs, with applications to inference of genetic networks. Biostatistics 7(2), 302–317 (2006)  
Meinshausen, N., B ̈uhlmann, P.: High-dimensional graphs and variable selection with the lasso. The Annals of Statistics (2006)


