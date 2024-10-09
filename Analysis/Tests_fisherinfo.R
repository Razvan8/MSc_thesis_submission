
libs_path<-file.path("..","libs")
source(file.path(libs_path,'Shim3way_class_sc_loss.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))



data<- create_basic_dataset()
X<- data$X
y<- data$y$obs
mean(y)

beta_true<- data$beta[-1,]
l1=8
l2=7
l3=6
#print(beta_true)

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_matrix<-get_theta_from_theta_vec_2way3(beta_2way,l1=l1,l2=l2, l3=l3)
beta_3way_table<-get_psi_from_psi_vec3(beta_3way,l1=l1,l2=l2, l3=l3)



beta_main_recovered<- get_all_beta(beta_main, l1=l1, l2=l2, l3=l3, threshold = 0)
beta_2way_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_matrix, l1=l1, l2=l2, l3=l3, threshold = 0), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_recovered<- get_psi_vec3( get_all_psi(beta_3way_table, l1=l1, l2=l2, l3=l3, threshold = 0) , l1=l1+1, l2=l2+1, l3=l3+1)


beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

y_centered<-y-mean(y)


####START LASSO ##MAKE IT CLEAN FOR BOTH CASES##########



cv_fit <- cv.glmnet(X, y_centered, alpha = 1, intercept= FALSE, standardize = FALSE)
best_lambda <- cv_fit$lambda.min
cat("best lambda: ",best_lambda)
best_lambda<-0.06## cv best lambda varies a bit at different runs. SET TO 0.06 for cv or 0.16 for TPR90%!!!
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=best_lambda)

coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)

predict_lasso<-predict(lasso_model, s=best_lambda, newx = X)
print(r2(y_centered, predict_lasso))

length(c(beta_main, beta_2way, beta_3way))



######################### RESULTS LASSO #########################################################
print(" results lasso without recovered")
all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)





##### RESULTS LASSO ON RECOVERED PARAMS ##########
beta_2way_lasso_matrix<-get_theta_from_theta_vec_2way3(beta_2way_lasso,l1=l1,l2=l2, l3=l3)
beta_3way_lasso_table<-get_psi_from_psi_vec3(beta_3way_lasso,l1=l1,l2=l2, l3=l3)

##hierarchy tests

test_hierarchy_layer12(beta_main_lasso,beta_2way_lasso_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_lasso_matrix, beta_3way_lasso_table, strong = TRUE)

beta_main_lasso_recovered<- get_all_beta(beta_main_lasso, l1=l1, l2=l2, l3=l3, threshold = 0)
beta_2way_lasso_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_lasso_matrix, l1=l1, l2=l2, l3=l3, threshold = 0), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_lasso_recovered<- get_psi_vec3( get_all_psi(beta_3way_lasso_table, l1=l1, l2=l2, l3=l3, threshold = 0) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results lasso recovered")
all_beta_functions(beta_main_recovered, beta_main_lasso_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_lasso_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_lasso_recovered)

test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)










##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL #########

lambda_beta<- 0.12#0.2 #0.12
lambda_gamma<- 0.6  #0.8 #0.6
lambda_delta<-1  #9 #1

my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)

#Initial large search
#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.01,0.1, 0.5), lambda_values_2way=c( 0.1, 0.8,1.5), lambda_delta=1, split_percentage = 0.6, k=3)

#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.15, 0.2, 0.4), lambda_values_2way=c( 0.6,0.8), lambda_delta=1, split_percentage = 0.6, k=3)


#Final search
#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.5, 0.1, 0.05, 0.01), lambda_values_2way=c( 1,0.5,0.1,0.05), lambda_delta=1, split_percentage = 0.6, k=3)

#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.05,0.1,0.5,1), lambda_values_2way=c( 1,0.8,0.6), lambda_delta=1, split_percentage = 0.6, k=3)


#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.12,0.08), lambda_values_2way=c( 0.8), lambda_delta=1, split_percentage = 0.6, k=3)

#cv_good<-cross_val_pipeline( X=X, y=y, lambda=best_lambda, lambda_values_main=c( 0.1,0.15, 0.2, 0.25 ), 
                             #lambda_values_2way=c(0.1,0.3,0.6,0.9), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=3, l1=l1, l2=l2, l3=l3)

cv_good<-cross_val_pipeline( X=X, y=y, lambda=best_lambda, lambda_values_main=c( 0.08, 0.1 ), 
                             lambda_values_2way=c(0.8), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=3, l1=l1, l2=l2, l3=l3)


lambda_beta<- cv_good$best_lambda1#0.08
lambda_gamma<- cv_good$best_lambda2   #0.8
lambda_delta<-1   #1
fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=1e-2)



my_shim$R2_score(self=fitted, X_new=X, y_true=y )




beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]



all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)

##hierarchy tests
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_shim,beta_2way_shim_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = TRUE)



### SHIM FOR RECOVERED DATA #############
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(beta_main_shim, l1=l1, l2=l2, l3=l3, threshold = 0)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 0), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 0) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results shim recovered")
all_beta_functions(beta_main_recovered, beta_main_shim_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_shim_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_shim_recovered)

sum(beta_main_recovered!=0)
sum(beta_2way_recovered!=0)
sum(beta_3way_recovered!=0)

test_hierarchy_layer12(beta_main_shim_recovered, get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)














####### SAME ANALYSIS BUT DIFFERENT LAMBDA (smaller - get aroung 80% of zeros) ###################
####START LASSO



cv_fit <- cv.glmnet(X, y_centered, alpha = 1)

best_lambda<-0.08
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=best_lambda)

coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)

predict_lasso<-predict(lasso_model, s=best_lambda, newx = X)
print(r2(y_centered, predict_lasso))

length(c(beta_main, beta_2way, beta_3way))

######################### RESULTS LASSO #########################################################
print(" results lasso without recovered")
all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)


##hierarchy tests

beta_2way_lasso_matrix<-get_theta_from_theta_vec_2way3(beta_2way_lasso,l1=l1,l2=l2, l3=l3)
beta_3way_lasso_table<-get_psi_from_psi_vec3(beta_3way_lasso,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_lasso,beta_2way_lasso_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_lasso_matrix, beta_3way_lasso_table, strong = TRUE)


##### RESULTS LASSO ON RECOVERED PARAMS ##########

beta_main_lasso_recovered<- get_all_beta(beta_main_lasso, l1=l1, l2=l2, l3=l3, threshold = 1)
beta_2way_lasso_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_lasso_matrix, l1=l1, l2=l2, l3=l3, threshold = 1), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_lasso_recovered<- get_psi_vec3( get_all_psi(beta_3way_lasso_table, l1=l1, l2=l2, l3=l3, threshold = 1) , l1=l1+1, l2=l2+1, l3=l3+1)


#beta_main_lasso
#beta_main_lasso_recovered

print("results lasso recovered")
all_beta_functions(beta_main_recovered, beta_main_lasso_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_lasso_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_lasso_recovered)


test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)



##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL #########

lambda_beta<-50
lambda_gamma<-900
lambda_delta<-1000


my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=1e-3)

my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]



all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)

##hierarchy tests
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_shim,beta_2way_shim_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = TRUE)



### SHIM FOR RECOVERED DATA #############
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(beta_main_shim, l1=l1, l2=l2, l3=l3, threshold = 1)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 1), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 1) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results shim recovered")
all_beta_functions(beta_main_recovered, beta_main_shim_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_shim_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_shim_recovered)

test_hierarchy_layer12(beta_main_shim_recovered, get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)






############# USE PIPELINES


##### results pipeline Lasso ######
lambda<-0.16 ##0.27 is best lambda
results_lasso<-results_pipeline_lasso(X=X, y=y_centered, lambda=lambda, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
                       beta_3way_recovered, threshold = 0, strong = TRUE)


beta_main_lasso<- results_lasso$main
beta_2way_lasso<-results_lasso$'2way'
beta_3way_lasso<-results_lasso$'3way'



############# PLS ############
n_comp=6
results_pls<-results_pipeline_pls(X=X, y=y_centered, n_comp=n_comp, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
                                  beta_3way_recovered, threshold = 0, strong = TRUE)
beta_main_pls<- results_pls$main
beta_2way_pls<-results_pls$'2way'
beta_3way_pls<-results_pls$'3way'






lambda_beta<-200
lambda_gamma<-2000
lambda_delta<-500


## results pipeline shim



#shim with lasso init
results_pipeline_shim(X=X, y=y, beta_main_init = beta_main_lasso, beta_2way_init = beta_2way_lasso, beta_3way_init = beta_3way_lasso, 
                      lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, l1=l1, l2=l2, l3=l3,
                      beta_main = beta_main, beta_2way = beta_2way, beta_3way = beta_3way , beta_main_recovered = beta_main_recovered, 
                      beta_2way_recovered = beta_2way_recovered, beta_3way_recovered = beta_3way_recovered, tol=1e-2, threshold = 0)




#shim with pls init


my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_main_pls, gamma_init = beta_2way_pls, delta_init = beta_3way_pls, l1=l1, l2=l2, l3=l3, scale = FALSE)
cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(200, 300,500,1000), lambda_values_2way=c( 200,500,1000,2500), lambda_delta=3000, split_percentage = 0.6)


lambda_beta<-200
lambda_gamma<-2000
lambda_delta<-500

results_pipeline_shim(X=X, y=y, beta_main_init = beta_main_pls, beta_2way_init = beta_2way_pls, beta_3way_init = beta_3way_pls, 
                      lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, l1=l1, l2=l2, l3=l3,
                      beta_main = beta_main, beta_2way = beta_2way, beta_3way = beta_3way , beta_main_recovered = beta_main_recovered, 
                      beta_2way_recovered = beta_2way_recovered, beta_3way_recovered = beta_3way_recovered, tol=1e-3, threshold = 0)




abline(0,1)



















