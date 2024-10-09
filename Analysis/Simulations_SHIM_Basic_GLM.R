libs_path<-file.path("..","libs")
source(file.path(libs_path,'Shim3way_GLM.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))



data<- create_basic_dataset_bern()


X<- data$X
y<- data$y$obs
y<-array(y, dim=(c(length(y),1)) )

hist(y)





#res_lasso<-irlasso.cb(X=X, Y=y, lambda=0.001, w.lambda=NULL, beta0=NULL,
#  centering=FALSE, scaling=FALSE, intercept=T,
#  maxit=10, tol=0.0545, sd.tol=1e-6,
#  verbose=TRUE)

#res_lasso$beta

#y_pred<-kappa1(X%*%array(res_lasso$BETA[-1], dim=c(length(res_lasso$BETA),1) )[-1]  + res_lasso$BETA[1]*1 )
#lasso_init<-res_lasso$BETA[1]

beta_true<- data$beta[-1,]
l1=8
l2=8
l3=4
#print(beta_true)

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_matrix<-get_theta_from_theta_vec_2way3(beta_2way,l1=l1,l2=l2, l3=l3)
beta_3way_table<-get_psi_from_psi_vec3(beta_3way,l1=l1,l2=l2, l3=l3)



beta_main_recovered<- get_all_beta(beta_main, l1=l1, l2=l2, l3=l3, threshold = 1e-4)
beta_2way_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-4), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_recovered<- get_psi_vec3( get_all_psi(beta_3way_table, l1=l1, l2=l2, l3=l3, threshold = 1e-4) , l1=l1+1, l2=l2+1, l3=l3+1)


beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0


source(file.path(libs_path,'tools.R'))
rez<-plsglm.cb(X=X, Y=y, ncomp=6, beta0=NULL,
                      centering=FALSE, scaling=FALSE, intercept=TRUE,
                      maxit=10, tol=0.0545,
                      verbose=FALSE, clip=0.01)
coefs_all_pls<-rez$BETA[, ,6]
####START LASSO


###crossval lasso 

cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 1e-2, 5e-2,1e-3, 5e-3, 1e-4, 5e-4, 1e-5), l1=l1,l2=l2,l3=l3, split_percentage = 0.6)

lambda=1e-4 #lambda chose by cross val

res_lasso<-irlasso.cb(X=X, Y=y, lambda=lambda, w.lambda=NULL, beta0=NULL,
                      centering=FALSE, scaling=FALSE, intercept=T,
                      maxit=10, tol=0.0545, sd.tol=1e-6,
                      verbose=TRUE)
res_lasso$beta

coefs_lasso<-array(res_lasso$beta[-1,1,1])
interc_init<-res_lasso$beta[1,1,1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case

beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso_without_gamma*gamma_hat, l1=l1, l2=l2, l3=l3, only_beta = TRUE) #maybe better for shim init
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

predict_lasso<-kappa1(X%*%array(coefs_lasso, dim=c(length(coefs_lasso),1) )  + interc_init  ) #no intercept

print(r2(y, predict_lasso))

plot(y,predict_lasso)


delta_true

#length(c(beta_main, beta_2way, beta_3way))



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



##PREPARE FOR SHIM ###################################################################################################

#beta_hat<-beta_main_lasso
#gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
#gamma_hat[is.nan(gamma_hat)]<-0
#gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
#delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
#delta_hat[!is.finite(delta_hat)]<-0
#delta_hat[is.nan(delta_hat)]<-0



##USE SHIM MODEL #########

lambda_beta<-1e-5
lambda_gamma<-1e-5
lambda_delta<-1e-4



#delta_hat=delta_hat+rnorm(length(delta_hat))

#delta_true

#source(file.path(libs_path,'Shim3way_GLM.R'))
my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(1e-3,1e-4,1e-5), lambda_values_2way=c(1e-4, 1e-5, 1e-3), lambda_delta=1e-4,
                              intercept=interc_init, split_percentage = 0.6)

fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, 
                    w_gamma = 1, w_delta = 1, tol=1e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE, bind_C = FALSE)
#fitted
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





##### USE PIPELINES #####
#pipeline lasso

lambda<-1e-3
results_lasso<-results_pipeline_lasso_GLM(X=X, y=y, lambda=lambda, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
                                      beta_3way_recovered, threshold = 0, strong = TRUE, use_intercept=TRUE)


beta_main_lasso<- results_lasso$main
beta_2way_lasso<-results_lasso$'2way'
beta_3way_lasso<-results_lasso$'3way'
intercept_lasso<-results_lasso$intercept





### pipeline pls
source(file.path(libs_path,'Shim3way_GLM.R'))
results_pls<-results_pipeline_pls_GLM(X=cbind(1,X), y=y, n_comp=6, l1=l1, l2=l2, l3=l3,
                                   beta_main, beta_2way, beta_3way, beta_main_recovered, 
                                   beta_2way_recovered, beta_3way_recovered, threshold = 0, strong = TRUE)

intercept_pls<-results_pls$intercept
beta_main_pls<- results_pls$main
beta_2way_pls<-results_pls$'2way'
beta_3way_pls<-results_pls$'3way'



  


## pipeline shim

#shim with lasso init
lambda_beta<-1e-3
lambda_gamma<-4e-3
lambda_delta<-1e-3

#beta_main_lasso ==beta_hat
results_pipeline_shim_GLM(X=X, y=y, beta_main_init = beta_main_lasso, beta_2way_init = beta_2way_lasso, beta_3way_init = beta_3way_lasso, 
                      lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, l1=l1, l2=l2, l3=l3,
                      beta_main = beta_main, beta_2way = beta_2way, beta_3way = beta_3way , beta_main_recovered = beta_main_recovered, 
                      beta_2way_recovered = beta_2way_recovered, beta_3way_recovered = beta_3way_recovered, tol=1e-2,intercept=intercept_lasso,
                      use_intercept=TRUE, bind_C=FALSE, threshold=0)



#shim with pls init ########## TAKE CARE MODIFIED INTERCEPT
#beta_main_lasso ==beta_hat

coefs_pls<-coefs_all_pls[-1]
intercept_pls<-coefs_all_pls[1]
beta_main_pls<-coefs_pls[range_main]
beta_2way_pls<-coefs_pls[range_theta]
beta_3way_pls<-coefs_pls[range_psi]

all_beta_functions(beta_main, beta_main_pls)
all_beta_functions(beta_2way, beta_2way_pls)
all_beta_functions(beta_3way, beta_3way_pls)



results_pipeline_shim_GLM(X=X, y=y, beta_main_init = beta_main_pls, beta_2way_init = beta_2way_pls, beta_3way_init = beta_3way_pls, 
                          lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, l1=l1, l2=l2, l3=l3,
                          beta_main = beta_main, beta_2way = beta_2way, beta_3way = beta_3way , beta_main_recovered = beta_main_recovered, 
                          beta_2way_recovered = beta_2way_recovered, beta_3way_recovered = beta_3way_recovered, tol=1e-2,intercept=intercept_pls,
                          use_intercept=TRUE, bind_C=FALSE, threshold=0)


