libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
source(file.path(libs_path,'Shim3way_GLM.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))


library(glmnet)

data<-load_deoxi_flourination()

print(max(data$yield))
hist(data$prob)

which(data$yield_corrected==100)
which(data$yield>100)




plot_yield_histo(data)


levels(data$a)
levels(data$b)
levels(data$s)

y_true<-data$prob*100
X_data<-data[c('a','b','s')]
colnames(X_data)<-c("alcohol", 'base', "sulfonyl-fluoride")


X_data <- lapply(X_data, as.factor)
X_data <- as.data.frame(X_data)  # Convert back to data frame


# Create the dummy matrix without intercept
X_dummy <- model.matrix(~ . - 1, data = X_data, contrasts.arg = lapply(X_data, contrasts, contrasts = FALSE))

# Find rows where the 37th column has value 1
rows_to_replace1 <- X_dummy[,37] == 1
# Replace values in columns 1 to 36 with -1 for these rows
X_dummy[rows_to_replace1, 1:36] <- -1


rows_to_replace2 <- X_dummy[,41] == 1
X_dummy[rows_to_replace2, 38:40] <- -1


rows_to_replace3 <- X_dummy[,46] == 1
X_dummy[rows_to_replace3, 42:45] <- -1


print(dim(X_dummy))
print(colSums(X_dummy))

##only once
X_dummy<-X_dummy[,-c(37,41,46)]
print(dim(X_dummy))
print(colSums(X_dummy))
###

new_colnames <- c(paste0("alcohol ", 1:36), paste0("base ", 1:3), paste0("sulfonyl-fluoride", 1:4))
colnames(X_dummy) <- new_colnames

X_dummy<-as.data.frame(X_dummy, col.names =TRUE)


for (i in 1:36) {
  for (jk in 37:43) {
    
    X_dummy[, paste0(colnames(X_dummy)[i], ':', colnames(X_dummy)[jk])] <-X_dummy[,i]*X_dummy[,jk]
  }
}

for (j in c(37:39))
{for (k in c(40:43))
{
  X_dummy[,paste0((colnames(X_dummy)[j]),':', colnames(X_dummy[k]) ) ] <- X_dummy[,j]*X_dummy[,k]
}
}




for (i in c(1:36))
{for (j in c(37:39))
{for (k in c(40:43))
{
  X_dummy[,paste0( (colnames(X_dummy)[i]), ':' , colnames(X_dummy[j]), ":",  colnames(X_dummy[k]) ) ] <- X_dummy[,i]*X_dummy[,j]*X_dummy[,k]
}
}}

X<-as.matrix(X_dummy)
y<-y_true


################################################################ USE LASSO AND MY CLASS ON DATA #################################


l1=36
l2=3
l3=4
#print(beta_true)

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

#y_centered<-y-mean(y)
y<-y/100 ############# only once

####START LASSO


#cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 1e-3, 6e-4, 3e-4, 1e-4, 6e-5), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)
cv<-cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 1e-3, 6e-4, 3e-4, 1e-4, 6e-5), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 5e-4), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

##CV FOR PLOT
cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 3e-4), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

###correct crossval for SHIM

lambda<-cv$best_lambda
lambda<-3e-4#what cv choses (keep this to run faster)

#initial cv
#cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c( 1e-5, 1e-4,1e-3), 
                             #lambda_values_2way=c( 1e-7,1e-6,1e-5,1e-4,1e-3,1e-2), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c( 1e-4, 5e-4, 1e-3), 
                             lambda_values_2way=c(1e-6), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

#cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c(  5e-4), 
                             #lambda_values_2way=c(5e-4, 1e-4, 1e-6), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

###JUST IF PLOT IS NEEDED
cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c(  5e-4), 
                             lambda_values_2way=c( 1e-5), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)


lambda<-cv$best_lambda
lambda=3e-4 #this is the value of the cross val
res_lasso<-irlasso.cb(X=X, Y=y, lambda=lambda, w.lambda=NULL, beta0=NULL,
                      centering=FALSE, scaling=FALSE, intercept=T,
                      maxit=10, tol=0.0545, sd.tol=1e-6,
                      verbose=TRUE)


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
plot(predict_lasso, y,  xlab = "Predicted Yield", ylab = "True Yield", main = "Predicted vs True Yield  ")
abline(a = 0, b = 1, col = "red")



sum(beta_main_lasso==0)/length(beta_main_lasso)
sum(beta_2way_lasso==0)/length(beta_2way_lasso)
sum(beta_3way_lasso==0)/length(beta_3way_lasso)

sum(beta_main_lasso!=0)
sum(beta_2way_lasso!=0)
sum(beta_3way_lasso!=0)


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
test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)


sum(beta_main_lasso_recovered==0)/length(beta_main_lasso_recovered)
sum(beta_2way_lasso_recovered==0)/length(beta_2way_lasso_recovered)
sum(beta_3way_lasso_recovered==0)/length(beta_3way_lasso_recovered)

sum(beta_main_lasso_recovered!=0)
sum(beta_2way_lasso_recovered!=0)
sum(beta_3way_lasso_recovered!=0)



########################### SHIM ##################################################



##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0



##USE SHIM MODEL #########

lambda_beta<-cv_good$best_lambda1
lambda_gamma<-cv_good$best_lambda2

##if cv_good was not run
lambda_beta<-5e-4
lambda_gamma<-1e-5
lambda_delta<-1e-4#have to run this



#delta_hat=delta_hat+rnorm(length(delta_hat))

#delta_true

source(file.path(libs_path,'Shim3way_GLM.R'))
my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)


fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, 
                    w_gamma = 1, w_delta = 1, tol=1e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE)
#fitted
my_shim$R2_score(self=fitted, X_new=X, y_true=y )

sum(abs(fitted$beta_all))
sum(abs(coefs_lasso))

sum(fitted$beta_all[range_main]==0)/length(range_main)
sum(fitted$beta_all[range_theta]==0)/length(range_theta)
sum(fitted$beta_all[range_psi]==0)/length(range_psi)

sum(fitted$beta_all[range_main]!=0)
sum(fitted$beta_all[range_theta]!=0)
sum(fitted$beta_all[range_psi]!=0)


beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(fitted$beta_all[range_theta],l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(fitted$beta_all[range_psi],l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(fitted$beta_all[range_main], l1=l1, l2=l2, l3=l3, threshold = 0)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 0), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 0) , l1=l1+1, l2=l2+1, l3=l3+1)



sum(beta_main_shim_recovered==0)/length(beta_main_shim_recovered)
sum(beta_2way_shim_recovered==0)/length(beta_2way_shim_recovered)
sum(beta_3way_shim_recovered==0)/length(beta_3way_shim_recovered)

sum(beta_main_shim_recovered!=0)
sum(beta_2way_shim_recovered!=0)
sum(beta_3way_shim_recovered!=0)



order(abs(fitted$beta_all), decreasing=TRUE)
order(abs(coefs_lasso), decreasing = TRUE)


min(sign(fitted$beta_all*coefs_lasso))














##################################################################################################


##same thing but higher lmd

####################################################################################################



#cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 1e-3, 6e-4, 3e-4, 1e-4, 6e-5), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)
cv<-cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 1e-3, 6e-4, 3e-4, 1e-4, 6e-5), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 5e-4), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

##CV FOR PLOT
cross_validation_irlasso.cb( X=X, y=y, lambda_values=c( 8e-4), l1=l1,l2=l2,l3=l3, split_percentage = 0.6, k=4)

###correct crossval for SHIM


lambda<-8e-4#what cv choses (keep this to run faster)

#initial cv
#cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c( 1e-5, 1e-4,1e-3), 
#lambda_values_2way=c( 1e-7,1e-6,1e-5,1e-4,1e-3,1e-2), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c( 1e-4, 5e-4, 1e-3), 
                             lambda_values_2way=c(1e-6), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

#cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c(  5e-4), 
#lambda_values_2way=c(5e-4, 1e-4, 1e-6), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

###JUST IF PLOT IS NEEDED
cv_good<-cross_val_pipeline( X=X, y=y, lambda=lambda, lambda_values_main=c(  1e-3), 
                             lambda_values_2way=c( 1e-4), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)


lambda<-cv$best_lambda
lambda=8e-4 #this is the value of the cross val
res_lasso<-irlasso.cb(X=X, Y=y, lambda=lambda, w.lambda=NULL, beta0=NULL,
                      centering=FALSE, scaling=FALSE, intercept=T,
                      maxit=10, tol=0.0545, sd.tol=1e-6,
                      verbose=TRUE)


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
plot(predict_lasso, y,  xlab = "Predicted Yield", ylab = "True Yield", main = "Predicted vs True Yield  ")
abline(a = 0, b = 1, col = "red")



sum(beta_main_lasso==0)/length(beta_main_lasso)
sum(beta_2way_lasso==0)/length(beta_2way_lasso)
sum(beta_3way_lasso==0)/length(beta_3way_lasso)

sum(beta_main_lasso!=0)
sum(beta_2way_lasso!=0)
sum(beta_3way_lasso!=0)


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
test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)


sum(beta_main_lasso_recovered==0)/length(beta_main_lasso_recovered)
sum(beta_2way_lasso_recovered==0)/length(beta_2way_lasso_recovered)
sum(beta_3way_lasso_recovered==0)/length(beta_3way_lasso_recovered)

sum(beta_main_lasso_recovered!=0)
sum(beta_2way_lasso_recovered!=0)
sum(beta_3way_lasso_recovered!=0)



########################### SHIM ##################################################



##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0



##USE SHIM MODEL #########

lambda_beta<-cv_good$best_lambda1
lambda_gamma<-cv_good$best_lambda2

##if cv_good was not run
lambda_beta<-1e-3
lambda_gamma<-1e-4
lambda_delta<-1e-4#have to run this



#delta_hat=delta_hat+rnorm(length(delta_hat))

#delta_true

source(file.path(libs_path,'Shim3way_GLM.R'))
my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)


fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, 
                    w_gamma = 1, w_delta = 1, tol=1e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE)
#fitted
my_shim$R2_score(self=fitted, X_new=X, y_true=y )

sum(abs(fitted$beta_all))
sum(abs(coefs_lasso))

sum(fitted$beta_all[range_main]==0)/length(range_main)
sum(fitted$beta_all[range_theta]==0)/length(range_theta)
sum(fitted$beta_all[range_psi]==0)/length(range_psi)

sum(fitted$beta_all[range_main]!=0)
sum(fitted$beta_all[range_theta]!=0)
sum(fitted$beta_all[range_psi]!=0)


beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(fitted$beta_all[range_theta],l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(fitted$beta_all[range_psi],l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(fitted$beta_all[range_main], l1=l1, l2=l2, l3=l3, threshold = 0)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 0), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 0) , l1=l1+1, l2=l2+1, l3=l3+1)



sum(beta_main_shim_recovered==0)/length(beta_main_shim_recovered)
sum(beta_2way_shim_recovered==0)/length(beta_2way_shim_recovered)
sum(beta_3way_shim_recovered==0)/length(beta_3way_shim_recovered)

sum(beta_main_shim_recovered!=0)
sum(beta_2way_shim_recovered!=0)
sum(beta_3way_shim_recovered!=0)



order(abs(fitted$beta_all), decreasing=TRUE)
order(abs(coefs_lasso), decreasing = TRUE)


min(sign(fitted$beta_all*coefs_lasso))



















