libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
source(file.path(libs_path,'Shim3way_class_sc_loss.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))


library(glmnet)

data<-load_deoxi_flourination()

print(max(data$yield))
hist(data$prob)

which(data$yield_corrected==100)
which(data$yield>100)

data


plot_yield_histo(data)


levels(data$a)
levels(data$b)
levels(data$s)

y_true<-data$prob*100
X_data<-data[c('a','b','s')]

X_data <- lapply(X_data, as.factor)
X_data <- as.data.frame(X_data)  # Convert back to data frame

unique(X_data$a)
unique(X_data$b)
unique(X_data$s)


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

new_colnames <- c(paste0("a.", 1:36), paste0("b.", 1:3), paste0("s.", 1:4))
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

y_centered<-y-mean(y)

class(y_centered)
####START LASSO

cv_fit <- cv.glmnet(X, y_centered, alpha = 1, intercept =FALSE, standardize=FALSE)

best_lambda <- cv_fit$lambda.min
cat("best lambda: ",best_lambda)
#best_lambda<-0.03
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=best_lambda)
y_pred<-predict(lasso_model, newx=X)
plot((y_pred+mean(y))/100, (y_centered+mean(y))/100, xlab = "Predicted Yield", ylab = "True Yield", main = "Predicted vs True Yield")
abline(a = 0, b = 1, col = "red")
r2(y_centered, y_pred)

max(y_pred+mean(y))
mean((y_centered-y_pred)^2)
#init for shim
coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)



print("Percentages of 0s for main 2way and 3way LASSO")
sum(beta_main_lasso==0)/length(beta_main_lasso)
sum(beta_2way_lasso==0)/ length(beta_2way_lasso)
sum(beta_3way_lasso==0)/ length(beta_3way_lasso)

predict_lasso<-predict(lasso_model, s=best_lambda, newx = X)
print(r2(y_centered, predict_lasso))

########################### SHIM ##################################################



beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0





my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
#Initial large search
cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.01,0.1), lambda_values_2way=c( 0.1, 1), lambda_delta=1, split_percentage = 0.6, k=3)

#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c(0.15, 0.2, 0.4), lambda_values_2way=c( 0.6,0.8), lambda_delta=1, split_percentage = 0.6, k=3)


#Final search
#cv<-my_shim$cross_validation( X=X, y=y, lambda_values_main=c( 0.01, 0.001), lambda_values_2way=c( 0.001,0.1), lambda_delta=1, split_percentage = 0.6, k=3)


####GOOD CROSS VAL PIPELINE##########################

cv_good<-cross_val_pipeline( X=X, y=y, lambda=best_lambda, lambda_values_main=c( 1e-4, 5e-4, 1e-3), 
                             lambda_values_2way=c(1e-6), lambda_delta=1e-4, split_percentage = 0.6, verbose=TRUE, k=4, l1=l1, l2=l2, l3=l3)

########################################################

lambda_beta<-cv$best_lambda1
lambda_gamma<-cv$best_lambda2
lambda_delta<-1e-4



fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=1e-2)

###PLOT FOR THESIS
y_scaled<-y/100
y_pred<-my_shim$predict(self=fitted, X_new=X, mean_y = mean(y))/100
plot(y_pred, y_scaled, xlab = "Predicted yield", ylab = "True yield", main = "Predicted vs true yield")
abline(a = 0, b = 1, col = "red")
####

my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]

print("Percentages of 0s for main 2way and 3way SHIM")
sum(beta_main_shim==0)/length(beta_main_shim)
sum(beta_2way_shim==0)/ length(beta_2way_shim)
sum(beta_3way_shim==0)/ length(beta_3way_shim)

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



beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1)
test_hierarchy_layer12(beta_main_shim_recovered,beta_2way_shim_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = TRUE)







