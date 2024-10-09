####Sparse hierarchical simulations####

library(glmnet)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'Create_synthetic_datasets.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))



data=create_hier_dataset_medium()
X<-data$X
y<-data$y$obs
beta_trues<-data$beta[-1,]  ##without intercept


#5,4,3
l1=5
l2=4
l3=3
#y
#beta_trues

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

colnames(X)[range_main]

beta_true<-beta_trues[range_main]
theta_true<-beta_trues[range_theta]
psi_true<-beta_trues[range_psi]

print(dim(X))
print(dim(y))


dim(X)


X_only_main<-X[,range_main]
X_2way<-X[,range_theta]
X_3way<-X[,range_psi]
y_all<-as.vector(y[,1])






################################### ANALYSIS WITH LASSO ##############################################################

lasso_model <- glmnet(X, y_all, alpha = 1, lambda=0.5)
coefs<-coefficients(lasso_model)[-1]
coefs
sum(coefs==0)

all_beta_functions(beta_true, coefs[range_main])
all_beta_functions(theta_true, coefs[range_theta])
all_beta_functions(psi_true, coefs[range_psi])




#####

beta_init_lasso<-coefs[range_main]
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)
coeffs_main<-coefs[range_main]

#theta_init<-matrix(0, ncol = dim(X_only_main)[2]







########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################


source(file.path(libs_path,'WeakHierNetSeq12_3way.R'))
r_1way<-y_all - mean(y_all)-X[,range_main]%*%coeffs_main

#print(sum(r_1way^2))
#print(sum(r_2way^2))
beta_bound<-2*coeffs_main
lmd<-10
t<-2e-5
theta_init<-get_theta_from_theta_vec_2way3(coefs[range_theta],l1=l1,l2=l2, l3=l3)

myseq12<-WeakHierNet_seq_2way3(X=X_2way, theta_init = theta_init, y=r_1way, beta_bound=beta_bound, lambda=lmd, t=t, tol=1e-7, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3, scale=FALSE)
fitted<-myseq12$fit(X=X_2way, theta_init= theta_init, y=r_1way, beta_bound=beta_bound, lambda=lmd, t=t, tol=1e-7, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3)

all_beta_functions(theta_true, c(fitted$vec_theta_hat))

print(coeffs_main)
print(fitted$theta_hat)

print(length(range_theta))













###M ANALYSIS WITH my library######################################################################
lmd<-30
t<-1e-3

myWeakHierNet<-WeakHierNetUnscaled (X=X_only_main, Beta_plus_init=abs(min(beta_init_lasso))+beta_init_lasso, Beta_minus_init=matrix(abs(min(beta_init_lasso)),nrow = dim(X_only_main)[2], ncol = 1), 
                                    Theta_init=theta_init, y=y_all, lambda=lmd, t=t, tol=1e-9, 
                                    max_iter=10000, eps=1e-8)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X_only_main,Beta_plus_init=abs(min(beta_init_lasso))+beta_init_lasso, Beta_minus_init=matrix(abs(min(beta_init_lasso)),nrow = dim(X_only_main)[2], ncol = 1), 
                         Theta_init=theta_init, y=y_all, lambda=lmd, t=t, tol=1e-9, 
                         max_iter=10000, eps=1e-8,l1=l1,l2=l2,l3=l3 )

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
#print(fitted$vec_theta_hat)
print(fitted$Theta_hat)
sum(fitted$vec_theta_hat==0)
sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus==0)
print(myWeakHierNet$R2_score(fitted,X_only_main,y_all))



###RESULTS FOR BETA AND THETA#####
all_beta_functions(beta_true, fitted$Beta_hat_plus-fitted$Beta_hat_minus)
all_beta_functions(theta_true, fitted$vec_theta_hat)



###HIERNET LIBRARY
#print("-----Hiernet library-----")

fit=hierNet(X_only_main, y_all, lam=10, diagonal = FALSE, stand.main=FALSE,tol=1e-7)
fit$th
predicted=predict(fit,X_only_main)
print(r2(y_all, predicted))

#fit$bp-fit$bn
#print(sum( get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)==0))

#all_beta_functions(beta_true, c(fit$bp-fit$bn) )
#all_beta_functions(theta_true, get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3))


#######

#Residuals
predicted<-myWeakHierNet$predict(fitted, X_only_main)
r_2way<-y_all-predicted


##Fit in on residuals

lasso_model <- glmnet(X_3way, r_2way, alpha = 1, lambda=0.05)
coefficients(lasso_model)

psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model)[-1]),l1=l1,l2=l2,l3=l3) ## psi, get intercept out
print(dim(psi_init))
print(dim(X_3way))
theta_bound<-fitted$Theta_hat +t(fitted$Theta_hat)*10 ##bound
lambda<-40


myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
                                    l1=l1,l2=l2,l3=l3, scale = FALSE)

# Fit the model
fitted=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
                            eps=1e-8,l1=l1,l2=l2,l3=l3)

r2(r_2way,myWeakHierNet_seq3$predict(fitted,X_3way))
#print(fitted$vec_psi_hat)

sum(fitted$vec_psi_hat==0)

###RESULTS FOR PSI#####
all_beta_functions(psi_true, c(fitted$vec_psi_hat))
#######





########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################


