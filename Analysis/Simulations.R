####Sparse hierarchical simulations####

library(glmnet)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'Create_synthetic_datasets.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))



data=create_sparse_hier_dataset()
X<-data$X
y<-data$y
beta_trues<-data$beta



beta_true<-beta_trues[1:16]
theta_true<-beta_trues[17:91]
psi_true<-beta_trues[92:199]

print(dim(X))
print(dim(y))

l1=4
l2=9
l3=3

colnames(X)[118]
range_main<-c(1:16)
range_theta<-c(17:91)
range_psi<-c(92:199)

X_only_main<-X[,c(1:16)]
X_3way<-X[,c(92:199)]
X_2way<-X[,range_theta]
y_all<-as.vector(y[,1])
#y_all


###HIERNET LIBRARY
print("-----Hiernet library-----")
##On train
fit=hierNet(X_only_main, y_all, lam=10, diagonal = FALSE, stand.main=FALSE,tol=1e-10)
#fit
#y_pred_train=predict(fit,X_only_main)

get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)

###RESULTS FOR BETA AND THETA#####
all_beta_functions(beta_true, c(fit$bp-fit$bn) )
all_beta_functions(theta_true, get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3))






#####
####TRUTH
#beta[1]=1 #A1
#beta[4]=4 #A4
#beta[9]=9 #B5
#beta[13]=-13 #B9
#beta[16]=-1 #C3

#2 way
#beta[21]=2 #A1:B5 i.e.5 
#beta[25]=-2 #A1:B9 i.e. 9
#beta[91]=9 #B9:C3  i.e. ultimul theta

#3way
#beta[106]=2  #A1B5C3  sau 106-91=15
#beta[105]=1  #A1B5C2             =14 
#beta[118]=-3 #A1B9C3     =25



#####
lmd<-7
t<-5e-3

myWeakHierNet<-WeakHierNetUnscaled (X=X_only_main, Beta_plus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), 
                                    Theta_init=matrix(0, ncol = dim(X_only_main)[2], nrow = dim(X_only_main)[2]), y=y_all, lambda=lmd, t=t, tol=1e-6, 
                                    max_iter=10000, eps=1e-8)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X_only_main, Beta_plus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), 
                         Theta_init=matrix(0, ncol = dim(X_only_main)[2], nrow = dim(X_only_main)[2]), y=y_all, lambda=lmd, t=t, tol=1e-6, 
                         max_iter=10000, eps=1e-8,l1=4,l2=9,l3=3 )

#print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
#print(fitted$vec_theta_hat)
#print(fitted$Theta_hat)
#sum(fitted$vec_theta_hat==0)
#sum(fitted$Beta_hat_minus-fitted$Beta_hat_minus==0)
#print(myWeakHierNet$R2_score(fitted,X_only_main,y_all))


###RESULTS FOR BETA AND THETA#####
all_beta_functions(beta_true, c(fitted$Beta_hat_plus-fitted$Beta_hat_minus))
all_beta_functions(theta_true, c(fitted$vec_theta_hat))

#######

#Residuals
predicted<-myWeakHierNet$predict(fitted, X_only_main)
r_2way<-y_all-predicted


##Fit in on residuals

lasso_model <- glmnet(X_3way, r_2way, alpha = 1, lambda=0.05)
coefficients(lasso_model)

psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model)[-1]),l1=4,l2=9,l3=3) ## psi, get intercept out
theta_bound<-fitted$Theta_hat +t(fitted$Theta_hat)*10 ##bound
lambda<-10


myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
                                    l1=4,l2=9,l3=3, scale = FALSE)

# Fit the model
fitted=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
                            eps=1e-8,l1=4,l2=9,l3=3)

r2(r_2way,myWeakHierNet_seq3$predict(fitted,X_3way))
print(fitted$vec_psi_hat)

sum(fitted$vec_psi_hat==0)

###RESULTS FOR PSI#####
all_beta_functions(psi_true, c(fitted$vec_psi_hat))
#######


################################### ANALYSIS WITH LASSO ##############################################################

lasso_model <- glmnet(X, y_all, alpha = 1, lambda=0.02)
coefs<-coefficients(lasso_model)[-1]

coefs

sum(coefs==0)

all_beta_functions(beta_true, coefs[range_main])
all_beta_functions(theta_true, coefs[range_theta])
all_beta_functions(psi_true, coefs[range_psi])

coeffs_main<-coefs[range_main]

print(colnames(X))

########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################


source(file.path(libs_path,'WeakHierNetSeq12_3way.R'))
r_1way<-y_all - mean(y_all)-X[,range_main]%*%coeffs_main

#print(sum(r_1way^2))
#print(sum(r_2way^2))
beta_bound<-2*coeffs_main
lmd<-10
t<-4e-4
theta_init<-get_theta_from_theta_vec_2way3(coefs[range_theta],l1=l1,l2=l2, l3=l3)

myseq12<-WeakHierNet_seq_2way3(X=X_2way, theta_init = theta_init, y=r_1way, beta_bound=beta_bound, lambda=lmd, t=t, tol=1e-8, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3, scale=FALSE)
fitted<-myseq12$fit(X=X_2way, theta_init= theta_init, y=r_1way, beta_bound=beta_bound, lambda=lmd, t=t, tol=1e-8, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3)

all_beta_functions(theta_true, c(fitted$vec_theta_hat))


print(length(range_theta))



