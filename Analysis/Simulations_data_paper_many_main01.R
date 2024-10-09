####Sparse hierarchical simulations####


library(glmnet)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'generate_dummy01.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'helper_functions.R'))



l1=6
l2=5
l3=4


data=generate_data01(nlevels=c(l1,l2,l3), error_sd=0, intercept=0)
X<-as.matrix(data$X)
y<-data$y
beta_trues<-data$beta  ##without intercept

mean(y)







#y
#beta_trues

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

colnames(X)[range_main]

beta_true<-beta_trues[range_main]
theta_true<-beta_trues[range_theta]
psi_true<-beta_trues[range_psi]



zeros_beta_true<-sum(beta_true==0)
zeros_theta_true<-sum(theta_true==0)
zeros_psi_true<-sum(psi_true==0)
total_zeros<-zeros_beta_true+zeros_theta_true+zeros_psi_true

print(zeros_beta_true)
print(zeros_theta_true)
print(zeros_psi_true)
print(total_zeros)
print(length(beta_trues))

dim(X)
# Perform singular value decomposition
svd_decomposition <- svd(X)

# Extract singular values
singular_values <- svd_decomposition$d

# Calculate the condition number
condition_number <- max(singular_values) / min(singular_values)

# Print the condition number
print(condition_number)


X_only_main<-X[,range_main]

X_2way<-X[,range_theta]
X_3way<-X[,range_psi]
y_all<-as.vector(y[,1])






################################### ANALYSIS WITH LASSO ##############################################################

####Lasso only on main effects
#lasso_model <- glmnet(X[,range_main], y_all, alpha = 1, lambda=0.3)
#coefs<-coefficients(lasso_model)[-1]
#coefs
#sum(coefs==0)

#all_beta_functions(beta_true, coefs[range_main])



lmd<-cv.glmnet(X, y_all, alpha = 1,standardize=FALSE)$lambda.min
lmd
lasso_model <- glmnet(X, y_all, alpha = 1, lambda=lmd, standardize = FALSE)
coefs<-coefficients(lasso_model)[-1]

#coefs
#sum(coefs==0)

#print(coefs[range_theta])

all_beta_functions(beta_true, coefs[range_main])
all_beta_functions(theta_true, coefs[range_theta])
all_beta_functions(psi_true, coefs[range_psi])


#coefs[range_main]
#coefs[range_theta]
#beta_true

test_hierarchy_layer12(coefs[range_main], get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3) )
test_hierarchy_layer23( get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3), get_psi_from_psi_vec3(coefs[range_psi],l1=l1,l2=l2,l3=l3) )


sum(coefs[c(range_main,range_theta)]==0)



r2(y_all, X%*%coefs + coefficients(lasso_model)[1])
r2(y_all, X%*%beta_trues )

print(mean(y_all))


#####

beta_init_lasso<-coefs[range_main]
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)

print(theta_init)
#theta_init<-matrix(0, ncol = dim(X_only_main)[2]




#######ANALYSIS SEQ1-2 + SEQ 2-3 #######################




#####

beta_init_lasso<-coefs[range_main]
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)
coeffs_main<-coefs[range_main]

#theta_init<-matrix(0, nrow=dim(X_only_main)[2], ncol = dim(X_only_main)[2])






###### ANALYSIS WITH my library   #########################################

lmd<-5
t<-1e-2

beta_init_lasso_plus<- beta_init_lasso
beta_init_lasso_plus[beta_init_lasso_plus<0]<-0

beta_init_lasso_minus<- -beta_init_lasso
beta_init_lasso_plus[beta_init_lasso_plus<0]<-0






source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
myWeakHierNet<-WeakHierNetUnscaled (X=X_only_main, Beta_plus_init=beta_init_lasso_plus*0, Beta_minus_init=beta_init_lasso_minus*0, 
                                    Theta_init=theta_init*0, y=y_all, lambda=lmd, t=t, tol=1e-7, 
                                    max_iter=10000, eps=1e-8, center = FALSE, standard_form=FALSE)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X_only_main,Beta_plus_init=beta_init_lasso_plus*0, Beta_minus_init = beta_init_lasso_minus*0, 
                         Theta_init=theta_init*0, y=y_all, lambda=lmd, t=t, tol=1e-7, 
                         max_iter=10000, eps=1e-8,l1=l1,l2=l2,l3=l3 )

r2(y_all,myWeakHierNet$predict(fitted, X_only_main))

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
beta_true


print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)

print(fitted$Beta_hat_plus)
print(fitted$Beta_hat_minus)

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
#print(fitted$vec_theta_hat)
print(fitted$Theta_hat)
print(get_vec_theta_hat3(fit$th,l1=l1,l2=l2,l3=l3))

#sum(fitted$vec_theta_hat==0)
#sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus==0)
#print(myWeakHierNet$R2_score(fitted,X_only_main,y_all))

print(fitted$vec_theta_hat)
theta_true
print(beta_trues[range_theta])
print(get_vec_theta_hat3(fit$th,l1=l1,l2=l2,l3=l3))

print(fit$th)

###RESULTS FOR BETA AND THETA#####
print("My library")
all_beta_functions(beta_true, fitted$Beta_hat_plus-fitted$Beta_hat_minus)
all_beta_functions(theta_true, fitted$vec_theta_hat)

test_hierarchy_layer12(fitted$Beta_hat_plus-fitted$Beta_hat_minus, fitted$Theta_hat, strong = FALSE)


sum(fitted$vec_theta_hat==0)+sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus ==0)
sum(beta_init_lasso_plus-beta_init_lasso_minus ==0) +sum(coefs[range_theta]==0)
coefs[range_theta]==get_vec_theta_hat3(theta_init,l1=l1,l2=l2,l3=l3)

print(fitted$Theta_hat)
fitted$vec_theta_hat
beta_trues[range_theta]

###HIERNET LIBRARY##############################################
#print("-----Hiernet library-----")

fit=hierNet(X_only_main, y_all, lam=15, diagonal = FALSE, stand.main=FALSE,tol=1e-10, center = FALSE)
#fit$th
predicted_lib=predict(fit,X_only_main)
print(r2(y_all, predicted_lib))

sum(get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)==0) + sum( (fit$bp-fit$bn) ==0)

fit$bp - fit$bn

get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)

#print(fit$th)

#fit$bp-fit$bn
#print(sum( get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)==0))
print("hiernet")
all_beta_functions(beta_true, c(fit$bp-fit$bn) )
all_beta_functions(theta_true, get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3))

#get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)[1:40]
#theta_true[1:40]
fit$th





#######


#Residuals their library
r_2way<-y_all-predicted_lib

#Residuals my WHN
predicted<-myWeakHierNet$predict(fitted, X_only_main)
r2(y_all, predicted )
r_2way<-y_all-predicted


##Fit in on residuals

lasso_model_3way <- glmnet(X_3way, r_2way, alpha = 1, lambda=0.05)
#coefficients(lasso_model)

psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model_3way)[-1]),l1=l1,l2=l2,l3=l3) ## psi, get intercept out
print(dim(psi_init))
print(dim(X_3way))
theta_bound<- (fitted$Theta_hat +t(fitted$Theta_hat))*5 ##bound
lambda<-75

source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))
t<-5e-2
myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
                                    l1=l1,l2=l2,l3=l3, scale = FALSE)

# Fit the model
fitted3=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
                            eps=1e-8,l1=l1,l2=l2,l3=l3)

r2(r_2way,myWeakHierNet_seq3$predict(fitted3,X_3way))
#print(fitted$vec_psi_hat)

sum(fitted3$vec_psi_hat==0)
dim(fitted3$vec_psi_hat)

print(sum(beta_true==0))
print(sum(theta_true==0))
print(sum(psi_true==0))

fitted3$psi_hat[8,15,22]

sum(get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)==0) + sum( (fit$bp-fit$bn) ==0) + sum(fitted3$vec_psi_hat==0)

print("r2 for WHN + SEQ2-3")
r2(y_all,predicted_lib + myWeakHierNet_seq3$predict(fitted3,X_3way))

sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus==0)+ sum(fitted$vec_theta_hat==0)+sum(fitted3$vec_psi_hat==0)

###RESULTS FOR PSI#####


all_beta_functions(psi_true, c(fitted3$vec_psi_hat))
test_hierarchy_layer23(fitted$Theta_hat, fitted3$psi_hat)
all_beta_functions(psi_true,coefs[range_psi])

#######

fitted$vec_psi_hat



########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################

#lasso_main <- glmnet(X_only_main, y_all, alpha = 1, lambda=0.005)

#coefs_lasso_main<-coefficients(lasso_main)[-1]
coefs_lasso_main<-coefs[range_main]
all_beta_functions(beta_true, coefs_lasso_main)#lasso main

beta_bound<-coefs_lasso_main*3
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)*1

print(beta_bound)

lmd<-30
t<-5e-6
r_main<-y_all- X_only_main%*%coefs_lasso_main -coefficients(lasso_main)[1]
source(file.path(libs_path,'WeakHierNetSeq12_3way.R'))
seq12<-WeakHierNet_seq_2way3(X=X_2way, theta_init=theta_init, y=r_main, beta_bound=beta_bound, lambda=lmd, t=t, tol=1e-8, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3, scale=FALSE)
fit12<-seq12$fit( X=X_2way, theta_init=theta_init, y=r_main, lambda=lmd,beta_bound=beta_bound, t=t, tol=1e-8, max_iter=5000, eps=1e-8,l1=l1,l2=l2,l3=l3)


all_beta_functions(theta_true, fit12$vec_theta_hat)#lasso all

test_hierarchy_layer12(beta_bound, fit12$theta_hat)

sum(fit12$vec_theta_hat==0)

print(fit12$vec_theta_hat)
print(theta_true)
#for (i in c(5:9) )
#{for (j in c(4:9))
#{print(fit12$theta_hat[i,j] + print(fit12$theta_hat[j,i]))
  #print(fit12$theta_hat[i,4] +print(fit12$theta_hat[4,i]))
  #print(fit12$theta_hat[4,j] +print(fit12$theta_hat[j,4]))}}


### COMPARE WITH LASSO ON ALL
all_beta_functions(beta_true, coefs[range_main]) #lasso all
all_beta_functions(theta_true, coefs[range_theta])#lasso all
all_beta_functions(psi_true, coefs[range_psi])#lasso all

##PREPARE FOR SEQ23###############################################3
r_2way<-r_main-seq12$predict(fit12, X_2way )
psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model)[-1]),l1=l1,l2=l2,l3=l3)*0 ## psi, get intercept out
theta_bound<- (fit12$theta_hat +t(fit12$theta_hat))*5 ##bound
lambda<-60
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))
t<-5e-2
myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
                                       l1=l1,l2=l2,l3=l3, scale = FALSE)

# Fit the model
fitted33=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
                               eps=1e-8,l1=l1,l2=l2,l3=l3)

test_hierarchy_layer23(theta_bound, fitted33$psi_hat)

r2(r_2way,myWeakHierNet_seq3$predict(fitted33,X_3way))
#print(fitted$vec_psi_hat)

sum(fitted33$vec_psi_hat==0)
length(fitted33$vec_psi_hat)
dim(fitted33$vec_psi_hat)

print(sum(beta_true==0))
print(sum(theta_true==0))
print(sum(psi_true==0))

fitted33$psi_hat[8,15,22]
all_beta_functions(psi_true, fitted33$vec_psi_hat)#lasso all

print("r2 for LASSO + SEQ1-2 + SEQ2-3")
r2(y_all, coefficients(lasso_main)[1] + X_only_main%*%coefs_lasso_main +X_2way%*%fit12$vec_theta_hat + X_3way%*%fitted33$vec_psi_hat)
print("Number of 0s in Lasso, seq 1-2 and seq2-3")
print(sum(coefficients(lasso_main)[-1] == 0) + sum(fit12$vec_theta_hat==0 )  +sum(fitted33$vec_psi_hat==0)   )

sum(coefficients(lasso_main)[-1] == 0)
sum(fit12$vec_theta_hat==0 )
sum(fitted33$vec_psi_hat==0)
    



