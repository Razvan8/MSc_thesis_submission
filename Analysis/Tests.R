###DATA PREP



library(Metrics)

library(hierNet)
library(caret)
library(dplyr)
library(Metrics)


libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
data<-load_deoxi_flourination()


print(levels(data$a))# 37
print(levels(data$b))# 4
print(levels(data$s))# 5



x.all<-data[,c('a','b','s')]
y.all<-data$p



library(caret)

dummy <- dummyVars(" ~ .", data=x.all)
X.all<- as.matrix(data.frame(predict(dummy, newdata = x.all))) *1

dim(X.all)
X.all

index <- createDataPartition(y = y.all, p = 0.7, list = FALSE)
# Separate X_train, X_test, y_train, y_test
X_train <- X.all[index, ]
y_train <- y.all[index]
X_test <- X.all[-index,]
y_test <- y.all[-index]




print("-----Hiernet library-----")
##On train
fit=hierNet(X.all, y.all, lam=1e-8, diagonal = FALSE, stand.main=FALSE, tol=1e-11)
#fit
y_pred_train=predict(fit,X.all)
y_pred_test=predict(fit, X_test)
sum(fit$bp+fit$bn)
sum(abs(fit$th))
print(fit$th)
#fit $bn

print("RMSE on only 2 mains with library: ")
rmse(y_test, y_pred_test)
print("R2 on only 2 mains with library: ")
r2(y.all,y_pred_train)
r2(y_test, y_pred_test)
y_pred_test


###MY library
print("-----MY HIERNET-----")

t<-5e-4
#t<-0.001
lmd<-0
tol<-1e-9
#colnames(X.all)


source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
myWeakHierNet<-WeakHierNetUnscaled (X=X.all, Beta_plus_init=matrix(0,nrow = dim(X.all)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X.all)[2], ncol = 1), 
                                    Theta_init=matrix(0, ncol = dim(X.all)[2], nrow = dim(X.all)[2]), y=y.all, lambda=lmd, t=t, tol=tol, 
                                    max_iter=10000, eps=1e-8, center=FALSE, standard_form=FALSE)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X.all, Beta_plus_init=matrix(0,nrow = dim(X.all)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X.all)[2], ncol = 1), 
                         Theta_init=matrix(0, ncol = dim(X.all)[2], nrow = dim(X.all)[2]), y=y.all, lambda=lmd, t=t, tol=tol, 
                         max_iter=10000, eps=1e-8)

#fitted$Beta_hat_plus
sum(fitted$Beta_hat_minus+fitted$Beta_hat_plus)
sum(abs(fitted$Theta_hat))

print(fitted$Theta_hat)

print(lmd)

# Make predictions
new_X <- X.all
print("R2 score on all data")
myWeakHierNet$R2_score(self=fitted, new_X= as.matrix(X.all), y_true = y.all, verbose = TRUE)

myWeakHierNet$predict(fitted,X.all)

print(fitted$Theta_hat)



print(x.all)
print(levels(X$a))
print(levels(X$b))
print(levels(X$s))
