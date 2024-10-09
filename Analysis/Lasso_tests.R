library(glmnet)
#library(hierNet)

r2 <- function(actual, predicted) {

  mean_actual <- mean(actual)

  total_sum_squares <- sum((actual - mean_actual)^2)

  residual_sum_squares <- sum((actual - predicted)^2)
  
  # Calculate R-squared
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}

# Set seed for reproducibility
set.seed(123)

# Define the number of rows
n <- 100

# Create the matrix with specified means
# Since we want to create columns with specific means, we generate random data and then adjust the means
X <- matrix(rnorm(n * 8), n, 8)

# Adjust the means of each column
#X[, 1] <- X[, 1] - mean(X[, 1]) + 1
#X[, 2] <- X[, 2] - mean(X[, 2]) + 2
#X[, 3] <- X[, 3] - mean(X[, 3]) + 3
#X[, 4] <- X[, 4] - mean(X[, 4]) + 4
#X[, 5] <- X[, 5] - mean(X[, 1]) + 1
#X[, 6] <- X[, 6] - mean(X[, 2]) + 2
#X[, 7] <- X[, 7] - mean(X[, 3]) + 3
#X[, 8] <- X[, 8] - mean(X[, 4]) + 4


pairwise_products <- NULL
for (i in 1:ncol(X)) {
  for (j in i:ncol(X)) {
    new_col <- X[, i] * X[, j]
    pairwise_products <- cbind(pairwise_products, new_col)
  }
}

# Combine the original matrix, pairwise products, and response variable into a data frame
X <- data.frame(X, pairwise_products)



eigs=svd(X)$d
print(max(eigs)/min(eigs))
  eigs
# Verify the means
colMeans(X)

# Create the response variable y
y <- 3+ 1 * X[, 1] + 4 * X[, 2] + 7 * X[, 3] - 3 * X[, 4] +5*X[,5]-7*X[,6] + 3*X[,12] -4*X[,13] +5*X[,14] -7*X[,23] +8*X[,24]  +10*X[,30] -11*X[,31]+8*X[,32]
X[,4] 

apply(X,2, mean)
lasso_model <- glmnet(X,y , alpha = 1, lambda=0.5, standardize =TRUE)
coefs<-coefficients(lasso_model)[-1]
interc<-coefficients(lasso_model)[1]
coefficients(lasso_model)


print(r2(y, as.matrix(X)%*%matrix(coefs, ncol=1)+interc))


###HIERNET ON UNCENTERED

fit=hierNet(as.matrix(X), y, lam=20, diagonal = FALSE, stand.main=FALSE,tol=1e-10)
#fit$th
predicted_lib=predict(fit,X)
r2(y, predicted_lib)

fit$bp-fit$bn



### SAME BUT CENTERED X

n <- 100

# Create the matrix with specified means
# Since we want to create columns with specific means, we generate random data and then adjust the means
X <- matrix(rnorm(n * 8), n, 8)

# Adjust the means of each column
X[, 1] <- X[, 1] - mean(X[, 1]) 
X[, 2] <- X[, 2] - mean(X[, 2]) 
X[, 3] <- X[, 3] - mean(X[, 3]) 
X[, 4] <- X[, 4] - mean(X[, 4]) 
X[, 5] <- X[, 5] - mean(X[, 1]) 
X[, 6] <- X[, 6] - mean(X[, 2]) 
X[, 7] <- X[, 7] - mean(X[, 3]) 
X[, 8] <- X[, 8] - mean(X[, 4]) 


pairwise_products <- NULL
for (i in 1:ncol(X)) {
  for (j in i:ncol(X)) {
    new_col <- X[, i] * X[, j]
    pairwise_products <- cbind(pairwise_products, new_col)
  }
}

# Combine the original matrix, pairwise products, and response variable into a data frame
X <- data.frame(X, pairwise_products)



# Verify the means
colMeans(X)
eigs=svd(X)$d
print(max(eigs)/min(eigs))

# Create the response variable y
y <- 1 * X[, 1] + 4 * X[, 2] + 7 * X[, 3] - 3 * X[, 4] +5*X[,5]-7*X[,6] + 3*X[,12] -4*X[,13] +5*X[,14] -7*X[,23] +8*X[,24]  +10*X[,30] -11*X[,31]+8*X[,32]



lasso_model <- glmnet(X,y , alpha = 1, lambda=1)
X[,1]
coefs<-coefficients(lasso_model)[-1]
coefficients(lasso_model)


