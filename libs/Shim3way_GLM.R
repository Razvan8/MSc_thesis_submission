library(glmnet)
libs_path<-file.path("..","libs")
library(lars)
library(polynom)
library(pls)
source(file.path(libs_path,"Create_synthetic_datasets.R"))
source(file.path(libs_path,"helper_functions.R"))

# Define the function
split_data_safe <- function(X, y, specified_columns, additional_percentage = 0.3, seed=1) {
  
  set.seed(seed)
  # Ensure at least one sample with 1 from each specified column is in the training set
  initial_train_indices <- unique(unlist(lapply(specified_columns, function(col) {
    sample(which(X[, col] == 1), 1)
  })))
  print(length(initial_train_indices))
  
  # Randomly add some percentage of additional samples with value 1 to the training set
  additional_indices <- setdiff((1:dim(X)[1]), initial_train_indices)
  num_additional <- round(length(additional_indices) * additional_percentage)
  
  set.seed(seed)
  additional_train_indices <- sample(additional_indices, num_additional)
  
  # Combine the indices to form the final train indices
  train_indices <- unique(c(initial_train_indices, additional_train_indices))
  
  
  # Combine indices to form the final train and test sets
  
  test_indices <- setdiff( (1:nrow(X) ), train_indices )
  
  # Create train and test sets
  X_train <- X[train_indices, ]
  X_test <- X[test_indices, ]
  y_train <- y[train_indices]
  y_test <- y[test_indices]
  
  # Return the results as a list
  list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test)
}


# Define the function
split_data_basic <- function(X, y, p = 0.8) {
  set.seed(42)  # For reproducibility
  split <- createDataPartition(y, p = p, list = FALSE)
  
  # Create train and test sets
  X_train <- X[split, ]
  X_test <- X[-split, ]
  y_train <- y[split]
  y_test <- y[-split]
  
  # Return the results as a list
  list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test)
}



psi_value_from_table_position<-function (table,i,j,k)
{return( (table[i,j,k] + table[i,k,j] + table [j,i,k] +table[j,k,i] + table[k,i,j] + table[k,j,i] )/6)}


table_position_to_vector_index3<- function(position_tuple, l1=36,l2=3,l3=4) ## takes into account / works only for possible combinations!!!!
{
  
  l12<-l1+l2
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  
  assert(x<=l1, "x should be <=l1")
  assert(y<=l1+l2, "y should be <=l1+l2")
  assert(y>l1, "y should be >l1")
  assert(z>l1+l2, 'z should be >l1+l2')
  
  position_psi<-(x-1)*l2*l3 + (y-l1-1)*l3 + (z-l12) 
  
  return(position_psi)
  
}






get_theta_vec_2way3<-function(Theta_hat,l1=36,l2=3,l3=4)
{range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
counter<-0
vec_theta<-array(0,l1*(l2+l3)+l2*l3)

## case 1 :a with b or s}
for (i in range1)
{for (j in c(range2,range3))
{counter<-counter+1
vec_theta[counter]<-(Theta_hat[i,j] +Theta_hat[j,i])/2
}}

## case 2: b with s
for (i in range2)
{for (j in range3)
{counter<-counter+1
vec_theta[counter]<-(Theta_hat[i,j] +Theta_hat[j,i])/2
}}

return(vec_theta)
}





get_psi_vec3<-function(psi,l1=36,l2=3,l3=4)
{
  assert(all(dim(psi)==l1+l2+l3), "Dimensions are not ok")
  
  psi_vec<-array(0, dim=c(l1*l2*l3) )
  
  
  for (i in c(1:l1)) { #alc
    for (j in c((l1+1):(l1+l2) ) ) { 
      for (k in c( (l1+l2+1): (l1+l2+l3) ) ) {  
        psi_vec[table_position_to_vector_index3(c(i,j,k),l1=l1,l2=l2,l3=l3)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  return(psi_vec)
}






get_psi_from_psi_vec3<-function(psi_vec,l1=21,l2=14, l3=2)
{ counter<-1
psi<-array(0,dim=(c(l1+l2+l3,l1+l2+l3,l1+l2+l3)))

for (i in c(1:l1)) { #alcohol
  for (j in c((l1+1):(l1+l2) ) ) { #base
    for (k in c( (l1+l2+1): (l1+l2+l3) ) ) { #s
      #cat("i:", i, ", j:", j, ", k:", k ,'\n')
      psi[i,j,k]<-psi_vec[counter]/6
      psi[i,k,j]<-psi_vec[counter]/6
      psi[k,i,j]<-psi_vec[counter]/6
      psi[k,j,i]<-psi_vec[counter]/6
      psi[j,i,k]<-psi_vec[counter]/6
      psi[j,k,i]<-psi_vec[counter]/6
      counter<-counter+1
    }}} 

#print(counter)
return(psi)
}





get_theta_from_theta_vec_2way3<-function(vec_theta,l1=21,l2=14, l3=2)
{ counter<-1
range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
Theta_hat<-matrix(0,nrow=l1+l2+l3, ncol=l1+l2+l3)

## case 1 :a with b or s}
for (i in range1)
{for (j in c(range2,range3))
{
  Theta_hat[i,j]<-vec_theta[counter]
  Theta_hat[j,i]<-vec_theta[counter] 
  counter<-counter+1
}}

## case 2: b with s
for (i in range2)
{for (j in range3)
{
  Theta_hat[i,j]<-vec_theta[counter] 
  Theta_hat[j,i]<-vec_theta[counter] 
  counter<-counter+1
}}
assert(counter==l1*l2+l2*l3+l3*l1+1, 'smth wrong with counter')
return(Theta_hat)
}




## FUNCTIONS USED FOR THE CLASS

r2 <- function(actual, predicted) {
  # Calculate the mean of the actual values
  mean_actual <- mean(actual)
  
  # Calculate the total sum of squares
  total_sum_squares <- sum((actual - mean_actual)^2)
  
  # Calculate the residual sum of squares
  residual_sum_squares <- sum((actual - predicted)^2)
  
  # Calculate R-squared
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}

assert <- function(condition, message) {
  if (!condition) stop(message)
}

get_ranges<-function(l1,l2,l3)
{range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
return(list(range_main, range_theta, range_psi))}


get_range3<- function(x,l1=36,l2=3,l3=4)
{if (x<=l1)
{return(c(1:l1))}
  if (x<=l1+l2)
  {return(c( (l1+1) : (l1+l2) ))}
  
  return(c( (l1+l2+1) : (l1+l2+l3) ))
  
}


Soft_thresholding <- function(c, lambda) {
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  assert(lambda >= 0, "lambda cannot be negative.")
  
  # Apply soft thresholding component-wise
  result <- sign(c) * pmax(abs(c) - c(lambda), 0)
  
  return(result)
}


lasso_1d_closed_form<-function(X, y, lambda, w=1, scaled=FALSE)
{ xty<-sum(X*y)
  xtx<-sum(X*X)
  if (xtx ==0)
  {return(0)}
  c<-lambda*w
  if (scaled ==TRUE)
  {c<-c*length(X)}
  #cat(" xtx: ", xtx, " xty: ", xty, " c: ", c )
  #cat(" (xty-c) / xtx:", (xty-c)/xtx, xty/xtx -c/xtx)
  
  result<-Soft_thresholding(xty, c/2)/xtx
  return(result)
}




# Kappa functions
kappa0 <- function(x) {
  k0 <- rep(0,length(x))
  tt <- which((x<=700)&(x!=0))
  k0[tt] <- log((exp(x[tt])-1)/x[tt])
  tt <- which(x>700)
  k0[tt] <- x[tt] - log(x[tt])
  return(k0)
}

nu=10
yg=0.9
exp(yg*nu-kappa0(nu))

kappa1 <- function(x) {
  k1 <- rep(1/2,length(x))
  tt <- which(abs(x)<=0.0001)
  k1[tt] <- 1/2 + x[tt]/12 - x[tt]^3/720 + x[tt]^5/30240
  tt <- which(abs(x)>0.0001)
  k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
  return(k1)
}

kappa2 <- function(x) {
  k2 <- rep(1/12,length(x))
  tt <- which(abs(x)<=0.015)
  k2[tt] <- 1/12 - x[tt]^2/240 + x[tt]^4/6048
  tt <- which(abs(x)>0.015)
  k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
  return(k2)
}


g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN=function(v) min(max(v,0.001),0.999))
  g <- 3.5*tan(pi*(2*tt-1)/2) 
  return(g)
}



irlasso.cb <- function(X, Y, lambda, w.lambda=NULL, beta0=NULL,
                       centering=TRUE, scaling=TRUE, intercept=TRUE,
                       maxit=10, tol=0.0545, sd.tol=1e-6,
                       verbose=FALSE,C=0){
  
  if (verbose) print("Performing IRLASSO-NEW")
  
  # CB Likelihood
  kappa0 <- function(x) {
    k0 <- rep(0,length(x))
    tt <- which((x<=700)&(x!=0))
    k0[tt] <- log((exp(x[tt])-1)/x[tt])
    tt <- which(x>700)
    k0[tt] <- x[tt] - log(x[tt])
    return(k0)
  }
  
  kappa1 <- function(x) {
    k1 <- rep(1/2,length(x))
    tt <- which(abs(x)<=0.0001)
    k1[tt] <- 1/2 + x[tt]/12 - x[tt]^3/720 + x[tt]^5/30240
    tt <- which(abs(x)>0.0001)
    k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
    return(k1)
  }
  
  kappa2 <- function(x) {
    k2 <- rep(1/12,length(x))
    tt <- which(abs(x)<=0.015)
    k2[tt] <- 1/12 - x[tt]^2/240 + x[tt]^4/6048
    tt <- which(abs(x)>0.015)
    k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
    return(k2)
  }
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    #mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    #Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), sd.tol)})
    #sd.Y <- apply(Y, 2, function(v){max(sd(v), sd.tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), sd.tol)})
    #Y <- apply(Y, 2, function(v){v/max(sd(v), sd.tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X)
  }
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Initialization
  if (is.null(beta0)){ beta0 <- rep(0, p) }
  
  beta.old <- array(beta0, dim=c(p, 1, length(lambda)))
  Z.old <- array(0, dim=c(n, 1, length(lambda)))
  W.old <- array(0, dim=c(n, n, length(lambda)))
  
  beta <- array(beta0, dim=c(p, 1, length(lambda)))
  Z <- array(0, dim=c(n, 1, length(lambda)))
  W <- array(0, dim=c(n, n, length(lambda)))
  
  R <- list()
  
  nc <- 1:length(lambda) 
  it.stop <- rep(0, length(lambda))
  
  # Run until convergence or stop
  counter <- 0
  
  repeat{
    if (verbose) print(paste("Performing IRLASSO iter", counter, sep=" "))
    
    # Keep old variables
    beta.old <- beta
    Z.old <- Z
    W.old <- W
    
    # Compute W and Z
    eta <- array(0, dim=c(n, 1, length(lambda)))
    k.p <- array(0, dim=c(n, 1, length(lambda)))
    k.pp <- array(0, dim=c(n, 1, length(lambda)))
    
    if (length(nc)==0) {
      if (verbose) print("No lambda left...")
      break
    }
    
    if (verbose) print(paste("Lambda left", paste(nc, collapse=" ")))
    for (m in nc) {
      
      eta[,,m] <- X%*%beta[,,m] + C
      
      k.p[,,m] <- kappa1(eta[,,m])
      
      k.pp[,,m] <- kappa2(eta[,,m])
      k.pp[,,m][which(k.pp[,,m]<0.005)] <- 0.005
      
      W[,,m] <- diag(as.vector(k.pp[,,m]))
      
      Z[,,m] <- eta[,,m] + diag(1/as.vector(k.pp[,,m]))%*%(Y-as.vector(k.p[,,m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[,,m])%*%X)
      Z.W <- as.matrix(sqrt(W[,,m])%*%Z[,,m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept
      if (is.null(w.lambda)) w.lambda <- rep(1,ncol(X.W))
      fit.lasso <- glmnet(x=X.W, y=Z.W, family="gaussian", alpha=1, lambda=lambda[m],
                          standardize=FALSE, intercept=FALSE, penalty.factor=w.lambda)
      beta[,,m] <- as.numeric(fit.lasso$beta)
      
      # Compute model selection matrix
      if (intercept) {
        s.lasso <- which(as.numeric(fit.lasso$beta)[-1]!=0) # not the intercept, at most 0 <= s.lasso <= p-1
        R[[m]] <- matrix(0, nrow=p, ncol=length(s.lasso))
        R[[m]][1,1] <- 1                                    # always take intercept
        if (length(s.lasso)>0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- 1 + s.lasso[s]
            R[[m]][i.lasso,s] <- 1
          }
        }
      } else {
        s.lasso <- which(as.numeric(fit.lasso$beta)!=0)     # no intercept, at most 0 <= s.lasso <= p
        R[[m]] <- matrix(0, nrow=p, ncol=length(s.lasso))
        if (length(s.lasso)>0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- s.lasso[s]
            R[[m]][i.lasso,s] <- 1
          }
        }
      }
      
    }
    
    epsilon <- sqrt(apply((beta-beta.old)^2, 3, sum)/apply((beta.old)^2, 3, sum) )
    print(paste("Min Divergence", min(epsilon[nc]), sep=" "))
    
    log.like <- apply(beta, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    log.like.ratio <- log.like - apply(beta.old, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    print(paste("Min Loglike ratio", min(log.like.ratio[nc]), sep=" "))
    
    if (sum(is.nan(epsilon[nc]))>0) {
      nan.stop <- which(is.nan(epsilon))
      if (verbose) print(paste("Divergence NaN comps", paste(nan.stop, collapse=" ")))
      for (m in nc) {
        beta[,,m] <- beta.old[,,m]
        Z[,,m] <- Z.old[,,m]
        W[,,m] <- W.old[,,m]
      }
      nc <- setdiff(nc, nan.stop)
    }
    
    if ((min(epsilon[nc])<tol)|(min(log.like.ratio[nc])<tol)) { 
      nc.stop <- which((epsilon<tol)|(log.like.ratio<tol))
      it.stop[nc.stop] <- counter
      if (verbose) print(paste("Divergence/Loglike stop comps", paste(nc.stop, collapse=" ")))
      nc <- setdiff(nc, nc.stop)
    }
    
    if (counter==maxit) { 
      if (verbose) print("Maximum iterarion, no convergence...")
      it.stop[which(it.stop==0)] <- counter
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta*0
  
  for (m in 1:length(lambda)) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta[,,m][2:p]/sd.X                                         # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta[,,m][1:p]/sd.X                                       # no intercept, core
      }
    }
  }
  
  return(list(BETA=beta.tilde, 
              beta=beta, 
              Z=Z, W=W, 
              R=R,
              it=it.stop))
}



cross_validation_irlasso.cb <- function(X, y, lambda_values, l1, l2, l3, k = 3, split_percentage = 0.5) {
  best_lambda <- NULL
  best_R2score <- -Inf
  lambda_scores <- numeric(length(lambda_values))
  
  # Perform k random splits
  for (i in 1:k) {
    # Split the data
     # Change seed for each split to ensure different splits
    split_result <- split_data_safe(X=X, y=y, additional_percentage=split_percentage, specified_columns=unlist(get_ranges(l1=l1, l2=l2, l3=l3)[2]), seed=i*10)
    X_train <- split_result$X_train
    y_train <- split_result$y_train
    X_test <- split_result$X_test
    y_test <- split_result$y_test
    j=0
    for (lambda in lambda_values) {
      j<-j+1
      print(lambda)
      
      res_lasso <- irlasso.cb(X=X_train, Y=y_train, lambda=lambda, w.lambda=NULL, beta0=NULL,
                              centering=FALSE, scaling=FALSE, intercept=T,
                              maxit=10, tol=0.0545, sd.tol=1e-6,
                              verbose=F)
      
      coefs_lasso <- array(res_lasso$beta[-1,1,1])
      interc <- res_lasso$beta[1,1,1]
      pred_test <- kappa1(interc + X_test %*% coefs_lasso)
      plot(pred_test, y_test, xlab = "Predicted Yield", ylab = "True Yield", main = "Predicted vs True Yield")
      #plot(pred_test, y_test, xlab = "Predictions", ylab = "True Values", main = "Predictions vs True Values")
      abline(a = 0, b = 1, col = "red")
      

   
      
      
      # Compute the R2 score
      R2 <- r2(y_test, pred_test)
      
      # Store the R2 score for the current lambda
      lambda_scores[j] <- lambda_scores[j] + R2
    }
  }
  
  # Average the scores across the k splits
  lambda_scores <- lambda_scores / k
  
  # Find the best lambda
  best_index <- which.max(lambda_scores)
  best_lambda <- lambda_values[best_index]
  best_R2score <- lambda_scores[best_index]
  print(lambda_values)
  print(lambda_scores)
  
  return(list("best_lambda" = best_lambda, "best_R2score" = best_R2score))
}




poly_min_sign<-function(coefs, positive) #coefs should be from c0 to c4
{assert(length(coefs)==5, "A 4th deg poly should have 5 coefs")
  

  c4<-coefs[5]
  c3<-coefs[4]
  c2<-coefs[3]
  c1<-coefs[2]
  c0<-coefs[1]
  coefs<-array(coefs)
  poly <- function(x) {
    x_vect<-c(1,x,x^2,x^3,x^4)
    return(sum(x_vect*coefs))
  }
  d3<-4*c4
  d2<-3*c3
  d1<-2*c2
  d0<-c1
  coefs_deriv<-c(d0,d1,d2,d3)
  #print(coefs_deriv)
  roots_deriv <- polyroot(coefs_deriv)
  real_roots<-roots_deriv[abs(Im(roots_deriv))<=1e-14]
  real_roots<-Re(real_roots)
  if (positive == TRUE)
  {if(length(real_roots[real_roots>0])==0)
  {real_roots<-list(1e-10)} #add one as backup if it does not find a root
    real_roots<-real_roots[real_roots>0]}
  if(positive == FALSE)
  {if(length(real_roots[real_roots<0])==0)
  {real_roots<-list(-1e-10)} # add one as backup if it does not find a root
    real_roots<-real_roots[real_roots<0]}
  
  
  # Evaluate the polynomial at the real roots
  values <- sapply(real_roots, function(x) poly(Re(x)))
  #print(values)
  # Find the minimum value and the corresponding root
  min_value <- min(values)
  min_root <- real_roots[which.min(values)]
  
  # Print the results
  #print("root for min:")
  #print(min_root)  # The x-value at which the minimum occurs
  #print("Min value of function")
  #print(min_value)  # The minimum value of the polynomial
  
  return(c(min_root, min_value))
  
  
  }

#poly_min_sign(c(4,4,-3,-2,1), positive = FALSE)



poly_lasso_min<-function(coefs, lambda, old_x=0){ #finds the min for a< a=0, a>0 gets the min out of all
  ##long term take care if 2 with same min smart decision
  ##coefs[1]=a0
  ##Minimum using poly and derivative 0
  coefs_without_lambda<-coefs
  
  #case 1 sign positive
  coefs_poz<-coefs_without_lambda
  coefs_poz[2]<-coefs_poz[2]+lambda
  poz<-poly_min_sign(coefs=coefs_poz, positive = TRUE)
  x_min_poz<-poz[1]
  val_min_poz<-poz[2]
  
  #print("starts neg")

  #case 2 sign negative
  coefs_neg<-coefs_without_lambda
  coefs_neg[2]<-coefs_neg[2]-lambda
  neg<-poly_min_sign(coefs=coefs_neg, positive = FALSE)
  x_min_neg<-neg[1]
  val_min_neg<-neg[2]
  
  #print("start 0")
  
  #case 3 sign is 0
  x_zeros<-0
  val_min_zeros<-coefs[1]
  
  all_mins<-c(val_min_neg, val_min_poz, val_min_zeros)
  all_x_min<-c(x_min_neg, x_min_poz, x_zeros)#
  idx_min<-which.min(all_mins)
  x_min<-all_x_min[idx_min]
  x_min<-unlist(x_min)
  poly_min_val<-min(unlist(all_mins))
  
  ### MINIMUM WITH optimize
  fct<-function(x)
  {res<-coefs[1] + x*(coefs[2] +lambda*sign(x)) + x^2 * (coefs[3]) +x^3 *(coefs[4]) +x^4 *(coefs[5])
  }
  result_optimize <- optimize(fct, interval = c(min(-old_x/2 -1e-5, 5*old_x/2 -1e-5), max(-old_x/2 +1e-5, 3*old_x/2 + 1e-5 ) ))
  x_min_optimize<-result_optimize$minimum
  val_min_optimize<-result_optimize$objective
  
  ###x_min_final between optimize and poly
  final_mins<-c(unlist(poly_min_val), unlist(val_min_optimize))
  final_x_min<-c(x_min, x_min_optimize)
  final_idx_min<-which.min(final_mins)
  x_min<-final_x_min[final_idx_min]
  x_min<-unlist(x_min)

  return(x_min)
  
}









get_coef_from_xyz<-function(x,y,z)
{c4<-sum(z*z)
 c3<-sum(2*sum(x*z))
 c2<-sum(x*x)-2*sum(y*z)
 c1<--2*sum(x*y)
 c0<-sum(y*y)
 return(c(c0,c1,c2,c3,c4))}


#get_coef_from_xyz(x=array(0, dim=10), y=array(1, dim=10), z=c(1,2,1,2,1,2,1,2,1,2))

##position in matrix form to position in vector form 2way
matrix_position_to_vector_index_2way<- function(position_tuple, l1,l2,l3) ## takes into account / works only for possible combinations!!!!
{ x<-position_tuple[1]
  y<-position_tuple[2]
  
  range_x<-get_range3(x,l1=l1,l2=l2,l3=l3)
  range_y<-get_range3(y,l1=l1,l2=l2,l3=l3)
  

  assert(x<= l1+l2, "x should be <=l1+l2")
  assert(x<y, "x<y")
  assert(y>l1, 'y should be >l1')

  
  if( all(range_x == c(1:l1)) ==TRUE ) #ab or ac
  { 
    position_vector<- (x-1)*(l2+l3) +(y-l1)  }
  
  
  if( all ( range_x == c( (l1+1): (l1+l2) ) ) == TRUE )  #bc
  {position_vector<-l1*(l2+l3) + (x-l1-1)*l3 + y- (l1+l2)  } 
  return(position_vector)
  
}


get_positions_2way<-function(ls_positions, l1, l2, l3){
  all_positions<-c()
  for (tuple in ls_positions)
  {pos<-matrix_position_to_vector_index_2way(position_tuple = tuple, l1=l1, l2=l2, l3=l3)
  all_positions<-c(all_positions, pos)}
  return(all_positions)}
  
#ls_pos<-list(c(1,3), c(1,4), c(2,4), c(2,3), c(3,6))  
#get_positions_2way(ls_pos, 2,2,2)




#position in 3 dim table to vector index: works only for possible combinations
table_position_to_vector_index3<- function(position_tuple, l1,l2,l3) ## takes into account / works only for possible combinations!!!!
{
  
  l12<-l1+l2
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  
  assert(x<=l1, "x should be <=l1")
  assert(y<=l1+l2, "y should be <=l1+l2")
  assert(y>l1, "y should be >l1")
  assert(z>l1+l2, 'z should be >l1+l2')
  
  position_psi<-(x-1)*l2*l3 + (y-l1-1)*l3 + (z-l12) 
  
  return(position_psi)
  
}

get_positions_3way<-function(ls_positions, l1, l2, l3)
{all_positions<-c()
for (tuple in ls_positions)
{pos<-table_position_to_vector_index3(position_tuple = tuple, l1=l1, l2=l2, l3=l3)
#print(pos)
all_positions<-c(all_positions, pos)}
return(all_positions)}

#ls_pos<-list(c(1,3,6), c(1,3,5), c(1,4,7), c(2,4,7))
#get_positions_3way(ls_pos, l1=2, l2=2, l3=3)

get_xx.all<-function(X,l1,l2,l3)
{xx.all<- matrix(0, nrow=dim(X)[1], ncol=l1+l2+l3+l1*(l2+l3) +l2*l3 )
xx.all[,c(1:(l1+l2+l3))]<-X
counter<-l1+l2+l3+1

for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2+l3) ) ){
    #print(counter)
    xx.all[,counter]<-X[,i]*X[,j]
    counter<-counter+1}}

for (i in c((l1+1): (l1+l2))){ #bc
  for (j in  c( (l1+l2+1): (l1+l2+l3) ) ){
    xx.all[,counter]<-X[,i]*X[,j]
    counter<-counter+1
    }}

assert(counter== l1+l2+l3+ l1*(l2+l3)+l2*l3+1)
return(xx.all)
}

get_xxx.all<-function(X,l1,l2,l3)
{xxx.all<- matrix(0, nrow=dim(X)[1], ncol=l1+l2+l3+l1*(l2+l3) +l2*l3 + l1*l2*l3)
xxx.all[,c(1:(l1+l2+l3+ l1*(l2+l3)+l2*l3))]<-get_xx.all(X=X, l1=l1, l2=l2, l3=l3)
counter<-l1+l2+l3+l1*(l2+l3)+l2*l3+1

for (i in c(1:l1)){ #abc
  for (j in c ( (l1+1): (l1+l2) ) ){
    for (k in c ( (l1+l2+1): (l1+l2+l3) ) ){
      
      xxx.all[,counter]<-X[,i]*X[,j]*X[,k]
      counter<-counter+1}}}
assert(counter==l1+l2+l3+l1*(l2+l3)+l2*l3+l1*l2*l3+1)
return(xxx.all)
}



get_beta_vec_2way<-function(beta,l1,l2,l3, gamma, only_beta = FALSE) #beta is beta_main
{beta_vec2way<- array(0, dim = l1*(l2+l3) +l2*l3 )
counter<-1
  
  for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2+l3) ) ){
    beta_vec2way[counter]<-beta[i]*beta[j]
    counter<-counter+1}}

  for (i in c((l1+1): (l1+l2))){ #bc
    for (j in  c( (l1+l2+1): (l1+l2+l3) ) ){
      beta_vec2way[counter]<-beta[i]*beta[j]
      counter<-counter+1}}

assert(counter==l1*(l2+l3)+l2*l3+1)
if (only_beta==FALSE)
{beta_vec2way<-beta_vec2way*gamma}
return(beta_vec2way)
}


get_beta_vec_3way<-function(beta_2way,l1,l2,l3, delta, only_beta=FALSE) ## # beta_2way should be final beta_2way; #only_beta means product of beta_2ways (gamma included) without delta
{beta_vec3way<- array(0, dim = l1*l2*l3 )
counter<-1

#Iterate over possible positions
for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2) ) ){
    for (k in c ( (l1+l2+1): (l1+l2+l3) ) ){
    
      beta_vec3way[counter]<-beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(i,j),l1=l1,l2=l2,l3=l3)]*
                           beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(i,k),l1=l1,l2=l2,l3=l3)]*
                           beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(j,k),l1=l1,l2=l2,l3=l3)]
                           
    counter<-counter+1}}}

if (only_beta == FALSE)
{beta_vec3way<-beta_vec3way*delta}


assert(counter==l1*l2*l3+1)

return(beta_vec3way)

}



mains_contribution<-function(X, beta_main, l1,l2,l3)
{ range_main<-unlist(get_ranges(l1,l2,l3)[1])
  mains_contrib<-X[,range_main]%*%beta_main
  return(mains_contrib)}

two_ways_contribution<-function(X, gamma_vec, beta_vec_2way,l1,l2,l3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X, beta vec is only beta
{if (already_multiplied==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))}
  range_2ways<-unlist(get_ranges(l1,l2,l3)[2])
#print(dim())
 two_ways_contrib<- X[,range_2ways]%*%(beta_vec_2way*gamma_vec) ##last multiplication should be elementwise
 return(two_ways_contrib)}

three_ways_contribution<-function(X, delta_vec, beta_vec_3way, l1,l2,l3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{if (already_multiplied==TRUE)
{delta_vec<-array(1, dim=length(delta_vec))}
  range_3ways<-unlist(get_ranges(l1,l2,l3)[3])
three_ways_contrib<-X[,range_3ways]%*%(beta_vec_3way*delta_vec) ##last multiplication should be elementwise
return(three_ways_contrib)}



### g function

g_normal<-function(X, beta, gamma_vec, delta_vec,l1,l2,l3, already_multiplied=TRUE) #bet_2way is without gamma, beta_3way without delta only
  #already multiplied=True means beta already has gamma delta inside
{beta_main<-beta[unlist(get_ranges(l1,l2,l3)[1])]
 beta_2way<-beta[unlist(get_ranges(l1,l2,l3)[2])]
 beta_3way<-beta[unlist(get_ranges(l1,l2,l3)[3])]
 if (already_multiplied==TRUE)
 {gamma_vec<-array(1, dim=length(gamma_vec))
  delta_vec<-array(1, dim=length(delta_vec))}
 result<-mains_contribution(X=X, beta_main = beta_main, l1=l1, l2=l2, l3=l3)+ 
         two_ways_contribution(X=X, gamma_vec=gamma_vec, beta_vec_2way=beta_2way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)+
         three_ways_contribution(X=X, delta_vec = delta_vec, beta_vec_3way = beta_3way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)
 #cat("g:", result[1:10])
 return(result)
}

g_bern<-function(X, beta, gamma_vec, delta_vec,l1,l2,l3, already_multiplied=TRUE) #bet_2way is without gamma, beta_3way without delta only
  #already multiplied=True means beta already has gamma delta inside
{beta_main<-beta[unlist(get_ranges(l1,l2,l3)[1])]
beta_2way<-beta[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta[unlist(get_ranges(l1,l2,l3)[3])]
if (already_multiplied==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))
delta_vec<-array(1, dim=length(delta_vec))}
v<-mains_contribution(X=X, beta_main = beta_main, l1=l1, l2=l2, l3=l3)+ 
  two_ways_contribution(X=X, gamma_vec=gamma_vec, beta_vec_2way=beta_2way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)+
  three_ways_contribution(X=X, delta_vec = delta_vec, beta_vec_3way = beta_3way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)
result<-kappa1(v)
#cat("g:", result[1:10])
return(result)
}



##penalty for 1 vector
get_penalty<-function(vector, weights, lambda, already_weighted=TRUE){
  result=lambda*sum(abs(vector)*abs(weights))
  return(result)
}





##loss function normal- Q
Q_normal<-function(X,y, beta, gamma_vec, delta_vec, lambda_beta, lambda_gamma, lambda_delta, w_beta, w_gamma, w_delta,l1,l2,l3,
                   already_multiplied=TRUE, scaled=FALSE)
{ if (length(beta)== l1 +l2+l3)
{#print("Beta was given only main and computed for the rest")
  already_multiplied = TRUE
  beta_2way<-get_beta_vec_2way(beta = beta, l1=l1, l2=l2, l3=l3, gamma= gamma_vec, only_beta = FALSE )
  beta_3way<-get_beta_vec_3way(beta = beta_2way, l1=l1, l2=l2, l3=l3, delta = delta_vec, only_beta = FALSE)
  beta<-c(beta, beta_2way,beta_3way)}
  #should be before making it 1
  penalty_beta<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[1])], weights=w_beta, lambda = lambda_beta  )
  penalty_gamma<-get_penalty(vector=gamma_vec, weights=w_gamma, lambda = lambda_gamma  )
  penalty_delta<-get_penalty(vector=delta_vec, weights=w_delta, lambda = lambda_delta  )
  
  if (already_multiplied ==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))
delta_vec<-array(1, dim=length(delta_vec))}
 error<-sum((y-g_normal(X=X, beta=beta, gamma_vec = gamma_vec, delta_vec = delta_vec, l1=l1, l2=l2, l3=l3, already_multiplied = already_multiplied))**2)
  if(scaled==TRUE)
   {error<-error/(2*dim(X)[1])}
 loss<- error+penalty_beta+penalty_gamma+penalty_delta
 #cat("err,", error, '  ',penalty_beta,' ',penalty_gamma,' ',penalty_delta )
 return(loss)
}


##loss function normal- Q
Q_bern<-function(X,y, beta, gamma_vec, delta_vec, lambda_beta, lambda_gamma, lambda_delta, w_beta, w_gamma, w_delta,l1,l2,l3,
                   already_multiplied=TRUE, scaled=TRUE, intercept=0)
{ if (length(beta)== l1 +l2+l3)
{#print("Beta was given only main and computed for the rest")
  already_multiplied = TRUE
  beta_2way<-get_beta_vec_2way(beta = beta, l1=l1, l2=l2, l3=l3, gamma= gamma_vec, only_beta = FALSE )
  beta_3way<-get_beta_vec_3way(beta = beta_2way, l1=l1, l2=l2, l3=l3, delta = delta_vec, only_beta = FALSE)
  beta<-c(beta, beta_2way,beta_3way)}
  #should be before making it 1
  penalty_beta<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[1])], weights=w_beta, lambda = lambda_beta  )
  penalty_gamma<-get_penalty(vector=gamma_vec, weights=w_gamma, lambda = lambda_gamma  )
  penalty_delta<-get_penalty(vector=delta_vec, weights=w_delta, lambda = lambda_delta  )
  
  if (already_multiplied ==TRUE)
  {gamma_vec<-array(1, dim=length(gamma_vec))
  delta_vec<-array(1, dim=length(delta_vec))}
  
  
  #def log like: sum(y*(Xbeta)-k(Xbeta))
  #v=g_normal(X=X, beta=beta, gamma_vec = gamma_vec, delta_vec = delta_vec, l1=l1, l2=l2, l3=l3, already_multiplied = already_multiplied) #Xbeta
  v=X%*%beta+intercept
  #cat(y*v[1:5], '    aa   ', kappa0(v)[1:5], '   bbbb   ')
  log.like<-sum(y*v-kappa0(v))
  if(scaled==TRUE) ############# CHECK THIS ###########################
  {log.like<-log.like/(2*dim(X)[1])}
  loss<- -log.like+penalty_beta+penalty_gamma+penalty_delta
  #cat("log.like,", log.like, '  ',penalty_beta,' ',penalty_gamma,' ',penalty_delta )
  return(loss)
}

minimizer_Q_bern<-function(X,y, C, lambda, beta_old, weight=1, scaled=TRUE, intercept=0) #function to find the minimum for 1D beta update # interval should depend on old beta/gamma
{
 
 fct<-function(b)
 {#penalty<-get_penalty(vector=c(b), weights = c(weight), lambda=lambda)
 penalty<-abs(b)*lambda*weight
 v= X*b+C+intercept ##minimize for kappa1(Xbeta+C) =y
 #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
 log.like<-sum(y*v-kappa0(v)) 
 if(scaled==TRUE) ############# CHECK THIS ###########################
 {log.like<-log.like/(2*length(X))}
 loss<- -log.like+penalty
 #cat("log, pen: ", log.like, " ", penalty)
 #cat("loss", loss)
 return(loss)
 }
 #fct(1)
 interval<-c(min(-beta_old/2 -1e-5, 5*beta_old/2 -1e-5), max(-beta_old/2 +1e-5, 3*beta_old/2 + 1e-5 ) )
 #cat("interval",interval)
 result_optimize <- optimize(fct, interval = interval )
 minimum<-result_optimize$minimum
 f_0<-fct(0)
 if ( f_0 <= fct(minimum) & f_0 <=fct(beta_old))
 {return(0)}
 
 if (fct(beta_old)<=fct(minimum))
 {return(beta_old)}
 
 return(minimum)
 }




minimizer_Q_bern_beta<-function(X1, X2, y, C, lambda, beta_old, weight=1, scaled=TRUE, intercept=0) #function to find the minimum for 1D beta update # interval should depend on old beta/gamma
{
  
  fctb<-function(b)
  {penalty<-get_penalty(vector=c(b), weights = c(weight), lambda=lambda)
    #penalty<-abs(b)*lambda*weight
    C<-array(C)
    v= as.matrix(X1*b +X2*(b^2)+C+ intercept) ##minimize for kappa1(Xbeta+C) =y
    #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
    log.like<-sum(y*v-kappa0(v)) 
    #cat("dimX: ", length(X1),' pen:', penalty)
    if(scaled==TRUE) ############# CHECK THIS ###########################
    {log.like<-log.like/(2*length(X1))}
    loss<- -log.like+penalty
    #cat("log, pen: ", log.like, " ", penalty)
    #cat("loss", loss)
    return(loss)
  }
  #fct(1)
  interval<-c(min(-beta_old/2 -1e-5, 5*beta_old/2 -1e-5), max(-beta_old/2 +1e-5, 3*beta_old/2 + 1e-5 ) )
  #cat("interval",interval)
  result_optimize <- optimize(fctb, interval = interval )
  minimum<-result_optimize$minimum
  f_0<-fctb(0)
  if ( f_0 <= fctb(minimum) & f_0 <=fctb(beta_old))
  {return(0)}
  
  if (fctb(beta_old)<=fctb(minimum))
  {return(beta_old)}
  
  return(minimum)
}





###RELATIVE DIFFERENCE

compute_rel_dif<-function(Q_old, Q_new)
{rel_dif<- abs(Q_old-Q_new)/abs(Q_old)
return(rel_dif)}



##### UPDATE DELTA FUNCTION #####



update_intercept<-function(X, y, beta_all, intercept_old) #function to find the minimum for 1D beta update # interval should depend on old beta/gamma
{
  
  fctint<-function(intercept)
  {
  v= as.matrix(X%*%beta_all+ intercept) ##minimize for kappa1(Xbeta+C) =y
  #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
  log.like<-sum(y*v-kappa0(v)) 
  loss<- -log.like #scale does not matter
  return(loss)
  }
  #fct(1)
  interval<-c(min(-intercept_old/2 -1e-5, 5*intercept_old/2 -1e-5), max(-intercept_old/2 +1e-5, 3*intercept_old/2 + 1e-5 ) )
  #cat("interval",interval)
  result_optimize <- optimize(fctint, interval = interval )
  minimum<-result_optimize$minimum
  f_0<-fctint(0)
  if ( f_0 <= fctint(minimum) & f_0 <=fctint(intercept_old))
  {return(0)}
  
  if (fctint(intercept_old)<=fctint(minimum))
  {return(intercept_old)}
  
  return(minimum)
}





update_delta<-function(X, y,beta_hat, gamma_hat, delta_hat, lambda_delta, l1, l2, l3, intercept=0, bind_C=TRUE) 
  ##beta_hat is only for mains
{beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
 beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = TRUE) #This is with gamma WITHOUTdelta
 #y_tilde <- y - mains_contribution(X=X, beta_main = beta_hat, l1=l1,l2=l2, l3=l3) 
           #-two_ways_contribution(X=X, gamma_vec = gamma_hat, beta_vec_2way = beta_2way, l1=l1, l2=l2, l3=l3, already_multiplied = TRUE )
 y_tilde<-y ### IT IS Y IN GLM CASE
 X_3way<-X[,unlist(get_ranges(l1,l2,l3)[3])]
 if (var(beta_3way)==0) #lasso does not work if predictor has variance 0
 {return(beta_3way*0)}
 X_tilde<-matrix(rep(beta_3way, each = nrow(X_3way)), nrow = nrow(X_3way))*X_3way
 X_c<-X[,c( unlist(get_ranges(l1,l2,l3)[1]), unlist(get_ranges(l1,l2,l3)[2]) ) ] #Xc
 beta_c<- c(beta_hat, beta_2way) #beta2way is with delta
 C<-X_c%*%beta_c+intercept#add intercept to C
 X_tilde<-cbind(X_tilde,C)#add C to X_tilde
 #lambda_delta<-lambda_delta/(2*nrow(X)) ##scale lambda because in lasso we have 1/(2n)* loss
 #lasso_model <- glmnet(X_tilde, y_tilde, alpha = 1, lambda = lambda_delta, intercept = FALSE, standardize = FALSE)
 
 delta_hat_old<-delta_hat
 #print(lambda_delta)
 beta0<-0*delta_hat#init beta 0 good
 beta0<-c(beta0,1)
 Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                   lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta  , 
                   w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
 if (bind_C==TRUE){
 lasso_rez<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_delta, w.lambda=NULL, beta0=beta0,
            centering=FALSE, scaling=FALSE, intercept=F,
            maxit=10, tol=0.0545, sd.tol=1e-6,
            verbose=F)
 #print(dim(delta_hat))
 lasso_coef <- array(lasso_rez$BETA, dim= length(lasso_rez$BETA) )
 #cat("lasso coef: ", lasso_coef)
 delta_hat<- lasso_coef[-length(lasso_coef)]# last is coef for all other factors including intercept
 }
 else{
   lasso_rez<-irlasso.cb(X=X_tilde[,-dim(X_tilde)[2]], Y=y_tilde, lambda=lambda_delta, w.lambda=NULL, beta0=beta0,
                         centering=FALSE, scaling=FALSE, intercept=F,
                         maxit=10, tol=0.0545, sd.tol=1e-6,
                         verbose=F, C=C)
   #print(dim(delta_hat))
   delta_hat <- array(lasso_rez$BETA, dim= length(lasso_rez$BETA) )
 }
 
 Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                   lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta, 
                   w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
 #cat(" try delta: new- old: ", Q_new -Q_old)
 
 if ( Q_new-Q_old >abs(Q_old)*1e-10){
   delta_hat<-delta_hat_old #keep the odl one
   print("There might be numerical instability in update delta which was taken care of by using old delta. ")
   cat("Actually new-old = 0")
 }

 #print(delta_hat)
 return(delta_hat)
 
}



###test update delta
#data<- create_basic_dataset()
#X<- data$X
#y<- data$y$true

#beta_true<- data$beta[-1,]
#l1=8
#l2=8
#l3=4

#for (lambda_delta in c( c(0), c(100) ) ) {
  #for (noise in c(0,0.01,1)){

#beta_main<-beta_true[1:(l1+l2+l3)]
#beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
#beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
#beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
#gamma_hat<-beta_2way/beta_2way_without_gamma
#gamma_hat[is.nan(gamma_hat)]<-0
#lambda_delta<-100

#beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)

#delta_true<-beta_3way/beta_3way_without_gamma
#delta_true[is.nan(delta_true)]<-0

#delta_hat<-delta_true+rnorm(length(delta_true),0,noise)
#delta_pred <- update_delta(X=X, y=y,beta_hat=beta_main, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_delta=lambda_delta, l1=l1, l2=l2, l3=l3) 
 


#beta_3way_with_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = delta_hat, only_beta = FALSE)




#q_pred <- Q_normal(X=X,y=y, beta=beta_main, gamma_vec=gamma_hat, delta_vec=delta_pred, 
#         lambda_beta=1, lambda_gamma=1, lambda_delta=lambda_delta, 
#          w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)

#q_before<-Q_normal(X=X,y=y, beta=beta_main, gamma_vec=gamma_hat, delta_vec=delta_hat, 
#         lambda_beta=1, lambda_gamma=1, lambda_delta=lambda_delta, 
#         w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3,already_multiplied=TRUE)

#print("info: ")
#cat(" noise:",noise,"lmd:", lambda_delta)
#print(q_pred-q_before)

#if(q_pred-q_before>0)
#{print("numerical instability !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")}

#print(delta_pred)
#print("pred was now true")
#print(delta_true)

#}}

#print(delta_pred)



##### UPDATE GAMMA FUNCTION #####


update_gamma<-function(X, y,beta_hat, gamma_hat, delta_hat, lambda_gamma, l1, l2, l3, w=1, intercept=0) 
{range1<- c(1:l1)
 range2<-c((l1+1):(l1+l2))
 range3<-c( (l1+l2+1) : (l1+l2+l3) )
 X_main<-X[,c(1:(l1+l2+l3)) ]
 X_2way<-X[,c( (l1+l2+l3+1): (l1+l2+l3+l1*l2+l2*l3+l1*l3) )]
 X_3way<-X[,c( (l1+l2+l3+l1*l2+l2*l3+l1*l3+1):(l1+l2+l3+l1*l2+l2*l3+l1*l3+l1*l2*l3) )]
 beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
 beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
 
 if (w==1)
 {w=array(1, dim=length(gamma_hat))}
 
for(i in range1){
  for (j in range2){
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    discard_2way<-matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)
    ls_tuples_discard_3way<-list()
    
    for (elem3 in range3) {
      tuple <- c(i, j, elem3)
      ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
    }
    #print("ls tuples: ")
    #print(ls_tuples_discard_3way)
    
    discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
    
    #print("ls positions 3 way: ")
    #print(discard_3way)
    X_2way_kept<-X_2way[,-discard_2way]
    X_3way_kept<-X_3way[, -discard_3way]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    #print(X_2way_kept)
    #print(beta_2way_kept)
    #print(X_3way_kept)
    #print(beta_3way_kept)
    
    #y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    #print(y_tilde)
    y_tilde<-y
    three_ways=0
    #print("ok")
    discard_from_c_3way<-c()
    for (k in range3) #compute 3 ways contrib
    {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
                             gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *  
                             gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
                             delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
    discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
    
    }
 
    X_tilde<-X_2way[,discard_2way]*beta_hat[i]*beta_hat[j]+  three_ways
    
    
    X_c<-cbind(X_main,X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
    beta_c<-c(beta_hat, beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
    C<-X_c%*%beta_c
    #X_tilde<-cbind(X_tilde,C)
    
    Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                  lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                  w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
    
    
     
    #rez_lasso<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_gamma, w.lambda=NULL, 
                          #beta0=gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)],
                          #centering=FALSE, scaling=FALSE, intercept=F,
                          #maxit=10, tol=0.0545, sd.tol=1e-6,
                          #verbose=F)
    #cat("lasso1d beta: ",array(rez_lasso$BETA, dim=2), ' ')
    #gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] <- array(rez_lasso$BETA, dim=2)[-length(array(rez_lasso$BETA, dim=2))]
    gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] <- minimizer_Q_bern(
      X=X_tilde,y=y_tilde, C=C, lambda=lambda_gamma, beta_old=gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] 
      , weight=1, scaled=TRUE, intercept = intercept)

    
   
    
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
             lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
             w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3,already_multiplied=TRUE, intercept = intercept)
    
    #if (Q_new-Q_old >=0)
    #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
    #print("gamma")
    #print(Q_new-Q_old)
    if ( Q_new-Q_old >abs(Q_old)*1e-2){
      print("There might be numerical instability in update gamma.")
    }
    
  }
}
 

 print("stat ik")
 
 for(i in range1){
   for (k in range3){
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     discard_2way<-matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)
     ls_tuples_discard_3way<-list()
     
     for (elem3 in range2) {
       tuple <- c(i, elem3, k)
       ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
     }
     #print("ls tuples: ")
     #print(ls_tuples_discard_3way)
     
     discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
     
     #print("ls positions 3 way: ")
     #print(discard_3way)
     X_2way_kept<-X_2way[,-discard_2way]
     X_3way_kept<-X_3way[, -discard_3way]
     beta_2way_kept <- beta_2way[-discard_2way]
     beta_3way_kept <- beta_3way[-discard_3way]
     gamma_hat_kept <- gamma_hat[-discard_2way]
     delta_hat_kept <- delta_hat[-discard_3way]
     
     #y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
     y_tilde<-y
     three_ways=0
     discard_from_c_3way<-c()
     for (j in range2) #compute 3 ways contrib
     {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
       gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
       gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
       delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
     discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
     
     }
     
     X_tilde<-X_2way[,discard_2way]*beta_hat[i]*beta_hat[k]+  three_ways
     
     
     X_c<-cbind(X_main,X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
     beta_c<-c(beta_hat, beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
     C<-X_c%*%beta_c
     #X_tilde<-cbind(X_tilde,C)
     
     Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                     lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                     w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
     
     
     
     #lasso_1d_closed_form(X=X_tilde, y= y_tilde, lambda=lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)] )
     #print(lambda_gamma)
     #gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)]<-lasso_1d_closed_form(X=X_tilde, y= y_tilde, lambda=lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3) ] )
     
     #rez_lasso<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_gamma, w.lambda=NULL, beta0=NULL,
                           #centering=FALSE, scaling=FALSE, intercept=T,
                           #maxit=10, tol=0.0545, sd.tol=1e-6,
                           #verbose=F)$BETA

     
     #gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)] <- array(rez_lasso$BETA, dim=1)[-1]
     gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)] <- minimizer_Q_bern(
       X=X_tilde,y=y_tilde, C=C, lambda=lambda_gamma, beta_old=gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)] 
       , weight=1, scaled=TRUE, intercept = intercept)
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     Q_new <-  Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
     
     #if (Q_new-Q_old >=0)
     #print("gamma")
     #print(Q_new-Q_old)
     #cat(" new-old, old delta: ", Q_new-Q_old, " ", Q_old)
     if ( Q_new-Q_old >abs(Q_old)*1e-10){
       print("There might be numerical instability in update gamma.")
     }
     
   }
 }
 

 #cat("gamma : ", gamma_hat)

 
 for(j in range2){
   for (k in range3){
     
     discard_2way<-matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3)
     ls_tuples_discard_3way<-list()
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     
     for (elem3 in range1) {
       tuple <- c(elem3, j, k)
       ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
     }
     #print("ls tuples: ")
     #print(ls_tuples_discard_3way)
     
     discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
     
     #print("ls positions 3 way: ")
     #print(discard_3way)
     X_2way_kept<-X_2way[,-discard_2way]
     X_3way_kept<-X_3way[, -discard_3way]
     beta_2way_kept <- beta_2way[-discard_2way]
     beta_3way_kept <- beta_3way[-discard_3way]
     gamma_hat_kept <- gamma_hat[-discard_2way]
     delta_hat_kept <- delta_hat[-discard_3way]
     
     #y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
     y_tilde<-y
     three_ways=0
     discard_from_c_3way<-c()
     for (i in range1) #compute 3 ways contrib
     {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
       gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
       gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
       delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
     discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
     
     }
     
     X_tilde<-X_2way[,discard_2way]*beta_hat[j]*beta_hat[k]+  three_ways
     
     
     X_c<-cbind(X_main,X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
     beta_c<-c(beta_hat, beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
     C<-X_c%*%beta_c
     #X_tilde<-cbind(X_tilde,C)
     
     Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                     lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                     w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
     
     
     
      
    # rez_lasso<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_gamma, w.lambda=NULL, beta0=NULL,
                           #centering=FALSE, scaling=FALSE, intercept=T,
                           #maxit=10, tol=0.0545, sd.tol=1e-6,
                           #verbose=F)$BETA
     
     
     gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3)] <- minimizer_Q_bern(
       X=X_tilde,y=y_tilde, C=C, lambda=lambda_gamma, beta_old=gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3)] 
       , weight=1, scaled=TRUE, intercept = intercept)
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                     lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                     w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
     #if (Q_new-Q_old >=0)
     #print("gamma")
     #print( Q_new-Q_old)
     if ((Q_new-Q_old)> abs(Q_old)*1e-3 )
     {print("There might be numerical instability in gamma.")
       }
     
     
   }
 }
 
return(gamma_hat)
  
}


  
#data<- create_basic_dataset()
#X<- data$X
#y<- data$y$true

#beta_true<- data$beta[-1,]
#l1=8
#l2=8
#l3=4
#print(beta_true)

#beta_main<-beta_true[1:(l1+l2+l3)]
#beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
#beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
#beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
#beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)

#gamma_hat<-beta_2way/beta_2way_without_gamma
#gamma_hat[is.nan(gamma_hat)]<-0
#lambda_gamma<-1


#gamma_true<-beta_2way/beta_2way_without_gamma
#gamma_true[is.nan(gamma_true)]<-0

#delta_true<-beta_3way/beta_3way_without_gamma
#delta_true[is.nan(delta_true)]<-0

#delta_hat<-delta_true

#lambda_gamma<-2e3

#gamma_hat<-gamma_true+rnorm(length(gamma_hat), 0,0.1)
#gamma_hat[1]<-30
#gamma_hat
#gamma_pred <- update_gamma(X=X, y=y,beta_hat=beta_main, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_gamma=lambda_gamma,  l1=l1, l2=l2, l3=l3, w=1) 
#gamma_pred

#sum(gamma_pred==0)  
#sum(abs(gamma_pred))


#sum(abs(gamma_true))
#sum(gamma_true==0)




##### UPDATE BETA FUNCTION #####


update_beta <- function(X, y, beta_hat, gamma_hat, delta_hat, lambda_beta, l1, l2, l3, w=1, intercept=0)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  X_main<-  X[, c( 1:(l1 + l2 + l3 ) )]
  X_2way <- X[, c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l1 * l2 + l2 * l3 + l1 * l3))]
  X_3way <-X[, c((l1 + l2 + l3 + l1 * l2 + l2 * l3 + l1 * l3 + 1):(l1 + l2 + l3 + l1 *l2 + l2 * l3 + l1 * l3 + l1 * l2 * l3))]
  beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
  beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
  
  if (w == 1)
  {
    w = array(1, dim = length(gamma_hat))
  }
  
  ##CASE 1 Iterate for i
  for (i in range1) {
    discard_main <- c(i)
    #print(i)
    
    #discard 2way
    ls_pos_2way <- list()
    for (jk in c(range2, range3))
      #positions ij ik in this order
    {
      ls_pos_2way <- append(ls_pos_2way, list(c(i, jk)) )}
    #print("before get_poss")
    #print(ls_pos_2way)
    discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
    #print("after get_poss")
    
    
    #discard 3way
    ls_pos_3way<-list()
    for (j in range2) {
      for (k in range3) {
        ls_pos_3way <- append(ls_pos_3way, list(c(i, j, k)) )
        
      }
    }
    discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
    
    #print("X3way")
    #print(X_3way)
    #print("dicard_3way")
    #print(discard_3way)
    
    X_main_kept<-X_main[,-discard_main]
    X_2way_kept <- X_2way[, -discard_2way]
    X_3way_kept <- X_3way[,-discard_3way]
    
    beta_main_kept<-beta_hat[-discard_main]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    
    
    #y_tilde<-y - X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    #print("ytilde")
    #print(y_tilde)
    
    #Create X_tilde from X1_tilde X2_tilde
    two_ways<-0
    for (jk in c(range2, range3)) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(i,jk),l1=l1, l2=l2, l3=l3)]*(beta_hat[jk])*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,jk), l1=l1, l2=l2 ,l3=l3)] }
    
    three_ways<-0
    discard_from_c_3way<-c()
    for (j in range2) #compute 3 ways contrib
    {for(k in range3)
    { 
      three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[j]*beta_hat[k])^2)*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
      gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
      delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
      discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
    
    }}
    #cat("X[, discard main]: ", X[,discard_main], "dim2ways", two_ways)
    #print(class(X_main[,discard_main]))
    #print(class(two_ways))
    
    
    
    
  
    
    X1_tilde<-  X_main[,discard_main] + array(two_ways)
    X2_tilde<-array(three_ways)

    y_tilde<-y
    X_c<-cbind(X_main[,-discard_main],X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
    beta_c<-c(beta_hat[-discard_main], beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
    C<-X_c%*%beta_c
    
    #cat("X1dim: ", dim(X1_tilde), " X2dim: ", dim(X2_tilde))
    Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)

    
    #coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde) ) #coefs as c0 c1 c2...c4
    beta_hat_old<-beta_hat[i]
    #beta_hat[i]<-poly_lasso_min(coefs = coefs, lambda = lambda_beta, old_x = beta_hat_old ) #beta updated
    beta_hat[i]<-minimizer_Q_bern_beta(X1=X1_tilde, X2=X2_tilde, y=y_tilde, C=C, lambda=lambda_beta, beta_old=beta_hat_old,
                                                 weight=1, scaled=TRUE, intercept = intercept)

    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    
    Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
    if ( Q_new-Q_old >abs(Q_old)*1e-10){
      print("There might be numerical instability in update beta.")
      cat(" new-old beta: ",Q_new- Q_old,' ')
    }

    print(Q_new-Q_old) 
    
    
  }
  

  
  #print("start j")
  
  ##CASE 2 Iterate for j
  for (j in range2) {
    discard_main <- c(j)
    
    #discard 2way
    ls_pos_2way <- list()
    for (i in range1) #ij
    {ls_pos_2way <- append(ls_pos_2way, list(c(i,j)) )}
    for (k in range3) #ij
    {ls_pos_2way <- append(ls_pos_2way, list( c(j,k)) )}
    
    discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
    
    
    #discard 3way
    ls_pos_3way<-list()
    for (i in range1) {
      for (k in range3) {
        ls_pos_3way <- append(ls_pos_3way, list(c(i, j, k) ))
        
      }
    }
    #print(ls_pos_3way)
    discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
    
    X_main_kept<-X_main[,-discard_main]
    X_2way_kept <- X_2way[, -discard_2way]
    X_3way_kept <- X_3way[,-discard_3way]
    
    beta_main_kept<-beta_hat[-discard_main]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    y_tilde<-y - X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    

    
    #Create X_tilde from X1_tilde X2_tilde
    two_ways<-0
    for (i in range1) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(i,j),l1=l1, l2=l2, l3=l3)]*(beta_hat[i])*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] }
    for (k in range3) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(j,k),l1=l1, l2=l2, l3=l3)]*(beta_hat[k])*
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] }
   
    
    three_ways<-0
    discard_from_c_3way<-c()
    for (i in range1) #compute 3 ways contrib
    {for(k in range3)
    {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[k])^2)*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
      gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
      delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
    discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
    
    }}
    
    X1_tilde<- X[,discard_main] + two_ways
    X2_tilde<-array(three_ways)
    
    y_tilde<-y
    X_c<-cbind(X_main[,-discard_main],X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
    beta_c<-c(beta_hat[-discard_main], beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
    C<-X_c%*%beta_c
    
    
    
    Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)

    
    #coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde)) #coefs as c0 c1 c2...c4
    beta_hat_old<-beta_hat[j]
    beta_hat[j]<-minimizer_Q_bern_beta(X1=X1_tilde, X2=X2_tilde, y_tilde, C=C, lambda=lambda_beta, beta_old=beta_hat_old,
                                                 weight=1, scaled=TRUE, intercept = intercept)

    
    #update beta_23way
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
    if ( Q_new-Q_old >abs(Q_old)*1e-10){
      print("There might be numerical instability in update beta.")
      cat(" new-old beta: ",Q_new-Q_old,' ')
    }
    
    print(Q_new-Q_old) 
    
  }
    
    
    
    ##CASE 3 Iterate for k
    for (k in range3) {
      discard_main <- c(k)
      
      #discard 2way
      ls_pos_2way <- list()
      for (ij in c(range1, range2))
        #positions ij ik in this order
      {
        ls_pos_2way <- append(ls_pos_2way, list( c(ij, k)) )}
      discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
      
      
      #discard 3way
      ls_pos_3way<-list()
      for (i in range1) {
        for (j in range2) {
          ls_pos_3way <- append(ls_pos_3way, list( c(i, j, k)) )
          
        }
      }
      discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
      
      X_main_kept<-X_main[,-discard_main]
      X_2way_kept <- X_2way[, -discard_2way]
      X_3way_kept <- X_3way[,-discard_3way]
      
      beta_main_kept<-beta_hat[-discard_main]
      beta_2way_kept <- beta_2way[-discard_2way]
      beta_3way_kept <- beta_3way[-discard_3way]
      gamma_hat_kept <- gamma_hat[-discard_2way]
      delta_hat_kept <- delta_hat[-discard_3way]
      y_tilde<-y -X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
      
      #Create X_tilde from X1_tilde X2_tilde
      two_ways<-0
      for (ij in c(range1, range2)) #compute 3 ways contrib
      {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(ij,k),l1=l1, l2=l2, l3=l3)]*(beta_hat[ij])*
        gamma_hat[matrix_position_to_vector_index_2way(c(ij,k), l1=l1, l2=l2 ,l3=l3)] }
      
      three_ways<-0
      discard_from_c_3way<-c() ##might be redundant
      for (i in range1) #compute 3 ways contrib
      {for(j in range2)
      {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j])^2)*
        gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
        gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
        gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
        delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
      discard_from_c_3way<-c(discard_from_c_3way,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3))
      
      }}
      
      X1_tilde<- X[,discard_main] + two_ways
      X2_tilde<-array(three_ways)
      
      y_tilde<-y
      X_c<-cbind(X_main[,-discard_main],X_2way[,-discard_2way], X_3way[,-discard_from_c_3way])
      beta_c<-c(beta_hat[-discard_main], beta_2way[-discard_2way], beta_3way[-discard_from_c_3way])
      C<-X_c%*%beta_c
      
     
      
      
      Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma = 0, lambda_delta = 0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
      
      #coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde) ) #coefs as c0 c1 c2...c4
      beta_hat_old<-beta_hat[k]
      beta_hat[k]<-minimizer_Q_bern_beta(X1=X1_tilde, X2=X2_tilde, y_tilde, C=C, lambda=lambda_beta, beta_old=beta_hat_old,
                                                                    weight=1, scaled=TRUE, intercept = intercept)
      #update beta_23way
      beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
      beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
      
      
      
      Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma = 0, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
      if ( Q_new-Q_old >abs(Q_old)*1e-10){
        print("There might be numerical instability in update beta.")
        cat(" new-old beta: ",Q_new-Q_old,' ')
      }
      
      
    }

  #print(beta_hat)
  return(beta_hat)
  
  }
  
  
  



#data<- create_basic_dataset()
#X<- data$X
#y<- data$y$obs

#beta_true<- data$beta[-1,]
#l1=8
#l2=8
#l3=4
#print(beta_true)
#beta_main<-beta_true[1:(l1+l2+l3)]
#beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
#beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
#beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
#beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


#gamma_true<-beta_2way/beta_2way_without_gamma
#gamma_true[is.nan(gamma_true)]<-0
#delta_true<-beta_3way/beta_3way_without_gamma
#delta_true[is.nan(delta_true)]<-0


#delta_hat<-delta_true
#gamma_hat<-gamma_true

#beta_hat<-beta_main+rnorm(length(beta_main), 0,0.5)
#beta_hat[1]<-30

#beta_pred <- update_beta(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_beta=100,  l1=l1, l2=l2, l3=l3, w=1) 
#beta_hat[17:20]
#beta_pred[17:20]
#beta_main[17:20]

#sum(beta_pred==0)  
#sum(abs(beta_pred))


#sum(abs(gamma_true))
#sum(gamma_true==0)











## SHIM CLASS for 3 way
SHIM_3way<-function(X,y, beta_init, gamma_init, delta_init,l1=36,l2=3,l3=4, scale=FALSE)
  
{
  self = list()
  
  self$beta_hat <- beta_init #matrix form
  self$gamma_hat<- gamma_init
  self$delta_hat<-delta_init
  
  self$means_X<-colMeans(X)
  self$stds_X<-apply(X,2,sd)
  self$mean_y<-mean(y)
  self$scale=scale
  
  self$l1=l1
  self$l2=l2
  self$l3=l3
  
  
  
  fit<-function(X, y, lambda_beta, lambda_gamma, lambda_delta, w_beta=NULL, w_gamma=NULL, w_delta=NULL,  tol=1e-5,
                max_iter=50, compute_Q=Q_normal, intercept=0, use_intercept=FALSE, bind_C=TRUE)
  {## STEP 0 (STANDARDIZE)
    if (self$scale == TRUE)
    {      print('was scaled')
      X <- scale(X)}#standardize X
    #y <- scale(y, center = TRUE, scale = FALSE) #let y in 0,1
    
  
    ## STEP 1 (INIT BETA AND GAMMA AND DELTA)
    beta_hat<-self$beta_hat
    gamma_hat<-self$gamma_hat
    delta_hat<-self$delta_hat
    Q_old<-1e100
    
    beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    beta_all<-array(c(self$beta_hat, beta_2way, beta_3way))
    
    
    for (i in c(1:max_iter))  ###print smth IF LOSS DOES NOT DECREASE AFTER ONE ITER
    {    
      if (use_intercept==TRUE)
        print(intercept)
      { Q_old<- compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                         lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta,
                         w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta,l1=self$l1, l2=self$l2, l3=self$l3, intercept = intercept)
        intercept<-update_intercept(X=X, y=y, beta_all=beta_all, intercept_old=intercept)#### later do update(intercept)
      Q_new<- compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta,
                w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta,l1=self$l1, l2=self$l2, l3=self$l3, intercept = intercept)
      if(Q_new -Q_old > abs(Q_old)*1e-5) cat("Might be numerical instability in intercept ", Q_new-Q_old )
      }
      
      ## STEP 2 (UPDATE DELTA)
      delta_hat<- update_delta(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_delta=lambda_delta,
                               l1=self$l1, l2=self$l2, l3=self$l3, intercept=intercept, bind_C = bind_C)
      
      ## STEP 3 (UPDATE GAMMA)
      
      gamma_hat<- update_gamma(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_gamma=lambda_gamma,
                               l1=self$l1, l2=self$l2, l3=self$l3, w=1, intercept = intercept)

      ## STEP 4 (UPDATE BETA)
      beta_hat<- update_beta(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_beta=lambda_beta,  
                             intercept = intercept, l1=self$l1, l2=self$l2, l3=self$l3, w=1)

      ## STEP 5 (COMPUTE REL_DIF)
      Q_new<-compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta,
                       w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta,l1=self$l1, l2=self$l2, l3=self$l3, intercept = intercept)
      print("-----------------------Q_new---------------------------------------")
      #print(Q_new)
      if(Q_new ==Q_old) #rel dif 0 instead of nan
      {#cat("beta_hat: ", beta_hat)
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
        beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
        beta_all<-array(c(self$beta_hat, beta_2way, beta_3way))
        self$intercept<-intercept
        
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "beta_all" = beta_all,
                     "intercept" = self$intercept)) }
      rel_dif<-compute_rel_dif(Q_old=Q_old, Q_new=Q_new)

      
      if ((Q_new-Q_old)>abs(Q_old)* 1e-2 )
      {print("there is numerical instability overall. ")}
      if(i%%1==0)
      {cat("  Q: ",Q_new  )}
      
      
      
      ## STEP 6 (CHECK IF SMALL ENOUGH TO STOP)
      if (abs(rel_dif)<=tol){
        #cat("beta_hat: ", beta_hat)
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
        beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
        beta_all<-array(c(self$beta_hat, beta_2way, beta_3way))
        self$intercept<-intercept
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "beta_all" = beta_all,
                     "intercept" =self$intercept)) }
      Q_old<-Q_new #UPDATE Q_old
      }
      cat("It has not converged. The relative difference between last 2 residuals is:", compute_rel_dif(Q_old=Q_old,Q_new=Q_new) )
      self$beta_hat<-beta_hat
      self$gamma_hat<-gamma_hat
      self$delta_hat<-delta_hat
      beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
      beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
      beta_all<-array(c(self$beta_hat, beta_2way, beta_3way))
      self$intercept<-intercept
      return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "beta_all" = beta_all,
                   "intercept" = self$intercept)) 
      }

  
  
  
  predict<-function(self, X_new ,scale=FALSE)
  {if (scale ==TRUE)
  {X<-scale(X_new)
    print("Take care! X was scaled.")}
  beta_all<-self$beta_all
  v<-  X_new%*%beta_all+self$intercept  ### think should be always 0
  y_pred<-kappa1(v)
  return(y_pred)
  }
  
  R2_score<-function(self, X_new, y_true, scale=FALSE, verbose=TRUE)
  {y_pred<-predict(self, X_new, scale=scale)
  if (verbose == TRUE)
  {cat ("r2 score is ", r2(y_true, y_pred))
    #cat(length(y_true), ' ', length(y_pred))
    plot(y_pred, y_true, xlab = "Predictions", ylab = "True Values", main = "Predicted vs True Values")
    #plot(y_pred, y_true, xlab = "Predicted Yield", ylab = "True Yield", main = "Predicted vs True Yield  ")
    abline(a = 0, b = 1, col = "red")}
  
  
  return(r2(y_true, y_pred))}
  
  
  cross_validation <- function( X, y, lambda_values_main, lambda_values_2way, lambda_delta, intercept, split_percentage = 0.5, verbose=TRUE, k=3) {
    # Split the data
    #self is big model SHIM_GLM
    ranges<-get_ranges(l1=self$l1, l2=self$l2, l3=self$l3)
    range_main<- unlist(ranges[1])
    range_teta<- unlist(ranges[2])
    range_psi<- unlist(ranges[3])
    R2_scores <- matrix(0, nrow = length(lambda_values_main), ncol = length(lambda_values_2way))
    for (i in 1:k) {
    split_result <- split_data_safe(X=X, y=y, additional_percentage=split_percentage, specified_columns =unlist( get_ranges(l1=self$l1, l2=self$l2, l3=self$l3)[2]), seed=i )
    X_train <- split_result$X_train
    y_train <- split_result$y_train
    X_test <- split_result$X_test
    y_test <- split_result$y_test
    
    best_lambda <- NULL
    best_R2score <- -Inf
    
    for (j in 1:length(lambda_values_main)) {
      lambda1 <- lambda_values_main[j]
      for (l in 1:length(lambda_values_2way)) {
        lambda2 <- lambda_values_2way[l]

      # Create and fit the model with the current lambda
      fitted <- fit(X=X_train, y=y_train, lambda_beta=lambda1, lambda_gamma=lambda2, lambda_delta=lambda_delta, w_beta=1, w_gamma=1, w_delta=1,  tol=1e-2,
                          max_iter=20, compute_Q=Q_bern, intercept=intercept, use_intercept=TRUE)
      if (verbose==TRUE){
        cat("Percentages of zeros in fitted: ", sum(fitted$beta_all[range_main]==0)/length(range_main), ', ', 
            sum(fitted$beta_all[range_theta]==0)/ length(range_teta),', ' ,sum(fitted$beta_all[range_psi]==0)/length(range_psi), '   ')
      }

      
      # Compute the R2 score
      R2_scores[j,l] <- R2_scores[j,l]+  R2_score(self=fitted, X_new=X_test, y_true=y_test) }} }
      
    # Average the R2 scores over all k splits
    R2_scores <- R2_scores / k
    
    # Find the best lambda combination with the highest average R2 score
    best_index <- which(R2_scores == max(R2_scores), arr.ind = TRUE)
    best_lambda1 <- lambda_values_main[best_index[1]]
    best_lambda2 <- lambda_values_2way[best_index[2]]
    best_R2score <- R2_scores[best_index[1], best_index[2]]
    
    print(R2_scores)
    
    return(list("best_lambda1" = best_lambda1, "best_lambda2" = best_lambda2, "best_R2score" = best_R2score))
  }
  
  return(list( fit = fit, predict = predict, R2_score = R2_score, self = self, cross_validation = cross_validation))
    
}
  





cross_val_pipeline<-function( X, y, lambda, lambda_values_main, 
                             lambda_values_2way, lambda_delta, intercept, l1, l2, l3, split_percentage = 0.5, verbose=TRUE, k=3)
{# Split the data
  #self is big model SHIM_GLM
  ranges<-get_ranges(l1=l1, l2=l2, l3=l3)
  range_main<- unlist(ranges[1])
  range_teta<- unlist(ranges[2])
  range_psi<- unlist(ranges[3])
  R2_scores <- matrix(0, nrow = length(lambda_values_main), ncol = length(lambda_values_2way))
  for (i in 1:k) {
    split_result <- split_data_safe(X=X, y=y, additional_percentage=split_percentage, specified_columns =unlist( get_ranges(l1=l1, l2=l2, l3=l3)[2]), seed=10*i )
    X_train <- split_result$X_train
    y_train <- split_result$y_train
    X_test <- split_result$X_test
    y_test <- split_result$y_test
    
    best_lambda <- NULL
    best_R2score <- -Inf
    
    ####INIT FROM LASSO ####
    res_lasso<-irlasso.cb(X=X_train, Y=y_train, lambda=lambda, w.lambda=NULL, beta0=NULL,
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
    my_shim<-SHIM_3way(X=X_train, y=y_train, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
    
    ###INIT FINISHED---START SHIM
    
    for (j in 1:length(lambda_values_main)) {
      lambda1 <- lambda_values_main[j]
      for (l in 1:length(lambda_values_2way)) {
        lambda2 <- lambda_values_2way[l]
        
        # Create and fit the model with the current lambda
        fitted <- my_shim$fit(X=X_train, y=y_train, lambda_beta=lambda1, lambda_gamma=lambda2, lambda_delta=lambda_delta, w_beta=1, w_gamma=1, w_delta=1,  tol=1e-2,
                      max_iter=20, compute_Q=Q_bern, intercept=interc_init, use_intercept=TRUE)
        if (verbose==TRUE){
          cat("Percentages of zeros in fitted: ", sum(fitted$beta_all[range_main]==0)/length(range_main), ', ', 
              sum(fitted$beta_all[range_theta]==0)/ length(range_teta),', ' ,sum(fitted$beta_all[range_psi]==0)/length(range_psi), '   ')
        }
        
        
        # Compute the R2 score
        R2_scores[j,l] <- R2_scores[j,l]+  my_shim$R2_score(self=fitted, X_new=X_test, y_true=y_test)
        print(R2_scores)}} }
  
  # Average the R2 scores over all k splits
  R2_scores <- R2_scores / k
  
  # Find the best lambda combination with the highest average R2 score
  best_index <- which(R2_scores == max(R2_scores), arr.ind = TRUE)
  best_lambda1 <- lambda_values_main[best_index[1]]
  best_lambda2 <- lambda_values_2way[best_index[2]]
  best_R2score <- R2_scores[best_index[1], best_index[2]]
  
  print(R2_scores)
  
  return(list("best_lambda1" = best_lambda1, "best_lambda2" = best_lambda2, "best_R2score" = best_R2score))
}

























##### results pipeline Lasso ######

#use it on y_centered
results_pipeline_lasso_GLM<-function(X, y, lambda, l1, l2, l3,
                                 beta_main, beta_2way, beta_3way, beta_main_recovered, 
                                 beta_2way_recovered, beta_3way_recovered, threshold = 0, strong = TRUE, use_intercept=FALSE){
  
  
  
  range_main<-c(1: (l1+l2+l3) )
  range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
  range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
  
  

  lasso_rez<-irlasso.cb(X=X, Y=y, lambda=lambda, w.lambda=NULL, beta0=NULL,
                        centering=FALSE, scaling=FALSE, intercept=use_intercept,
                        maxit=10, tol=0.0545, sd.tol=1e-6,
                        verbose=F)
  coefs_all<-lasso_rez$beta
  if (use_intercept == F) {coefs_all<-c(0, coefs_all)}
  
  coefs_lasso<-coefs_all[-1]
  intercept_lasso<-coefs_all[1]
  beta_main_lasso<-coefs_lasso[range_main]
  beta_2way_lasso<-coefs_lasso[range_theta]
  beta_3way_lasso<-coefs_lasso[range_psi]
  beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
  beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)
  
  predict_lasso<- kappa1(intercept_lasso + X%*% coefs_lasso)
  print("r2 lasso")
  print(r2(y, predict_lasso))
  plot(predict_lasso, y)
  print(noquote(""))
  
  
  ######################### RESULTS LASSO #########################################################
  print("------------------ results lasso without recovered -------------------")
  all_beta_functions(beta_main, beta_main_lasso)
  print(noquote(""))
  all_beta_functions(beta_2way, beta_2way_lasso)
  print(noquote(""))
  all_beta_functions(beta_3way, beta_3way_lasso)
  print(noquote(""))
  
  
  ##hierarchy tests
  
  beta_2way_lasso_matrix<-get_theta_from_theta_vec_2way3(beta_2way_lasso,l1=l1,l2=l2, l3=l3)
  beta_3way_lasso_table<-get_psi_from_psi_vec3(beta_3way_lasso,l1=l1,l2=l2, l3=l3)
  
  test_hierarchy_layer12(beta_main_lasso,beta_2way_lasso_matrix, strong = strong)
  print(noquote(""))
  test_hierarchy_layer23(beta_2way_lasso_matrix, beta_3way_lasso_table, strong = strong)
  print(noquote(""))
  
  
  ##### RESULTS LASSO ON RECOVERED PARAMS ##########
  
  beta_main_lasso_recovered<- get_all_beta(beta_main_lasso, l1=l1, l2=l2, l3=l3, threshold = threshold)
  beta_2way_lasso_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_lasso_matrix, l1=l1, l2=l2, l3=l3, threshold = threshold), l1=l1+1, l2=l2+1, l3=l3+1)
  beta_3way_lasso_recovered<- get_psi_vec3( get_all_psi(beta_3way_lasso_table, l1=l1, l2=l2, l3=l3, threshold = threshold) , l1=l1+1, l2=l2+1, l3=l3+1)
  
  
  #beta_main_lasso
  #beta_main_lasso_recovered
  
  print("--------------------- results lasso recovered ------------------------")
  all_beta_functions(beta_main_recovered, beta_main_lasso_recovered)
  print(noquote(""))
  all_beta_functions(beta_2way_recovered, beta_2way_lasso_recovered)
  print(noquote(""))
  all_beta_functions(beta_3way_recovered, beta_3way_lasso_recovered)
  print(noquote(""))
  
  
  test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong = strong)
  print(noquote(""))
  test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                          get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong = strong)
  
  return (list("intercept" = intercept_lasso, "main"= beta_main_lasso , '2way' = beta_2way_lasso, '3way' = beta_3way_lasso))
  
  
}

#lambda<-0.08
#results_pipeline_lasso(X=X, y=y, lambda=lambda, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
#                       beta_3way_recovered, threshold = 1, strong = TRUE)







####pipeline SHIM
results_pipeline_shim_GLM<-function(X, y, beta_main_init, beta_2way_init, beta_3way_init, lambda_beta, lambda_gamma, lambda_delta, l1, l2, l3,
                                beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, beta_3way_recovered, tol=1e-3,  
                                w_beta = 1, w_gamma = 1, w_delta = 1,  scale=FALSE, strong= TRUE, threshold=0, intercept=0,
                                use_intercept=FALSE, bind_C=TRUE)

{
  ##PREPARE FOR SHIM ###################################################################################################\
  
  range_main<-c(1: (l1+l2+l3) )
  range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
  range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
  
  
  
  beta_main_recovered<- get_all_beta(beta_main, l1=l1, l2=l2, l3=l3, threshold = threshold)
  beta_2way_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_matrix, l1=l1, l2=l2, l3=l3, threshold = threshold), l1=l1+1, l2=l2+1, l3=l3+1)
  beta_3way_recovered<- get_psi_vec3( get_all_psi(beta_3way_table, l1=l1, l2=l2, l3=l3, threshold = threshold) , l1=l1+1, l2=l2+1, l3=l3+1)
  
  
  beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
  beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)
  
  
  gamma_true<-beta_2way/beta_2way_without_gamma
  gamma_true[is.nan(gamma_true)]<-0
  delta_true<-beta_3way/beta_3way_without_gamma
  delta_true[is.nan(delta_true)]<-0
  
  
  
  
  
  beta_hat<-beta_main_init
  beta_2way_without_gamma_init<-get_beta_vec_2way(beta_hat,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
  gamma_hat<- beta_2way_init/beta_2way_without_gamma_init
  gamma_hat[is.nan(gamma_hat)]<-0
  gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
  
  #beta_3way_without_delta_init<-get_beta_vec_3way(beta_2way = beta_2way_init, l1=l1, l2=l2, l3=l3, only_beta = TRUE)
  beta_3way_without_delta_init<- get_beta_vec_3way(beta_2way_without_gamma_init*gamma_hat, l1=l1, l2=l2, l3=l3, only_beta = TRUE) #maybe better for shim init
  
  delta_hat<- beta_3way_init/beta_3way_without_delta_init
  delta_hat[!is.finite(delta_hat)]<-0
  delta_hat[is.nan(delta_hat)]<-0
  
  

  
  my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = scale)
  fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                      lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = w_beta, w_gamma = w_gamma, w_delta = w_delta, tol=tol,
                      intercept=intercept, use_intercept=use_intercept, bind_C=bind_C, compute_Q = Q_bern)
  

  my_shim$R2_score(self=fitted, X_new=X, y_true=y )
  
  
  beta_all_shim<-fitted$beta_all
  beta_main_shim<-beta_all_shim[range_main]
  beta_2way_shim<-beta_all_shim[range_theta]
  beta_3way_shim<-beta_all_shim[range_psi]
  intercept_shim<-fitted$intercept
  
  
  print("--------------beta functions and hierarchy before recovery-----------")
  all_beta_functions(beta_main, beta_main_shim)
  print(noquote(""))
  all_beta_functions(beta_2way, beta_2way_shim)
  print(noquote(""))
  all_beta_functions(beta_3way, beta_3way_shim)
  print(noquote(""))
  
  ##hierarchy tests
  beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
  beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1,l2=l2, l3=l3)
  
  test_hierarchy_layer12(beta_main_shim,beta_2way_shim_matrix, strong = strong)
  print(noquote(""))
  test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = strong)
  print(noquote(""))
  
  
  
  ### SHIM FOR RECOVERED DATA #############
  beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
  beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1, l2=l2, l3=l3)
  
  beta_main_shim_recovered<- get_all_beta(beta_main_shim, l1=l1, l2=l2, l3=l3, threshold = threshold)
  beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = threshold), l1=l1+1, l2=l2+1, l3=l3+1)
  beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = threshold) , l1=l1+1, l2=l2+1, l3=l3+1)
  
  
  print("--------------beta functions and hierarchy after recovery-----------")
  all_beta_functions(beta_main_recovered, beta_main_shim_recovered)
  print(noquote(""))
  all_beta_functions(beta_2way_recovered, beta_2way_shim_recovered)
  print(noquote(""))
  all_beta_functions(beta_3way_recovered, beta_3way_shim_recovered)
  print(noquote(""))
  
  test_hierarchy_layer12(beta_main_shim_recovered, get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong = strong)
  print(noquote(""))
  test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                          get_psi_from_psi_vec3(beta_3way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong = strong)
  print(noquote(""))
  
  
  
}



#results_pipeline_shim(X=X, y=y, beta_main_init = beta_main_lasso, beta_2way_init = beta_2way_lasso, beta_3way_init = beta_3way_lasso, 
#                      lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, l1=l1, l2=l2, l3=l3,
#                      beta_main = beta_main, beta_2way = beta_2way, beta_3way = beta_3way , beta_main_recovered = beta_main_recovered, 
#                      beta_2way_recovered = beta_2way_recovered, beta_3way_recovered = beta_3way_recovered, tol=1e-3, threshold = 1)





#use it on y_centered
results_pipeline_pls_GLM<-function(X, y, n_comp, l1, l2, l3,
                               beta_main, beta_2way, beta_3way, beta_main_recovered, 
                               beta_2way_recovered, beta_3way_recovered, threshold = 0, strong = TRUE){
  
  
  
  range_main<-c(1: (l1+l2+l3) )
  range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
  range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
  
  
  pls_model <- plsr(g.link(y) ~ X, ncomp =  n_comp, scale = FALSE)
  
  print(coefficients((pls_model)))
  coefs_pls<-coefficients(pls_model)[-1]
  interc_pls<-coefficients(pls_model)[1]
  beta_main_pls<-coefs_pls[range_main]
  beta_2way_pls<-coefs_pls[range_theta]
  beta_3way_pls<-coefs_pls[range_psi]
  beta_2way_pls_without_gamma<-get_beta_vec_2way(beta_main_pls,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
  beta_3way_pls_without_delta<- get_beta_vec_3way(beta_2way_pls, l1=l1, l2=l2, l3=l3, only_beta = TRUE)
  
  predict_pls<- kappa1(predict(pls_model, newdata = X)[,,n_comp])
  cat("R2 pls: ", r2(y, predict_pls))
  plot(predict_pls, y)
  print(noquote(""))
  
  
  ######################### RESULTS pls #########################################################
  print("------------------ results pls without recovered -------------------")
  all_beta_functions(beta_main, beta_main_pls)
  print(noquote(""))
  all_beta_functions(beta_2way, beta_2way_pls)
  print(noquote(""))
  all_beta_functions(beta_3way, beta_3way_pls)
  print(noquote(""))
  
  
  ##hierarchy tests
  
  beta_2way_pls_matrix<-get_theta_from_theta_vec_2way3(beta_2way_pls,l1=l1,l2=l2, l3=l3)
  beta_3way_pls_table<-get_psi_from_psi_vec3(beta_3way_pls,l1=l1,l2=l2, l3=l3)
  
  test_hierarchy_layer12(beta_main_pls,beta_2way_pls_matrix, strong = strong)
  print(noquote(""))
  test_hierarchy_layer23(beta_2way_pls_matrix, beta_3way_pls_table, strong = strong)
  print(noquote(""))
  
  
  ##### RESULTS PLS ON RECOVERED PARAMS ##########
  
  beta_main_pls_recovered<- get_all_beta(beta_main_pls, l1=l1, l2=l2, l3=l3, threshold = threshold)
  beta_2way_pls_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_pls_matrix, l1=l1, l2=l2, l3=l3, threshold = threshold), l1=l1+1, l2=l2+1, l3=l3+1)
  beta_3way_pls_recovered<- get_psi_vec3( get_all_psi(beta_3way_pls_table, l1=l1, l2=l2, l3=l3, threshold = threshold) , l1=l1+1, l2=l2+1, l3=l3+1)
  
  
  #beta_main_pls
  #beta_main_pls_recovered
  
  print("--------------------- results pls recovered ------------------------")
  all_beta_functions(beta_main_recovered, beta_main_pls_recovered)
  print(noquote(""))
  all_beta_functions(beta_2way_recovered, beta_2way_pls_recovered)
  print(noquote(""))
  all_beta_functions(beta_3way_recovered, beta_3way_pls_recovered)
  print(noquote(""))
  
  
  test_hierarchy_layer12(beta_main_pls_recovered, get_theta_from_theta_vec_2way3( beta_2way_pls_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong = strong)
  print(noquote(""))
  test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_pls_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                          get_psi_from_psi_vec3(beta_3way_pls_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong = strong)
  
  return (list("intercept" = interc_pls, "main"= beta_main_pls , '2way' = beta_2way_pls, '3way' = beta_3way_pls))
  
  
}


#n_comp=6
#results_pipeline_pls(X=X, y=y_centered, n_comp=n_comp, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
#                     beta_3way_recovered, threshold = 0, strong = TRUE)








#### TEST SHIM CLASS ##########
  
  
#data<- create_basic_dataset()
#X<- data$X
#y<- data$y$obs

#beta_true<- data$beta[-1,]
#l1=8
#l2=8
#l3=4
#print(beta_true)

#range_main<-c(1: (l1+l2+l3) )
#range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
#range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

#beta_main<-beta_true[1:(l1+l2+l3)]
#beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
#beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
#beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
#beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


#gamma_true<-beta_2way/beta_2way_without_gamma
#gamma_true[is.nan(gamma_true)]<-0
#delta_true<-beta_3way/beta_3way_without_gamma
#delta_true[is.nan(delta_true)]<-0

#y_centered<-y-mean(y)
#lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=0.1)
#coefs_lasso<-coefficients(lasso_model)[-1]
#beta_main_lasso<-coefs_lasso[range_main]
#beta_2way_lasso<-coefs_lasso[range_theta]
#beta_3way_lasso<-coefs_lasso[range_psi]
#beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
#beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)


#length(c(beta_main, beta_2way, beta_3way))
#print("lasso")
#all_beta_functions(beta_main, beta_main_lasso)
#all_beta_functions(beta_2way, beta_2way_lasso)
#all_beta_functions(beta_3way, beta_3way_lasso)



#beta_hat<-beta_main_lasso
#gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
#gamma_hat[is.nan(gamma_hat)]<-0
#gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
#delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
#delta_hat[!is.finite(delta_hat)]<-0
#delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL

#lambda_beta<-100
#lambda_gamma<-1200
#lambda_delta<-500


#my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
#fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    #lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=6e-3)

#my_shim$R2_score(self=fitted, X_new=X, y_true=y )


#beta_all_shim<-fitted$beta_all
#beta_main_shim<-beta_all_shim[range_main]
#beta_2way_shim<-beta_all_shim[range_theta]
#beta_3way_shim<-beta_all_shim[range_psi]



#all_beta_functions(beta_main, beta_main_shim)
#all_beta_functions(beta_2way, beta_2way_shim)
#all_beta_functions(beta_3way, beta_3way_shim)





##### results pipeline PLS ######

#use it on y_centered
results_pipeline_pls<-function(X, y, n_comp, l1, l2, l3,
                               beta_main, beta_2way, beta_3way, beta_main_recovered, 
                               beta_2way_recovered, beta_3way_recovered, threshold = 0, strong = TRUE){
  
  
  
  range_main<-c(1: (l1+l2+l3) )
  range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
  range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
  
  
  pls_model <- plsr(y ~ X - 1, ncomp =  n_comp, scale = FALSE)
  
  
  coefs_pls<-coefficients(pls_model)
  beta_main_pls<-coefs_pls[range_main]
  beta_2way_pls<-coefs_pls[range_theta]
  beta_3way_pls<-coefs_pls[range_psi]
  beta_2way_pls_without_gamma<-get_beta_vec_2way(beta_main_pls,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
  beta_3way_pls_without_delta<- get_beta_vec_3way(beta_2way_pls, l1=l1, l2=l2, l3=l3, only_beta = TRUE)
  
  predict_pls<- predict(pls_model, newdata = X)[,,n_comp]
  cat("R2 pls: ", r2(y, predict_pls))
  plot(predict_pls, y)
  print(noquote(""))
  
  
  ######################### RESULTS pls #########################################################
  print("------------------ results pls without recovered -------------------")
  all_beta_functions(beta_main, beta_main_pls)
  print(noquote(""))
  all_beta_functions(beta_2way, beta_2way_pls)
  print(noquote(""))
  all_beta_functions(beta_3way, beta_3way_pls)
  print(noquote(""))
  
  
  ##hierarchy tests
  
  beta_2way_pls_matrix<-get_theta_from_theta_vec_2way3(beta_2way_pls,l1=l1,l2=l2, l3=l3)
  beta_3way_pls_table<-get_psi_from_psi_vec3(beta_3way_pls,l1=l1,l2=l2, l3=l3)
  
  test_hierarchy_layer12(beta_main_pls,beta_2way_pls_matrix, strong = strong)
  print(noquote(""))
  test_hierarchy_layer23(beta_2way_pls_matrix, beta_3way_pls_table, strong = strong)
  print(noquote(""))
  
  
  ##### RESULTS PLS ON RECOVERED PARAMS ##########
  
  beta_main_pls_recovered<- get_all_beta(beta_main_pls, l1=l1, l2=l2, l3=l3, threshold = threshold)
  beta_2way_pls_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_pls_matrix, l1=l1, l2=l2, l3=l3, threshold = threshold), l1=l1+1, l2=l2+1, l3=l3+1)
  beta_3way_pls_recovered<- get_psi_vec3( get_all_psi(beta_3way_pls_table, l1=l1, l2=l2, l3=l3, threshold = threshold) , l1=l1+1, l2=l2+1, l3=l3+1)
  
  
  #beta_main_pls
  #beta_main_pls_recovered
  
  print("--------------------- results pls recovered ------------------------")
  all_beta_functions(beta_main_recovered, beta_main_pls_recovered)
  print(noquote(""))
  all_beta_functions(beta_2way_recovered, beta_2way_pls_recovered)
  print(noquote(""))
  all_beta_functions(beta_3way_recovered, beta_3way_pls_recovered)
  print(noquote(""))
  
  
  test_hierarchy_layer12(beta_main_pls_recovered, get_theta_from_theta_vec_2way3( beta_2way_pls_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong = strong)
  print(noquote(""))
  test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_pls_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                          get_psi_from_psi_vec3(beta_3way_pls_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong = strong)
  
  return (list("main"= beta_main_pls , '2way' = beta_2way_pls, '3way' = beta_3way_pls))
  
  
}


#n_comp=6
#results_pipeline_pls(X=X, y=y_centered, n_comp=n_comp, l1=l1, l2=l2, l3=l3, beta_main, beta_2way, beta_3way, beta_main_recovered, beta_2way_recovered, 
#                     beta_3way_recovered, threshold = 0, strong = TRUE)



