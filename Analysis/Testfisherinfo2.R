libs_path<-file.path("..","libs")
source(file.path(libs_path,'Shim3way_GLM.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))






set.seed(1)
create_fisher_dataset <- function(l1,l2,l3) {
  set.seed(1)
  x.3w <- dummy.matrix(NF = 3, NL = c(l1+1, l2+1, l3+1))
  print(colnames(x.3w))
  
  col_mains <- colnames(x.3w)[c(1:(l1 + l2 + l3))]
  #print(col_mains)
  col_theta_good <- c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))
    }
  }
  for (j in c(1:l2))
  {
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))
    }
  }
  
  
  col_psi_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      for (k in c(1:l3))
      {col_psi_good <- c(col_psi_good, paste0("A.", i, ":B.", j, ":C.", k))
      }}}
  
  
  
  
  
  #print(col_theta_good)
  col_all_good = c(col_mains, col_theta_good, col_psi_good)
  print(col_all_good)
  x.3w<-x.3w[,col_all_good]
  print(x.3w)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 2
  beta.max <- 5
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(3)
  
  
  
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max) * sample(c(1, -1), size =  n.3w, replace = TRUE)
  
  
  levs.true <-
    c(
      "A.1","A.2","A.3","A.4", "A.5", "A.10", "A.11", "A.12",
      "B.1","B.2","B.3",
      "C.1","C.2","C.3",
      "A.1:B.1","A.1:B.2","A.1:B.3", "A.1:C.1", "A.1:C.2", "A.1:C.3",
      "A.2:B.1","A.2:B.2", "A.2:C.1", "A.2:C.2", "A.2:C.3", 
      "A.10:B.1","A.10:B.2","A.10:B.3", "A.10:C.1", "A.10:C.2", "A.10:C.3",
      "A.11:B.1","A.11:B.2","A.11:B.3", "A.11:C.1", "A.11:C.2", "A.11:C.3",
      "B.1:C.1", "B.1:C.2",
      
      "A.1:B.1:C.1",
      "A.2:B.1:C.1",
      "A.10:B.1:C.1",
      "A.11:B.1:C.1"
    )
  
  
  beta.true$coeffs[which(!is.element(rownames(beta.true), levs.true))] <-0
  beta.true$coeffs[2:(2+l1+l2+l3)]<-beta.true$coeffs[2:(2+l1+l2+l3)]+sign(beta.true$coeffs[2:(2+l1+l2+l3)])*3 #mains are from U[4,8]
  beta.true$coeffs[1]<- -8
  #beta.true_new <-
  #beta.true[c("interc", col_all_good), , drop = FALSE]
  #print("new beta true")
  #print(beta.true_new)
  #rownames(beta.true_new) <- c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 0.5
  y.3w <- data.frame(row.names = rownames(x.3w))
  set.seed(111)
  print(mean(
    beta.true$coeffs[1] + as.matrix(x.3w) %*% as.vector(beta.true$coeffs)[-1]
  ))
  eta=beta.true$coeffs[1] + as.matrix(x.3w) %*% as.vector(beta.true$coeffs)[-1]
  prob=exp(eta)/(exp(eta)+1)
  y.3w$obs <- sample_contbern(p=prob)
  y.3w$true <-kappa1(beta.true$coeffs[1] + as.matrix(x.3w) %*% as.vector(beta.true$coeffs)[-1]) #get bern
  x.3w_new <- x.3w[, col_all_good , drop = FALSE]
  
  if (all(rownames(beta.true)[-1] == colnames(x.3w_new)) == TRUE)
  {
    print("Data loaded properly")
  }
  print("SNR:")
  print(var(as.vector(y.3w$true)) / (sigma.y ^ 2))
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list(
    'X' = as.matrix(x.3w_new),
    'y' = y.3w,
    'beta' = beta.true
  ))
}









data<- create_fisher_dataset(l1=10,l2=1,l3=1)


rownames(X)<- data$X

vals<-eigen(t(X)%*%X)$values
print(sort(vals))

y<- data$y$obs
y<-array(y, dim=(c(length(y),1)) )

hist(y)




#################################################################################



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



beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

gamma_matrix<-get_theta_from_theta_vec_2way3(gamma_true,l1=l1,l2=l2, l3=l3)/2
delta_table<-get_psi_from_psi_vec3(delta_true,l1=l1, l2=l2, l3=l3)
X_mains<-X[,range_main]



library(numDeriv)

# Define the model's predicted values Y_hat given beta, gamma, delta
predict_model <- function ( beta, gamma, delta, X) {
  beta_2way <- get_beta_vec_2way(beta, l1 = l1, l2 = l2, l3 = l3, gamma = gamma, only_beta = FALSE) ###This is with delta
  beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta, only_beta = FALSE) #This is with gamma WITH delta
  beta_all<-array(c(self$beta_hat, beta_2way, beta_3way))
  Y_hat<-  X%*%beta_all
  return(Y_hat)
}

# Define the log-likelihood function
log_likelihood <- function(beta, gamma, delta, X, Y, sigma2 = 1) {
  Y_hat <- predict_model(beta, gamma, delta, X)  # Compute predicted values
  residuals <- Y - Y_hat
  
  # Normal likelihood for the residuals
  log_likelihood_value <- -0.5 * sum(residuals^2 / sigma2)
  
  return(log_likelihood_value)
}

predict_model(beta_main, gamma_matrix, delta_table, X)
dim(X_mains)
lentgh(beta_main)

# Compute the numerical Hessian matrix
hessian_matrix <- hessian(beta_main, gamma_matrix, delta_table, log_likelihood, initial_params, X = X_mains, Y = y)

# The Fisher Information Matrix (FIM) is the negative of the Hessian
fisher_info <- -hessian_matrix

