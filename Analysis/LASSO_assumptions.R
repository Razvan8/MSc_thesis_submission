dummy.matrix <- function(NF = 2, NL = rep(2, NF)) {
  # Computes dummy matrix from number of factors NF and number of levels NL
  # NF is an integer between 2 and 4
  # NL is a vector of length NF having integer entries between 2 and 100
  
  # Factors and Levels
  fac.tags <- c("A", "B", "C", "D")
  fac.levs <- as.character(1:100)
  
  if (NF == 2) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    
    # Dummy design matrix 2-way
    n.2w <- L.12
    p.2w <- n.2w - 1
    x.2w <- data.frame(matrix(0, nrow = n.2w, ncol = p.2w))
    rownames(x.2w) <- c(fac.12)
    colnames(x.2w) <- c(fac.1[-L.1], fac.2[-L.2],
                        sort(as.vector(outer(
                          fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
                        ))))
    for (col in 1:ncol(x.2w)) {
      col.tags <- unlist(strsplit(colnames(x.2w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.2w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[1]]), col] <- 1
        x.2w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[2]]), col] <- 1
        x.2w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[1]]), col] <- -1
        x.2w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[2]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.2w) == col.tags[1])
        col.2 <- which(colnames(x.2w) == col.tags[2])
        x.2w[, col] <- x.2w[, col.1] * x.2w[, col.2]
      }
    }
    
    return(x.2w)
  }
  
  if (NF == 3) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    L.13 <- L.1 * L.3
    L.23 <- L.2 * L.3
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep = ":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep = ":")))
    
    # Three-ways
    L.123 <- L.1 * L.2 * L.3
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep = ":")))
    
    # Dummy design matrix 3-way
    n.3w <- L.123
    p.3w <- n.3w - 1
    x.3w <- data.frame(matrix(0, nrow = n.3w, ncol = p.3w))
    rownames(x.3w) <- c(fac.123)
    colnames(x.3w) <- c(
      fac.1[-L.1],
      fac.2[-L.2],
      fac.3[-L.3],
      sort(as.vector(outer(
        fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.1[-L.1], fac.3[-L.3], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.1[-L.1], sort(as.vector(outer(
          fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
        ))), paste, sep = ":"
      )))
    )
    for (col in 1:ncol(x.3w)) {
      col.tags <- unlist(strsplit(colnames(x.3w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[1]]), col] <- 1
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[2]]), col] <- 1
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[3]]), col] <- 1
        x.3w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[1]]), col] <- -1
        x.3w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[2]]), col] <- -1
        x.3w[grepl(paste(fac.tag[1], L.3, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[3]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.3w) == col.tags[1])
        col.2 <- which(colnames(x.3w) == col.tags[2])
        x.3w[, col] <- x.3w[, col.1] * x.3w[, col.2]
      }
      if (length(col.tags) == 3) {
        col.1 <- which(colnames(x.3w) == col.tags[1])
        col.2 <- which(colnames(x.3w) == col.tags[2])
        col.3 <- which(colnames(x.3w) == col.tags[3])
        x.3w[, col] <- x.3w[, col.1] * x.3w[, col.2] * x.3w[, col.3]
      }
    }
    
    return(x.3w)
  }
  
  if (NF == 4) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    L.4 <- NL[4]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep = ".")
    fac.4 <- paste(fac.tags[4], fac.levs[1:L.4], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    L.13 <- L.1 * L.3
    L.14 <- L.1 * L.4
    L.23 <- L.2 * L.3
    L.24 <- L.2 * L.4
    L.34 <- L.3 * L.4
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep = ":")))
    fac.14 <- sort(as.vector(outer(fac.1, fac.4, paste, sep = ":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep = ":")))
    fac.24 <- sort(as.vector(outer(fac.2, fac.4, paste, sep = ":")))
    fac.34 <- sort(as.vector(outer(fac.3, fac.4, paste, sep = ":")))
    
    # Three-ways
    L.123 <- L.1 * L.2 * L.3
    L.124 <- L.1 * L.2 * L.4
    L.134 <- L.1 * L.3 * L.4
    L.234 <- L.2 * L.3 * L.4
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep = ":")))
    fac.124 <- sort(as.vector(outer(fac.1, fac.24, paste, sep = ":")))
    fac.134 <- sort(as.vector(outer(fac.1, fac.34, paste, sep = ":")))
    fac.234 <- sort(as.vector(outer(fac.2, fac.34, paste, sep = ":")))
    
    # Four-ways
    L.1234 <- L.1 * L.2 * L.3 * L.4
    fac.1234 <-
      sort(as.vector(outer(fac.1, fac.234, paste, sep = ":")))
    
    # Dummy design matrix 4-way
    n.4w <- L.1234
    p.4w <- n.4w - 1
    x.4w <- data.frame(matrix(0, nrow = n.4w, ncol = p.4w))
    rownames(x.4w) <- c(fac.1234)
    colnames(x.4w) <-
      c(
        fac.1[-L.1],
        fac.2[-L.2],
        fac.3[-L.3],
        fac.4[-L.4],
        sort(as.vector(outer(
          fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], fac.3[-L.3], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], sort(as.vector(outer(
            fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], sort(as.vector(outer(
              fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
            ))),
            paste, sep =
              ":"
          ))), paste, sep = ":"
        )))
      )
    for (col in 1:ncol(x.4w)) {
      col.tags <- unlist(strsplit(colnames(x.4w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[1]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[2]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[3]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[4]]), col] <- 1
        x.4w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[1]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[2]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.3, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[3]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.4, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[4]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        x.4w[, col] <- x.4w[, col.1] * x.4w[, col.2]
      }
      if (length(col.tags) == 3) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        col.3 <- which(colnames(x.4w) == col.tags[3])
        x.4w[, col] <- x.4w[, col.1] * x.4w[, col.2] * x.4w[, col.3]
      }
      if (length(col.tags) == 4) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        col.3 <- which(colnames(x.4w) == col.tags[3])
        col.4 <- which(colnames(x.4w) == col.tags[4])
        x.4w[, col] <-
          x.4w[, col.1] * x.4w[, col.2] * x.4w[, col.3] * x.4w[, col.4]
      }
    }
    
    return(x.4w)
  }
  
}


X_df<-dummy.matrix(NF = 2, NL = c(4,2))
X_df
class(X_df)
X<-as.matrix(X_df)
class(X)
t(X)%*%X
eigen((t(X)%*%X))




# Load the nloptr package
library(nloptr)

# Define the objective function
objective_function <- function(params, X, S) {
  # Extract beta_S and beta_Sc from the params vector
  p_S <- length(S)
  beta_S <- params[S]
  beta_Sc <- params[-S]
  
  # Split X into X_S and X_Sc based on the indices S and S^c
  X_S <- X[, S]
  X_Sc <- X[, -S, drop = FALSE]
  
  # Compute the quadratic term ||X_S * beta_S - X_Sc * beta_Sc||_2^2
  diff <- X_S %*% beta_S - X_Sc %*% beta_Sc
  return(sum(diff^2))
}

# Equality constraint for ||beta_S||_1 = 1
equality_constraint <- function(params, X, S) {
  p_S <- length(S)
  beta_S <- params[S]
  return(sum(abs(beta_S)) - 1)  # This should equal 0 (for equality constraint)
}

# Inequality constraint for ||beta_Sc||_1 <= 4 (L = 4 is hardcoded)
inequality_constraint <- function(params, X, S) {
  p_S <- length(S)
  beta_Sc <- params[-S]
  return(sum(abs(beta_Sc))-4)  # This should be >= 0 (for inequality constraint)
}

# Main function to solve the optimization problem using nloptr
minimize_l1_constrained_nloptr <- function(X, S) {
  # Initial guess for beta_S and beta_Sc (random initialization)
  p_S <- length(S)
  p_Sc <- ncol(X) - p_S
  init_guess <- runif(p_S + p_Sc, -0.01, 0.01)
  
  # Lower and upper bounds for the variables
  lower_bounds <- rep(-4, p_S + p_Sc)
  upper_bounds <- rep(4, p_S + p_Sc)
  
  # Define the optimization problem using nloptr with ISRES
  result <- nloptr(
    x0 = init_guess,                            # Initial guess
    eval_f = objective_function,                # Objective function
    lb = lower_bounds,                          # Lower bounds
    ub = upper_bounds,                          # Upper bounds
    eval_g_eq = equality_constraint,            # Equality constraint function
    eval_g_ineq = inequality_constraint,        # Inequality constraint function
    opts = list("algorithm" = "NLOPT_GN_ISRES", # Use ISRES algorithm
                "xtol_rel" = 1e-4,              # Tolerance for convergence
                "maxeval" = 10000),              # Maximum number of iterations
    X = X, S = S                                # Pass the matrix and set of indices
  )
  
  # Extract the optimal solution
  beta_S_opt <- result$solution[S]
  beta_Sc_opt <- result$solution[-S]
  objective_value <- result$objective
  
  # Return the result
  return(list(
    beta_S = beta_S_opt,
    beta_Sc = beta_Sc_opt,
    objective_value = objective_value,
    status = result$status,
    message = result$message
  ))
}

# Example usage:
set.seed(123)
X <- matrix(rnorm(400), nrow = 20)  # A random matrix with 10 rows, 10 columns
S <- c(1, 3, 5)  # Example of indices in S

X_df<-dummy.matrix(NF = 3, NL = c(3,3,3))
X<-as.matrix(X_df)
# Call the function to solve the problem
#result <- minimize_l1_constrained_nloptr(X, S)

# Print the result
#print(result)
sum(abs(result$beta_Sc))
res<-eigen(t(X)%*%X)
# Step 3: Extract eigenvectors
eigenvectors <- res$vectors


# Step 4: Calculate and display the norm of each eigenvector
for (i in 1:ncol(eigenvectors)) {
  eigvec <- eigenvectors[, i]
  eigvec_norm2 <- norm(eigvec, type = "2")  # Type "2" for the Euclidean norm
  eigvec_norm1 <- sum(abs(eigvec)) # Type "1" for the Euclidean norm
  cat("Norm1 of eigenvector", i, ":", eigvec_norm1, "\n")
  #cat("Norm2 of eigenvector", i, ":", eigvec_norm2, "\n")
}





for (i in c(3:6)){
  for (j in c(3:6)){
X_df<-dummy.matrix(NF = 3, NL = c(i,i,j))
X<-as.matrix(X_df)
res<-eigen(t(X)%*%X)
#length(res$values)
vect<-res$vectors[,length((res$values))]
norm1<-sum(abs(vect))
norm2<-sqrt(sum(vect**2))
print(norm1/norm2)}}
