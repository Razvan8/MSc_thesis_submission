assert <- function(condition, message) {
  if (!condition)
    stop(message)
}



check_and_zero_diagonal <- function(matrix, strong = FALSE) {
  # Check if the matrix is square
  if (nrow(matrix) != ncol(matrix)) {
    stop("The matrix must be square.")
  }
  

  all_zero <- TRUE
  
  # Check the diagonal elements
  for (i in 1:nrow(matrix)) {
    if (matrix[i, i] != 0) {
      all_zero <- FALSE
      break
    }
  }
  

  if (!all_zero) {
    print(
      "Not all diagonal elements are zero. Setting them to zero. This might be a sign of numerical instability."
    )
  }
  
  # Set all diagonal elements to zero
  for (i in 1:nrow(matrix)) {
    matrix[i, i] <- 0
  }
  
  return(matrix)
}



test_hierarchy_layer12 <- function(beta, theta, strong = FALSE)
{
  hierarchy_break <- 0
  theta <- check_and_zero_diagonal(theta)
  p <- length(c(beta))
  assert(all(dim(theta) == c(p, p)), "The dimensions are not ok in the hierarchy test")
  
  zeros_positions <- which(beta == 0)
  #print(zeros_positions)
  
  if (strong == TRUE)
  {
    substract <- 0# substract because those with 2 zeros appear twice
    for (i in c(1:p))
    {
      for (j in zeros_positions)
      {
        if (theta[i, j] + theta[j, i] != 0)
        {
          hierarchy_break <- hierarchy_break + 1
          if (i %in% zeros_positions)
          {
            substract <- substract + 1/2
          }
        }
      }
    }
    
  hierarchy_break <-hierarchy_break - substract
    } ###combinations of 0s appear twice
  
  
  
  
  if (strong == FALSE)
  {
    for (i in zeros_positions)
    {
      for (j in zeros_positions)
      {
        if (theta[i, j] + theta[j, i] != 0)
        {
          hierarchy_break <- hierarchy_break + 1 / 2
        } ##because it appears 2 times
      }
    }
  }
  
  cat("Hierarchy is broken for", hierarchy_break, "times")
  
}

#beta <- c(1, 1, 2, 1)

#print("matrix")
#matrix <- matrix(4, nrow = 4, ncol = 4)
#matrix[1, 2] <- 1
#matrix[2, 1] <- 0

#print(matrix)

#test_hierarchy_layer12(beta = beta, theta = matrix, strong = TRUE)

test_hierarchy_layer23 <- function(theta, psi, strong = FALSE)
{
  hierarchy_break <- 0
  theta <- check_and_zero_diagonal(theta)
  p <- dim(theta)[1]
  assert(all(dim(theta) == c(p, p)), "The dimensions are not ok in the hierarchy test")
  assert(all(dim(psi) == c(p, p, p)), "The dimensions are not ok in the hierarchy test")
  
  zeros_positions <- list()
  #print(zeros_positions)
  
  if (strong == FALSE)
  {
    # substract because those with 2 zeros appear twice
    for (i in c(1:(p-2)) )
    {
      for (j in c((i + 1):(p-1)) )
      {
        for (k in c((j + 1):p)) {
          if ((theta[i, j] + theta[j, i] == 0) &&
              (theta[i, k] + theta[k, i] == 0) &&
              (theta[j, k] + theta[k, j] == 0) &&
              psi[i, j, k] + psi[i, k, j] + psi[j, i, k] + psi[j, k, i] + psi[k, i, j] + psi[k, j, i] != 0)
       { hierarchy_break <- hierarchy_break + 1 }    

        }
            }
          }
        }
      
    
    
    
    
  if (strong == TRUE)
  {
    # substract because those with 2 zeros appear twice
    for (i in c(1:(p-2)) )
    {
      for (j in c((i + 1): (p-1)) )
      {
        for (k in c((j + 1):p)) {
          if ( ( (theta[i, j] + theta[j, i] == 0) |
              (theta[i, k] + theta[k, i] == 0) |
              (theta[j, k] + theta[k, j] == 0) ) &&
              (psi[i, j, k] + psi[i, k, j] + psi[j, i, k] + psi[j, k, i] + psi[k, i, j] + psi[k, j, i] != 0) )
          { hierarchy_break <- hierarchy_break + 1 }    
          
        }
      }
    }
  }
  
    
    cat("Hierarchy is broken for", hierarchy_break, "times")
    
  }
  


#print("matrix")
#matrix <- matrix(0, nrow = 4, ncol = 4)
#matrix[1, 2] <- 1
#matrix[2, 3] <- 1

#psi<-array(1, dim=c(4,4,4))
#psi[1,2,3]<--1

#print(matrix)

#test_hierarchy_layer23(matrix,psi, strong = FALSE)





###########TESTS  for APROX HIERARCHY ###################


test_hierarchy_layer12_approx <- function(beta, theta, strong = FALSE, threshold_coef=1e-3)
{
  
  threshold<-threshold_coef*mean(abs(theta[theta!=0]))
  hierarchy_break <- 0
  theta <- check_and_zero_diagonal(theta)
  p <- length(c(beta))
  assert(all(dim(theta) == c(p, p)), "The dimensions are not ok in the hierarchy test")
  
  zeros_positions <- which(beta == 0)
  #print(zeros_positions)
  
  if (strong == TRUE)
  {
    substract <- 0# substract because those with 2 zeros appear twice
    for (i in c(1:p))
    {
      for (j in zeros_positions)
      {
        if (abs(theta[i, j] + theta[j, i] ) >=threshold)
        {
          hierarchy_break <- hierarchy_break + 1
          if (i %in% zeros_positions)
          {
            substract <- substract + 1/2
          }
        }
      }
    }
    
    hierarchy_break <-hierarchy_break - substract
  } ###combinations of 0s appear twice
  
  
  
  
  if (strong == FALSE)
  {
    for (i in zeros_positions)
    {
      for (j in zeros_positions)
      {
        if ( abs(theta[i, j] + theta[j, i])  >=threshold )
        {
          hierarchy_break <- hierarchy_break + 1 / 2
        } ##because it appears 2 times
      }
    }
  }
  
  paste("Hierarchy is broken for", hierarchy_break, "times")
  
}

beta <- c(0, 0, 0, 0)

print("matrix")
matrix <- matrix(1, nrow = 4, ncol = 4)
matrix[1, 2] <- 0.00001
matrix[2, 1] <- 0.00001

#print(matrix)

test_hierarchy_layer12_approx(beta = beta, theta = matrix, strong = TRUE)

test_hierarchy_layer23_approx <- function(theta, psi, strong = FALSE, threshold_coef=1e-3)
{
  threshold<-mean(abs(psi[psi!=0]))*threshold_coef
  hierarchy_break <- 0
  theta <- check_and_zero_diagonal(theta)
  p <- dim(theta)[1]
  assert(all(dim(theta) == c(p, p)), "The dimensions are not ok in the hierarchy test")
  assert(all(dim(psi) == c(p, p, p)), "The dimensions are not ok in the hierarchy test")
  
  zeros_positions <- list()
  #print(zeros_positions)
  
  if (strong == FALSE)
  {
    # substract because those with 2 zeros appear twice
    for (i in c(1:(p-2)) )
    {
      for (j in c((i + 1):(p-1)) )
      {
        for (k in c((j + 1):p)) {
          if ((theta[i, j] + theta[j, i] == 0) &&
              (theta[i, k] + theta[k, i] == 0) &&
              (theta[j, k] + theta[k, j] == 0) &&
              abs(psi[i, j, k]) + abs(psi[i, k, j]) + abs(psi[j, i, k]) + abs(psi[j, k, i]) + abs(psi[k, i, j]) + abs(psi[k, j, i]) >=threshold)
          { hierarchy_break <- hierarchy_break + 1 }    
          
        }
      }
    }
  }
  
  
  
  
  
  if (strong == TRUE)
  {
    # substract because those with 2 zeros appear twice
    for (i in c(1:(p-2)) )
    {
      for (j in c((i + 1): (p-1)) )
      {
        for (k in c((j + 1):p)) {
          if ( ( (theta[i, j] + theta[j, i] == 0) |
                 (theta[i, k] + theta[k, i] == 0) |
                 (theta[j, k] + theta[k, j] == 0) ) &&
               ( abs(psi[i, j, k]) + abs(psi[i, k, j]) + abs(psi[j, i, k]) + abs(psi[j, k, i]) + abs(psi[k, i, j]) + abs(psi[k, j, i]) >=threshold) )
          { hierarchy_break <- hierarchy_break + 1 }    
          
        }
      }
    }
  }
  
  
  paste("Hierarchy is broken for", hierarchy_break, "times")
  
}



#print("matrix")
matrix <- matrix(0, nrow = 4, ncol = 4)
matrix[1, 2] <- 1
matrix[2, 3] <- 2

matrix[matrix>0]

psi<-array(1, dim=c(4,4,4))
psi[1,3,4]<--0.3e-4
psi[1,4,3]<--0.3e-4
psi[3,4,1]<--0.3e-4
psi[3,1,4]<--0.3e-3
psi[4,1,3]<--0.3e-4
psi[4,3,1]<--0.3e-3

#print(matrix)

test_hierarchy_layer23_approx(matrix,psi, strong = FALSE)





