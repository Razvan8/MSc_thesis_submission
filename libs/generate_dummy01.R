set.seed(123)

assert <- function(condition, message) {
  if (!condition)
    stop(message)
}

library(caret)




generate_data01 <-
  function(nlevels,
           min_main = 1,
           max_main = 10,
           min_2way = 1,
           max_2way = 10,
           min_3way = 1,
           max_3way = 10, 
           intercept=0, error_sd=0.3)
  {
    l1 <- nlevels[1]
    l2 <- nlevels[2]
    l3 <- nlevels[3]
    
    #assert(length(coefficients_main)==l1+l2+l3, 'length coef main is wrong')
    #assert(length(coefficients_2way)==l1*l2+l1*l3+l2*l3, 'length coef main is wrong')
    #assert(length(coefficients_3way)==l1*l2*l3, 'length coef main is wrong')
    
    range1 <- 1:l1
    range2 <- 1:l2
    range3 <- 1:l3
    
    
    # Generate the Cartesian product of the ranges
    combinations <- expand.grid(range1, range2, range3)
    combinations <-
      combinations[order(combinations[, 1], combinations[, 2], combinations[, 3]), ]
    
    
    # Convert the data frame to a matrix
    matrix <- as.matrix(combinations)
    #matrix <- as.data.frame(matrix)
    
    
    matrix <- apply(matrix, 2, factor)
    
    
    # Print the column names of the data frame
    print(colnames(matrix))
    
    # Generate dummy variables
    dummy <- dummyVars(" ~ .", data = matrix)
    dummy_matrix <- predict(dummy, newdata = matrix)
    interaction_terms <-
      data.frame(matrix(ncol = 0, nrow = nrow(dummy_matrix)))
    
    # Create interaction terms and append them to the dummy matrix
    for (i in 1:l1) {
      #A
      for (j in (l1 + 1):(l1 + l2 + l3)) {
        # AB or AC
        if (j < l1 + l2 + 1) {
          new_col_name <- paste0("A.", i, ":B.", j - l1)
        }
        else
        {
          new_col_name <- paste0("A.", i, ":C.", j - l1 - l2)
        }
        interaction_terms[[new_col_name]] <-
          dummy_matrix[, i] * dummy_matrix[, j]
      }
    }
    
    
    
    # Create interaction terms and append them to the dummy matrix
    for (i in (l1 + 1):(l1 + l2)) {
      #B
      for (j in (l1 + l2 + 1):(l1 + l2 + l3)) {
        # BC
        new_col_name <- paste0("B.", i - l1, ":C.", j - l1 - l2)
        interaction_terms[[new_col_name]] <-
          dummy_matrix[, i] * dummy_matrix[, j]
      }
    }
    
    
    ###ADD 3way
    for (i in 1:l1) {
      for (j in (l1 + 1):(l1 + l2)) {
        for (k in (l1 + l2 + 1):(l1 + l2 + l3)) {
          new_col_name <- paste0("A.", i, ":B.", j - l1, ":C.", k - l1 - l2)
          interaction_terms[[new_col_name]] <-
            dummy_matrix[, i] * dummy_matrix[, j] * dummy_matrix[, k]
        }
      }
    }
    
    
    
    # Bind the interaction terms to the original dummy matrix
    dummy_matrix <- cbind(dummy_matrix, interaction_terms)
    
    
    # Good colnames
    new_colnames <- c(paste0("A.", 1:l1),
                      paste0("B.", 1:l2),
                      paste0("C.", 1:l3))
    
    # Assign new column names to the dummy matrix
    colnames(dummy_matrix)[1:(l1 + l2 + l3)] <- new_colnames
    
    # Print the updated dummy matrix
    #print("Updated dummy matrix with interaction terms:")
    #print(dummy_matrix)
    
    mains <- runif(l1+l2+l3, min=min_main, max=max_main)
    twoways<- runif(l1*l2+l1*l3+l2*l3, min=min_2way, max=max_2way)
    threeways<- runif(l1*l2*l3, min=min_3way, max=max_3way)
    coefs<-c(mains, twoways, threeways)
    coefs<-matrix(coefs, ncol=1)
    coefs<- sample(c(-1, 1), size = nrow(coefs), replace = TRUE)*coefs
    rownames(coefs)<-colnames(dummy_matrix)
    #print(rownames(coefs))
    
    #Set rest of the coefs to 0
    levs.true <- c("A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",  "A.1:B.1","A.1:B.2",
                   "A.2:B.1","A.2:B.2", "A.3:B.1","A.3:B.2",
                   "A.1:C.1","A.1:C.2","A.2:C.1", "A.2:C.2","A.3:C.1","A.3:C.2")
                  ## 3way interactions can also be added
    
    
    coefs[which(!is.element(rownames(coefs),levs.true))] <- 0

    
    y<-as.matrix(dummy_matrix)%*%coefs+rnorm(dim(dummy_matrix)[1], sd=error_sd)
    #print(coefs)
    return(list("X"=dummy_matrix, "y"=y, "beta"=coefs))

  }
X<-generate_data01(c(6, 5, 4))$X

eigs<-svd(X)$d
max(eigs)/min(eigs)
eigs


