### FUNCTIONS THAT I USE IN WEAK_HIERNET



## ASSERT FUNCTION
assert <- function(condition, message) {
  if (!condition) stop(message)
}


###USEFUL FUNCTION FOR LASSO LOSS AND WEAK HIER NET


stable_alpha<-function(x)
{if (x>=0)
{return((x+1e-5)*1.0001)}}


get_positions <- function(p, j) {
  # Generate all pairs
  all_pairs <- expand.grid(1:p, 1:p)
  #print(all_pairs)
  
  # Find positions where THE second column is equal to j
  positions <- which( all_pairs[, 2] == j)
  
  # Return the positions
  return(positions)
}



matrix_position_to_vector_index3<- function(position_tuple, l1=36,l2=3,l3=4) ## takes into account / works only for possible combinations!!!!
{
  

  x<-position_tuple[1]
  y<-position_tuple[2]
  
  
  
  range_x<-get_range3(x,l1=l1,l2=l2,l3=l3)
  range_y<-get_range3(y,l1=l1,l2=l2,l3=l3)
  
  #print(x)
  #print(y)

  assert(x<= l1+l2, "x should be <=l1+l2")
  assert(x<y, "x<y")
  assert(y>l1, 'y should be >l1+l2')
  #assert(any (get_range3(x,l1=l1,l2=l2,l3=l3)) != any(get_range3(y,l1=l1,l2=l2,l3=l3)), "x and y should be from different range")
  

   
  
  if( all(range_x == c(1:l1)) ==TRUE ) #ab or ac
  { 
   position_theta<- (x-1)*(l2+l3) +(y-l1)  }
    
    
  if( all ( range_x == c( (l1+1): (l1+l2) ) ) == TRUE )  #bc
  {position_theta<-l1*(l2+l3) + (x-l1-1)*l3 + y- (l1+l2)  } 
  #print(position_theta)
  return(position_theta)
  
}



get_cols_Xj3<-function(X,j,l1,l2,l3) ##assummes data in -1 0 1 form?
{ X_j<-matrix(0,nrow=nrow(X), ncol = l1+l2+l3)
  range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
  all_possible = c(range1, range2, range3)
  range_j<-get_range3(j,l1=l1,l2=l2,l3=l3)
  range_possible<-setdiff(all_possible,range_j)
  for (i in range_possible){
    #cat("i,j: ",i,j,'/n')
    X_j[,i]<-X[,matrix_position_to_vector_index3(sort(c(i,j)),l1=l1,l2=l2,l3=l3)] #coloanele care au main j in ele
  }
  return(X_j)
}



###SOFT THRESHOLDING OPERATOR
Soft_thresholding <- function(c, lambda) {
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  assert(lambda >= 0, "lambda cannot be negative.")
  
  # Apply soft thresholding component-wise
  result <- sign(c) * pmax(abs(c) - c(lambda), 0)
  
  return(result)
}




## RELU function
RELU <- function(x) {
  if (x>=0){
    return(x)}
  else {
    return (0)
  }
}




evaluate_knots_seq <- function(theta_tilda_j, lambda, f, t = 1,c=1) {
  p <- length(theta_tilda_j)
  
  #cat("thetatilda_jk", theta_tilda_j)
  
  set_elements <- 0##add also 0 to knots, it has positive value in evaluated
  
  for (k in 1:p)  { 
    if (abs(theta_tilda_j[k])>0){ ##discard impossible combinations
      element <- abs(theta_tilda_j[k]) / t - lambda / 2 #element to add in knot set
      set_elements <- c(set_elements, element)}
  }
  
  selected_elements <- set_elements[set_elements >= 0] ### cred ca e nevoie si de 0?!?!
  
  #selected_elements<-c(selected_elements,0) #### ADD O TO THE LIST ### IT SHOULD NOT INFLUENCE (TAKE CARE HERE!!!!!!!!!!)
  selected_elements<- sort(selected_elements) ## sort them increasing
  #print(selected_elements)
  
  result_vector <- sapply(selected_elements, f)## Evaluate f(p) for every p
  #cat("selected elements : ", selected_elements)
  #print("  ")
  #cat("result vector : ", result_vector)
  
  
  #print(list("knots" = selected_elements,"evaluated"=result_vector))
  return(list("knots" = selected_elements,"evaluated"=result_vector))
}





# Function to find adjacent knots and calculate alpha_hat - HELPING FUNCTION FOR ONEROW
find_adjacent_knots_and_alpha <- function(knots, evaluated) {
  # Find the index where evaluated > 0
  #print(evaluated)
  positive_index <- max(which(evaluated > 0))
  #print(which(evaluated > 0))
  #cat("positive index", positive_index)
  
  # Ensure there is a positive index and the next index exists
  if (!is.na(positive_index) && (positive_index + 1) <= length(knots)) {
    # Get corresponding knots
    p1 <- knots[positive_index]
    p2 <- knots[positive_index + 1]
    
    # Calculate alpha_hat
    #alpha_hat <- -evaluated[positive_index] * (evaluated[positive_index + 1] - evaluated[positive_index]) / (p2 - p1) ##initial one
    alpha_hat <- (p1*evaluated[positive_index + 1] - p2*evaluated[positive_index]) / (evaluated[positive_index + 1] - evaluated[positive_index]) #new
    #print(positive_index)
    return(alpha_hat)
  } else {
    # Return an informative message if no such adjacent knots are found
    stop("No adjacent knots found with evaluated > 0 and evaluated < 0 at consecutive positions.")
  }
}

final_return_seq_2way3<- function(theta_hat_j, lambda, t, alpha_hat) # returns theta_jk with all the elements inside (vector with 40 elements)
{ 
  alpha_hat<-stable_alpha(alpha_hat)
  theta_hat_j = Soft_thresholding(theta_hat_j, t*(lambda/2 + alpha_hat))
  return (theta_hat_j)
}

## FUNCTIONS ##########################################################

get_all_possible_kj3<-function(l1=36,l2=3,l3=4) ### all possible pair combinations (includes both i,j and j,i)
{range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
all_possible = c(range1, range2, range3)
possible_kj=list()

## ij 

for ( range in list(range1, range2, range3)  ){
  for (i in range){
    for (j in setdiff(all_possible,range) ){
      possible_kj <- c(possible_kj, list(c(i, j)))
    }}}
return(possible_kj)
}

#get_all_possible_kj3(l1=36,l2=3,l3=4)




get_possible_positions_2way3<- function(l1=36,l2=3,l3=4) #get possible positions only once and return a list of tuples of these positions
                                                         #This is the order of appereance in data
{
  
  
  list_possible_combinations <- list()
  range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
 
  ## case 1 :a with b or s
  for (i in range1)
  {for (j in c(range2,range3))
  {list_possible_combinations <- append(list_possible_combinations, list(c(i,j)))
  }}
  
  ## case 2: b with s
  for (i in range2)
  {for (j in range3)
  {list_possible_combinations <- append(list_possible_combinations, list(c(i,j)))
  }}
  
  
  return(list_possible_combinations)
}

#get_possible_positions_2way3()

#my_table <- array(1:64, dim = c(4, 4, 4))








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



# Create the 3-dimensional array
#array_data <- array(data = c(
# Values for the first layer
#0, 0.000000, 0.0000000, 0.000000,
#0, 0.000000, -3.3490767, 1.899415,
#0, 2.887088, 0.0000000, 3.100585,
#0, 1.206102, 0.1402117, 0.000000,

# Values for the second layer
#0.0000000, 0, 0.8724023, -0.07197379,
#0.0000000, 0, 0.0000000, 0.00000000,
#-3.7229584, 0, 0.0000000, 0.50213101,
#0.7168255, 0, 1.8879085, 0.00000000,

# Values for the third layer
#0.000000, 3.247535, 0, 0.3154444,
#2.306129, 0.000000, 0, 1.5593346,
#0.000000, 0.000000, 0, 0.0000000,
#-1.263373, -1.752465, 0, 0.0000000,

# Values for the fourth layer
#0.0000000, 1.4078473, 1.676665, 0,
#0.3514781, 0.0000000, 1.668121, 0,
#3.0134272, 0.1127374, 0.000000, 0,
#0.0000000, 0.0000000, 0.000000, 0
#), dim = c(4, 4, 4))

# Print the array
#print(array_data)


#theta<-  array(1, dim = c(4, 4, 4))
# 123 45 67 89
# Set some positions to 1
#theta[1, 2, 3] <- 6
#theta[1, 3, 2] <- 6
#theta[3, 1, 2] <- 6













sets_0s_theta_2way3<-function(theta,l1=21,l2=14,l3=2)### SET 0 in impossible positions
  
{range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))

## ij same range
for ( range in list(range1, range2, range3) ){
  for (i in  range) {
    for (j in range)  {   
      theta[i,j]<-0
      theta[j,i]<-0
    }}}

return(theta)
}



get_range3<- function(x,l1=36,l2=3,l3=4)
{if (x<=l1)
{return(c(1:l1))}
  if (x<=l1+l2)
  {return(c( (l1+1) : (l1+l2) ))}

  return(c( (l1+l2+1) : (l1+l2+l3) ))

}






#get_cols_Xkj(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7), nrow=4, ncol =7),3,4, 2,1,1,1)

#print(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7),nrow=4, ncol=7))


###############################################################################################################################################




# function for ONEROW
ONEROW_SEQ12 <- function(theta_tilda_j, lambda, beta_bound_j, j, t=1, c=1) {
  
  f <- function(alpha) {
    #print('alpha')
    #print(alpha)
    #cat(' suma ',sum(abs(Soft_thresholding( Theta_tilda_j, t*(lambda/2+alpha)  ) ) ))
    alpha<-stable_alpha(alpha)
    result<- sum( abs( Soft_thresholding( theta_tilda_j[-j], t*(lambda/2+alpha)  ) ) ) - c*abs(beta_bound_j) ### 1
    if (abs(result)<1e-10) ###this is new
    { if(abs(result)>0)
    {print("There might be numerical intsability")}
      return(0)}
    ## to solve numerical instability
    return(result)
  }
  if ( f(0)<=0)  ### 1 a)
  {#print('a')
    alpha_hat<-0
    return(final_return_seq_2way3(theta_hat_j = theta_tilda_j, lambda = lambda, t=t, alpha_hat = alpha_hat) ) ### 2)
  }   
  
  
  ### 1b) ### 1c)
  knots_evaluated <- evaluate_knots_seq( theta_tilda_j=theta_tilda_j, lambda=lambda, f=f, t = t,c=c)
  ### 1d)
  knots=knots_evaluated$knots
  evaluated=knots_evaluated$evaluated
  #print(knots_evaluated)
  
  zero_indices <- which(evaluated == 0)
  
  if (length(zero_indices) > 0) {
    #print('d')
    alpha_hat <- knots[zero_indices[1]] # i.e. alpha_hat =p s.t. f(p)=0
    #cat("alpha hat",alpha_hat)
    return (final_return_seq_2way3(theta_hat_j = theta_tilda_j, lambda =lambda, t=t, alpha_hat = alpha_hat) )}   ### 2)
  
  
  alpha_hat <- find_adjacent_knots_and_alpha(knots, evaluated)    ### 1e)
  
  #print('e')
  #cat("alpha hat",alpha_hat)
  
  ### STEP 2 
  return (final_return_seq_2way3(theta_hat_j = theta_tilda_j, lambda =  lambda, t=t, alpha_hat = alpha_hat)) ### 2)
  
}

#theta_tilda_j=c(3,2,0,0,0)
#lambda=6
#beta_bound_j=2
#t=1

#ONEROW_SEQ12(theta_tilda_j=theta_tilda_j, lambda=lambda, beta_bound_j=beta_bound_j,  t=t)



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





######################### WEAKHIERNET CLASS #########################################

# Define the WeakHierNet class
HierNet_seq_2way3 <- function(X, theta_init, y, beta_bound, lambda, t=1, tol=1e-6, max_iter=5000, eps=1e-8,l1=36,l2=3,l3=4, scale=TRUE) {
  
  
  
  self = list()
  self$theta_hat <- theta_init #matrix form
  self$beta_bound <- beta_bound
  self$vec_theta_hat<-get_theta_vec_2way3(theta_init,l1=l1,l2=l2,l3=l3)
  self$means_X<-colMeans(X)
  self$stds_X<-apply(X,2,sd)
  self$mean_y<-mean(y)
  self$scale=scale
  
  
  
  
  
  # Public method for fitting the model
  fit <- function( X, theta_init, y, lambda,beta_bound, t=1, tol=1e-2, max_iter=5000, eps=1e-8,l1=21,l2=14,l3=2, to_add_strong=FALSE, add_strong=0) {
    
    n_mains<-l1+l2+l3
    eps<-1e-8*lambda
    p <- ncol(X)
    n<-nrow(X)
    if(self$scale==TRUE){
      print('was scaled')
    X <- scale(X)}#standardize X
    y <- scale(y, center = TRUE, scale = FALSE) #center y # MAYBE ALREADY SCALED ACTUALLY???????????????????????
    
    # Initialize variables
    delta <- 1-t*eps
    Theta_hat <- sets_0s_theta_2way3(theta_init,l1=l1,l2=l2,l3=l3) #matrix form
    if(to_add_strong==FALSE)
    {add_strong<-0*Theta_hat}
    Omega_hat<- (Theta_hat + t(Theta_hat) )/2 ### INIT OMEGA HAT
    U_hat<-Omega_hat*0 ###### INIT U HAT
    r_hat <- matrix(-1, nrow = n, ncol = 1)
    
    gradient_stop=0
    for (it in 2:max_iter+1) {
      vec_theta_hat <- get_theta_vec_2way3(Theta_hat,l1=l1,l2=l2,l3=l3) ###select only good positions and take as sum /6
      r_hat_old <- r_hat
      r_hat <-y - X %*%   vec_theta_hat 
      r_hat<- -  r_hat ############################### TAKE CARE ! why like this? ?? why rhat/2 or not??? ########################
      if (it%%50==4)
      {cat(' Loss rhat',mean(r_hat^2))} ###/2 to be like in paper ????????????????????????????/ depends how is r_hat
      #cat("Theta", Theta_hat) 
      
      
      

      for (j in c(1:n_mains) ) { 
         
        Xj<-get_cols_Xj3(X=X,j=j,l1=l1,l2=l2,l3=l3)
        #cat("beta bound j", beta_bound[j])
        Theta_hat[,j] <-  ONEROW_SEQ12(delta * Theta_hat[,j] - t* t(Xj)%*%r_hat -t*add_strong[,j], lambda=lambda, t=t, beta_bound_j = beta_bound[j],j=j) 
      }
      if (it>=3)
      {if (mean((r_hat_old- r_hat)^2) <=tol) #TAKE CARE RHAT IS VECTOR  !!! Sum because already scaled !!!!!!
      {  cat("Converged at iteration ",it)
        self$theta_hat=Theta_hat
        self$vec_theta_hat=get_theta_vec_2way3(Theta_hat,l1=l1,l2=l2,l3=l3)
        return(self) }
      
        #if ( sum(r_hat_old^2)- sum(r_hat^2) <0) #TAKE CARE RHAT IS VECTOR  !!! Sum because already scaled !!!!!!
        #{ gradient_stop<-gradient_stop+1
        #print('decreased gradient')
        #t<-t/1.01} 
          
        if (gradient_stop>30){
          cat("Stopped because of gradient at iteration ",it)
          self$theta_hat=Theta_hat
          self$vec_theta_hat=get_theta_vec_2way3(Theta_hat,l1=l1,l2=l2,l3=l3)
          return(self) }
      }
      
    }
    
    self$theta_hat=Theta_hat
    self$vec_theta_hat=get_theta_vec_2way3(Theta_hat,l1=l1,l2=l2,l3=l3)
    print(self$vec_theta_hat)
    cat("It has not converged. The difference between last 2 residuals is:", abs(r_hat[it-1]- r_hat[it-2]))
    return(self) 
  }
  
  
  fitstrong<-function(X, beta_bound, Theta_init, y, lambda,
                      t=1, tol=1e-5, max_iter=1000, eps=1e-8, scale=FALSE,l1=0,l2=0,l3=0, rho=1, iter_strong=100)
  { 
    print("Fit strong") 
    Omega_hat<- (Theta_init + t(Theta_init) )/2 ### INIT OMEGA HAT
    U_hat<-Omega_hat*0 ###### INIT U HAT
    

    Theta_hat_old<-Theta_init
    
    for (it in 2:iter_strong+1) {
      #print("here")
      results<- fit(X=X, beta_bound = beta_bound, theta_init=Theta_hat_old, y=y, lambda=lambda,
                    t=t, tol=tol, max_iter=max_iter, eps=eps,l1=l1,l2=l2,l3=l3, add_strong = rho*(Theta_hat_old-Omega_hat) +U_hat, to_add_strong = TRUE  )

      Theta_hat_new<-results$theta_hat #update only Theta
      
      #cat("add strong",  rho*(Theta_hat_old-Omega_hat) +U_hat )
      Omega_hat<-(Theta_hat_new+t(Theta_hat_new))/2 + (U_hat + t(U_hat))/(2*rho)
      U_hat<-U_hat+rho*(Theta_hat_new-Omega_hat)
      
      
      
      if ( mean(abs(Theta_hat_new-Theta_hat_old)) <=  tol )
      {print("STRONG CONVERGED!")
        self$Theta_hat<-Theta_hat_new
        self$vec_theta_hat<-get_theta_vec_2way3(Theta_hat =Theta_hat_new, l1=l1, l2=l2, l3=l3)

        return (self)}
      
      if(it%%5==1)
      {cat("Theta dif:",  mean(abs(Theta_hat_new-Theta_hat_old)) )}
      
      Theta_hat_old<-Theta_hat_new
      
    }
    
    print(" STRONG ADMM did not converge.")
    self$Theta_hat<-Theta_hat_new
    self$vec_theta_hat<-get_theta_vec_2way3(Theta_hat =Theta_hat_new, l1=l1, l2=l2, l3=l3)
    return (self)
  } 
  
  # method for predicting
  
  predict <- function(self, new_X) {
    
    ## Later add also rescaling
    # Implement prediction code here
    p <- ncol(new_X)
    n<-nrow(new_X)
    if (self$scale == TRUE)
    {X <- (new_X-self$means_X)/self$stds_X} ### scale X
    y_pred= X %*%  self$vec_theta_hat + self$mean_y
    return(y_pred)
  }
  
  R2_score <- function(self, new_X, y_true, verbose= TRUE) { 
    y_pred = predict(self,new_X) 
    
    if (verbose == TRUE)
    {cat ("r2 score is ", r2(y_true, y_pred))
      plot(y_pred, y_true)}
    
    return(r2(y_true, y_pred))
  }
  
  
  return(list( fit = fit, fitstrong=fitstrong, predict = predict, R2_score = R2_score, self = self))
}

##### How to use the class ######

## First create a basic dataset###




## CREATE DATASET
# Set a seed for reproducibility
set.seed(123)


row1<-c(1,1,1,1)
row2<-c(1,1,1,-1)
row3<-c(1,1,-1,1)
row4<-c(1,1,-1,-1)

row5<-c(1,-1,1,1)
row6<-c(1,-1,1,-1)
row7<-c(1,-1,-1,1)
row8<-c(1,-1,-1,-1)

row9<-c(-1,1,1,1)
row10<-c(-1,1,1,-1)
row11<-c(-1,1,-1,1)
row12<-c(-1,1,-1,-1)

row13<-c(-1,-1,1,1)
row14<-c(-1,-1,1,-1)
row15<-c(-1,-1,-1,1)
row16<-c(-1,-1,-1,-1)

#X<-rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12, row13, row14, row15, row16)



###MAKE 2^4 times each matrix

#for (i in 1:5) {
#X <- rbind(X, X)}





# Create a beta vector of length 32
#beta <- numeric(4)
#for (i in 1:length(beta)) {
#  beta[i] <- i %% 5
#}

#print(beta)
#noise <- rnorm(512, mean = 0, sd = 0.1)
#print(dim(X))
# Generate response variable
#y <- X%*%beta+noise


#library(caret)

# Assuming X and y are your predictor variables and response variable, respectively

# Set the seed for reproducibility
set.seed(123)

# Split the data into training and testing sets
#train_index <- createDataPartition(y, p = 0.8, list = FALSE)
#X_train <- X[train_index, ]
#X_test <- X[-train_index, ]
#y_train <- y[train_index]
#y_test <- y[-train_index]


#print(length(y_train))
#print(dim(X_train))

#lm_model <- lm(y_train ~ ., data = as.data.frame(X_train))

#Predict on X_test
#predictions_test <- predict(lm_model, newdata = as.data.frame(X_test))
#predictions_train<-predict(lm_model, newdata =  as.data.frame(X_train))

#print(r2(y_train, predictions_train))
#print(r2(y_test,predictions_test))







############ use the class on syntetic dataset ###


#lambda <-0
#t<-0.0005
#eps=1e-8
#theta_init<-array(rnorm(4*4*4, mean = 0, sd = 2), dim = c(4,4,4))
#beta_bound<-matrix(2,nrow=4,ncol=4)
#beta_bound[c(1,2,3),c(1,2,3)]<-1






# Example usage:

# Create an instance of the WeakHierNet class
#myWeakHierNet_seq <- WeakHierNet_seq3(X=X_train, theta_init=theta_init, y=y_train, beta_bound=beta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
#                                    l1=1,l2=1,l3=1)

# Fit the model
#fitted=myWeakHierNet_seq3$fit(X=X_train, theta_init=theta_init, y=y_train, beta_bound=beta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
#                            eps=1e-8,l1=1,l2=1,l3=1)
#print(fitted$vec_theta_hat)
#print(beta)
#print(fitted$theta_hat)

# Make predictions
#new_X <- X_test

#print(dim(X_train))
#print(dim(X_test))
#print(length(colMeans(X_train)))
#(X_test-colMeans(X_train)) /apply(X_train, 2, sd)
#print(length(fitted$vec_theta_hat))
#pred_train<-scale(X_train)%*%fitted$vec_theta_hat+mean(y_train)
#pred_test<- ( (X_test-colMeans(X_train))/apply(X_train,2,sd) ) %*%fitted$vec_theta_hat+mean(y_train)
#print(r2(y_train,pred_train))
#print(r2(y_test,pred_test))

#myWeakHierNet_seq$R2_score(self=fitted, new_X= as.matrix(new_X), y_true = y_test, verbose = TRUE)

#fitted


### Use on another synthetic dataset is columns continous values from normal distribution

# Set the number of predictors
#num_predictors <- 68

# Set the beta values
#beta <- c(5, 3, 7, 9, 10, rep(1e-4, num_predictors - 5))

# Generate the predictor variables
#set.seed(123)  # for reproducibility
#predictors <- matrix(rnorm(num_predictors * 1000, mean = 0, sd = 5), ncol = num_predictors)

# Generate the noise
#noise <- rnorm(1000, mean = 0, sd = 0.3)

# Generate the response variable
#response <- predictors %*% beta + noise

# Combine predictors and response into a data frame
#data <- data.frame(response = response, predictors)


# Split predictors and response into training and testing sets
#set.seed(123)  # for reproducibility
#train_index <- sample(1:1000, 700)  # 70% for training
#test_index <- setdiff(1:1000, train_index)

#X_train <- predictors[train_index, ]
#X_test <- predictors[test_index, ]
#y_train <- response[train_index]
#y_test <- response[test_index]

#lm_model <- lm(y_train ~ ., data = as.data.frame(X_train))

# Make predictions on the testing data
#predictions_test <- predict(lm_model, newdata = as.data.frame(X_test))

# Calculate the R-squared value on the testing data
#r2(y_test, predictions_test)



##
#lambda <-3e2
#t<-8e-5
#eps=1e-8
#theta_init<-array(rnorm(11*11*11, mean = 0, sd = 2), dim = c(11,11,11))
#beta_bound<-matrix(1,nrow=11,ncol=11)






# Example usage:

# Create an instance of the WeakHierNet class
#myWeakHierNet_seq <- WeakHierNet_seq3(X=X_train, theta_init=theta_init, y=y_train, beta_bound=beta_bound, lambda=lambda, t=t, tol=1e-7, max_iter=5000, eps=1e-8,
#                                    l1=5,l2=2,l3=2)

# Fit the model
#fitted=myWeakHierNet_seq3$fit(X=X_train, theta_init=theta_init, y=y_train, beta_bound=beta_bound, lambda=lambda, t=t, tol=1e-7, max_iter=5000, 
#                             eps=1e-8,l1=5,l2=2,l3=2)
#print(fitted$vec_theta_hat/fitted$stds_X)
#print(beta)
#print(fitted$theta_hat)

# Make predictions
#new_X <- X_test
#print(dim(X_test))
#pred_train<-scale(X_train)%*%fitted$vec_theta_hat+mean(y_train)
#pred_test<-scale(X_test, center = FALSE)%*%fitted$vec_theta_hat+mean(y_train)
#pred_test2<- ((X_test-colMeans(X_train))/apply(X_train,2,sd) ) %*%fitted$vec_theta_hat+mean(y_train)
#print(r2(y_train,pred_train))
#print(r2(y_test,pred_test2))



#print('USING METHODS FROM MY CLASS')

#print(myWeakHierNet_seq$R2_score(fitted,X_test,y_test))




################ GANDURI ###################

## FA COLS CU MEAN 0 prin tehnica cu -1 in 
## Vf COD iar
### Fa pred si score methods


###TAKE CARE WHICH KNOTS I USE FOR e)
####
#### HOW TO DO SCALE AND MEAN IN PRED AND R2 SCORE
#### VEZI CUM E FECETUL LUI THETA BOUND DUPA SCALE SAU CAND E MIC
##why small numerical instability or variance from one run to another

