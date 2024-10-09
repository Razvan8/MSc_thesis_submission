

set_positioned_value_in_theta<-function(theta, val, i, j, threshold)
{  
  theta[i,j]<- set_0_threshold(val/2, threshold=threshold)
 theta[j,i]<-set_0_threshold(val/2, threshold = threshold)
 return(theta)}

set_positioned_value_in_psi<-function(psi, val, i, j, k, threshold)
{  
  psi[i,j,k]<- set_0_threshold(val, threshold=threshold) ##use it with the right val already in function
  psi[i,k,j]<-set_0_threshold(val, threshold = threshold)
  psi[j,i,k]<- set_0_threshold(val, threshold=threshold)
  psi[j,k,i]<-set_0_threshold(val, threshold = threshold)
  psi[k,i,j]<- set_0_threshold(val, threshold=threshold)
  psi[k,j,i]<-set_0_threshold(val, threshold = threshold)
  
  return(psi)}


get_val_psi<-function(psi, i,j,k){
  val<-psi[i,j,k] + psi[i,k,j] + psi[j,i,k] + psi[j,k,i]+ psi[k,i,j]+ psi[k,j,i]
  return(val)}

assert <- function(condition, message) {
  if (!condition)
    stop(message)
}

set_0_threshold<- function(x,threshold)
{ if (abs(x)<=threshold)
{return (0)}
  return (x)
}


#Function transform maineff+idx to position theta
get_position_in_table<-function(main, idx, l1, l2, l3)
{assert(main ==1 | main==2 | main ==3, " main should be 1 2 or 3, i.e. which main effect it is")
  assert(idx<=max(l1,l2,l3), "idx is too big")
  pre_position<-c(0,l1, l1+l2)
  position_theta<-pre_position[main]+ idx
  return(position_theta)
}

get_position_in_table(3,1,4,4,4)



get_all_beta<- function(beta, l1, l2, l3, threshold=1e-3)
{ assert(length(beta)==l1+l2+l3, "length of coefs is not right")
  beta_old<-beta
 beta_new<-array(0, dim=l1+l2+l3+3) # shape of new beta

 beta_new[1:l1]<-beta_old[1:l1]
 beta_new[l1+1]<- set_0_threshold(x=-sum(beta_old[1:l1]), threshold = threshold)
 
 beta_new[(l1+2):(l1+l2+1)]<-beta_old[(l1+1): (l1+l2)]
 beta_new[l1+l2+2]<- set_0_threshold(x=-sum(beta_old[(l1+1): (l1+l2)]), threshold = threshold) 
 
 beta_new[(l1+l2+3):(l1+l2+l3+2)]<-beta_old[(l1+l2+1): (l1+l2+l3)]
 beta_new[l1+l2+l3+3]<- set_0_threshold(x = -sum(beta_old[(l1+l2+1): (l1+l2+l3)]), threshold = threshold)
 
 return(beta_new)
}

#get_all_beta(beta=c(1,2,3,4,5,6,7),l1=3, l2=2, l3=2, threshold = 1e1)

get_all_theta <-function(theta, l1, l2, l3, threshold=1e-3)
{ theta_old<- (theta+t(theta))/2 #effect symmetric

range1_old<-c(1:l1)
range2_old<-c( (l1+1): (l1+l2)  )
range3_old<-c( (l1+l2+1) : (l1+l2+l3)   )
range1_new<-c(1:(l1+1))
range2_new<-c( (l1+2): (l1+l2+2)  )
range3_new<-c( (l1+l2+3) : (l1+l2+l3+3)   )


ls_range_old<-list(range1_old, range2_old, range3_old)
ls_range_new<-list(range1_new, range2_new, range3_new)

#Step1 create theta new
theta_new<- matrix(0, nrow = l1+l2+l3+3, ncol=l1+l2+l3+3) 

for (i in c(1:3))
{for (j in c(1:3))
{rangenew_i<-unlist(ls_range_new[i])
 rangenew_j<-unlist(ls_range_new[j])
 rangeold_i<-unlist(ls_range_old[i])
 rangeold_j<-unlist(ls_range_old[j])
    theta_new[rangenew_i[-length(rangenew_i)], rangenew_j[-length(rangenew_j)]] <- theta_old[rangeold_i, rangeold_j]}}


##STEP 2 Create matrix of coefs
theta<- (theta_new+t(theta_new) )/2

#STEP 3 Identify (ab)_{il1} ...
for (i in c(1:l1)){

  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum( t( theta_old[i, range2_old ]) +  theta_old[range2_old, i ] )
                                         , i= get_position_in_table(main=1, idx=i, l1=l1+1, l2=l2+1, l3=l3+1),
                                           j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) #(ab)_i l2+1
  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum( t(theta_old[i, range3_old ]) +  theta_old[range3_old, i ] )
                                       , i= get_position_in_table(main=1, idx=i, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold)} #(ac)_i l3+1


  
for (j in c(1:l2)){
  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum( t( theta_old[j+l1, range1_old ]) +  theta_old[range1_old, j+l1 ] )
                                       , i= get_position_in_table(main=2, idx=j, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) #(ab)_l1+1 j
  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum ( t(theta_old[j+l1, range3_old ]) +  theta_old[range3_old, j+l1 ] )
                                       , i= get_position_in_table(main=2, idx=j, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) }  #(bc)_j l3+1


for (k in c(1:l3)){
  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum( t( theta_old[k+l1+l2, range1_old ]) +  theta_old[range1_old, k+l1+l2 ] )
                                       , i= get_position_in_table(main=3, idx=k, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) #(ac)_l1+1 k
  theta<-set_positioned_value_in_theta(theta = theta, val = 
                                         -sum( t(theta_old[k+l1+l2, range2_old ]) +  theta_old[range2_old, k+l1+l2 ] )
                                       , i= get_position_in_table(main=3, idx=k, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) }  #(bc)_l2+1 k

  

  
  

  #STEP 4 IDENTIFY (ab)_{l1+ l2+1}...
theta<-set_positioned_value_in_theta(theta = theta, val = 
                                       -sum( t(theta[l1+l2+2, range1_new[-length(range1_new)]]) +  theta[range1_new[-length(range1_new)], l1+l2+2 ] ) #fixed l2
                                     , i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     threshold = threshold)  #l1+1 l2+1

theta<-set_positioned_value_in_theta(theta = theta, val = 
                                       -sum( t(theta[l1+l2+l3+3, range1_new[-length(range1_new)] ]) +  theta[range1_new[-length(range1_new)], l1+l2+l3+3 ] ) #fixed l3+`1`
                                     , i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     j= get_position_in_table(main=3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     threshold = threshold)   #l1+1 l3+1

theta<-set_positioned_value_in_theta(theta = theta, val = 
                                       -sum( t(theta[l1+l2+l3+3, range2_new[-length(range2_new)] ]) +  theta[range2_new[-length(range2_new)], l1+l2+l3+3 ] ) #fixed l3+1
                                     , i= get_position_in_table(main=3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     threshold = threshold)   #l2+1 l3+1
 





assert(all(theta==t(theta)), "theta is not symmetric but it should be ")
  
return(theta)  
  
}

#theta<-matrix(1, nrow=7,ncol=7)
#theta[c(1:2), c(1:2)]<-0
#theta[c(3:4), c(3:4)]<-0
#theta[c(5:7), c(5:7)]<-0
#theta[2,4]<-5
#print(theta)
#get_all_theta(theta, l1=2,l2=2,l3=3 )



get_all_psi<-function(psi, l1, l2, l3, threshold=1e-3)
{
 

  psi_old<- psi*0 #initialize
  
  range1_old<-c(1:l1)
  range2_old<-c( (l1+1): (l1+l2)  )
  range3_old<-c( (l1+l2+1) : (l1+l2+l3)   )
  range1_new<-c(1:(l1+1))
  range2_new<-c( (l1+2): (l1+l2+2)  )
  range3_new<-c( (l1+l2+3) : (l1+l2+l3+3)   )
  
  
  ls_range_old<-list(range1_old, range2_old, range3_old)
  ls_range_new<-list(range1_new, range2_new, range3_new)
  
  for (i in 1:dim(psi)[1]) {
    for (j in 1:dim(psi)[2]) {
      for (k in 1:dim(psi)[3]) {
        #Symmetric psi_old
        psi_old[i, j, k] <- (psi[i, j, k] + psi[i, k, j] + psi[j, i, k] + ##psi old already symmtric
                           psi[j, k, i] + psi[k, i, j] + psi[k, j, i]) / 6
      }
    }
  }
  
  #STEP 1 CREATE PSI NEW 
  psi_new<-array(0,dim=c(l1+l2+l3+3, l1+l2+l3+3,l1+l2+l3+3)  )
  
  for (i in c(1:3))
  {for (j in c(1:3))
  {for (k in c(1:3))
  {
    
  rangenew_i<-unlist(ls_range_new[i])
  rangenew_j<-unlist(ls_range_new[j])
  rangenew_k<-unlist(ls_range_new[k])
  rangeold_i<-unlist(ls_range_old[i])
  rangeold_j<-unlist(ls_range_old[j])
  rangeold_k<-unlist(ls_range_old[k])
  
  psi_new[rangenew_i[-length(rangenew_i)], rangenew_j[-length(rangenew_j)], rangenew_k[-length(rangenew_k)]] <- psi_old[rangeold_i, rangeold_j, rangeold_k] }}}
  
  #STEP 2 PSI SYMMTERIC 
  for (i in 1:dim(psi_new)[1]) {
    for (j in 1:dim(psi_new)[2]) {
      for (k in 1:dim(psi_new)[3]) {
        #Symmetric psi_old
        psi_new[i, j, k] <- (psi_new[i, j, k] + psi_new[i, k, j] + psi_new[j, i, k] +
                               psi_new[j, k, i] + psi_new[k, i, j] + psi_new[k, j, i]) / 6
      }
    }
  }
  psi<-psi_new
  
  
  ###STEP 3: one l+1
  for (i in c(1:l1)){
    for(j in c(1:l2)){
    
    psi<-set_positioned_value_in_psi(psi = psi, val = 
                                           -sum(  psi_old[i,l1+j, range3_old ] )
                                         , i= get_position_in_table(main=1, idx=i, l1=l1+1, l2=l2+1, l3=l3+1),
                                         j= get_position_in_table(main=2, idx=j, l1=l1+1, l2=l2+1, l3=l3+1),
                                         k=get_position_in_table(main = 3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                         threshold = threshold) }} # (abc)_{i j l3+1}
  

  for (i in c(1:l1)){
    for(k in c(1:l3)){
      
      psi<-set_positioned_value_in_psi(psi = psi, val = 
                                         -sum(  psi_old[i, range2_old, l1+l2+k] )
                                       , i= get_position_in_table(main=1, idx=i, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       k=get_position_in_table(main = 3, idx=k, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) }} # (abc)_{i l2+1 k}
  
  for (j in c(1:l2)){
    for(k in c(1:l3)){
      
      psi<-set_positioned_value_in_psi(psi = psi, val = 
                                         -sum(  psi_old[range1_old, l1+j, l1+l2+k ] )
                                       , i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=2, idx=j, l1=l1+1, l2=l2+1, l3=l3+1),
                                       k=get_position_in_table(main = 3, idx=k, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) }} # (abc)_{l1+1 j k}
  
  #print(psi)
  
  
  
  #STEP 4: 2 of l
  for (i in c(1:l1)){
      
      psi<-set_positioned_value_in_psi(psi = psi, val =   ## i and l2+1 fixed
                                         -sum( psi[ i, l1+l2+2, range3_new[-length(range3_new)]  ]),
                                        i= get_position_in_table(main=1, idx=i, l1=l1+1, l2=l2+1, l3=l3+1),
                                       j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       k=get_position_in_table(main = 3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                       threshold = threshold) } # (abc)_{i l2+1 l3+1}
  
  for (j in c(1:l2)){
    
    psi<-set_positioned_value_in_psi(psi = psi, val =   ## l1+1 and j fixed
                                       -sum( psi[ l1+1, j+l1+1, range3_new[-length(range3_new)]  ]),
                                     i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     j= get_position_in_table(main=2, idx=j, l1=l1+1, l2=l2+1, l3=l3+1),
                                     k=get_position_in_table(main = 3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     threshold = threshold) } # (abc)_{l1+1 j l3+1}
  
  for (k in c(1:l3)){
    
    psi<-set_positioned_value_in_psi(psi = psi, val =   ## l1+1 and k fixed
                                       -sum( psi[ l1+1, range2_new[-length(range2_new)], k+l1+l2+2  ]),
                                     i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                     k=get_position_in_table(main = 3, idx=k, l1=l1+1, l2=l2+1, l3=l3+1),
                                     threshold = threshold) } # (abc)_{l1+1 l2+1 k}
  
  
  ##STEP 5: 3 of l
  psi<-set_positioned_value_in_psi(psi = psi, val =   ## l1+1 and k fixed
                                      -sum( psi[ range1_new[-length(range1_new)], l1+l2+2, l1+l2+l3+3  ]),
                                    i= get_position_in_table(main=1, idx=l1+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                    j= get_position_in_table(main=2, idx=l2+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                    k=get_position_in_table(main = 3, idx=l3+1, l1=l1+1, l2=l2+1, l3=l3+1),
                                    threshold = threshold)  # (abc)_{l1+1 l2+1 l3+1}
  
  
  
  #print(psi_old)
  return(psi)
  
  
}



#psi<-array(1,dim=c(7,7,7))
#psi[c(1:3), c(1:3), ]<-0
#psi[c(1:3), , c(1:3) ]<-0
#psi[,c(1:3), c(1:3) ]<-0

#psi[c(4:5), c(4:5), ]<-0
#psi[c(4:5), , c(4:5) ]<-0
#psi[ , c(4:5), c(4:5) ]<-0


#psi[c(6:7), c(6:7), ]<-0
#psi[c(6:7), , c(6:7) ]<-0
#psi[ ,c(6:7), c(6:7) ]<-0

#psi[1,4,6]<-7
#psi_all<-get_all_psi(psi, l1=3, l2=2, l3=2)







