

###Create a load_data function

load_deoxi_flourination<-function(data_path=file.path("..","Data"))
  
{library(dplyr)
  # Constructing the path to the "libs" folder
  
  
  # Combine path with filename
  train_path <- file.path(data_path, "fluorination_train.csv")
  test_path<- file.path(data_path, "fluorination_test.csv")
  
  
  
  
  agd_palette <- c("#0063C3", "#00BCD0", "#CB7C85",  "#E8C9BB", "#252A5A", 
                   "#7D2C42", "#A9A9A9")
  
  # original data
  dat_1 <- read.csv(train_path, stringsAsFactors = T)
  dat_2 <- read.csv(test_path, stringsAsFactors = T)
  
  dat <- rbind(dat_1, dat_2)
  
  colnames(dat)
  
  dat <- dat %>%
    rename(s = sulfonyl.fluoride, a = alcohol, b = base) %>%
    mutate(s = sub("\\.", "-", s)) %>%
    mutate(s = gsub("3-", "4-", s)) %>%
    mutate(s = factor(s))
  
  dat[,c(1,2,3,length(colnames(dat)))]
  dat$yield[which(dat$yield>100)]
  
  # reactions data
  rxns <- dat %>% select(a, b, s, yield)
  
  # pre-process factor like features from the original publication
  da_orig <- dat %>% select("a", contains("alcohol")) %>% distinct()
  
  da_fact <- da_orig %>%
    select("a", where(is.integer)) %>%
    rename_with(~ sub("alcohol_", "", .x)) %>%
    rename_with(~ sub("\\.", "_", .x))
  
  # coalesce one-hot spread features that's not needed and only takes space
  order <- da_fact %>%
    select(primary, secondary, tertiary) %>%
    mutate_all(as.character) %>%
    mutate(secondary = recode(secondary, "1" = "2"),
           tertiary = recode(tertiary, "1" = "3")) %>%
    mutate(order = pmax(primary, secondary, tertiary)) %>%
    mutate(order = recode(order, "1" = "primary",
                          "2" = "secondary",
                          "3" = "tertiary")) %>%
    pull(order)
  
  ring_size <- da_fact %>%
    select(`4_membered_ring`, `5_membered_ring`,
           `6_membered_ring`, `7_membered_ring`) %>%
    mutate_all(as.character) %>%
    mutate(r4 = recode(`4_membered_ring`, "1" = "4", "0" = "NA"),
           r5 = recode(`5_membered_ring`, "1" = "5", "0" = "NA"),
           r6 = recode(`6_membered_ring`, "1" = "6", "0" = "NA"),
           r7 = recode(`7_membered_ring`, "1" = "7", "0" = "NA")) %>%
    mutate(ring_size = pmin(r4, r5, r6, r7)) %>%
    mutate(ring_size = recode(ring_size, "NA" = "0")) %>%
    pull(ring_size)
  
  da_fact <- da_fact %>%
    mutate(ring_size = ring_size, order = order) %>%
    select(-contains("membered_ring")) %>%
    select(-primary, -secondary, -tertiary, -cyclic) %>%
    mutate_all(as.factor)
  
  
  # read alcohol DFT features
  glb_feats <- c("number_of_atoms", "molar_mass", "electronegativity",
                 "electronic_spatial_extent", "hardness", "homo_energy",
                 "lumo_energy", "molar_volume", "dipole", 
                 "OC_length", "OC_L", "OC_B1", "OC_B5")
  c_feats <- c("C_APT_charge", "C_Mulliken_charge", "C_NMR_shift",
               "C_NPA_charge", "C_VBur", "C_angle", "C_PVBur")
  o_feats <- c("O_APT_charge", "O_Mulliken_charge", "O_NMR_shift",
               "O_NPA_charge", "O_VBur")
  
  # M062X calculation
  da_dft_m062x <- read.csv(file.path(data_path,"alcohols_Boltzmann_M062X_THF.csv")) %>%
    select(-inchi) %>%
    rename(a = name) %>%
    select("a", all_of(glb_feats), all_of(c_feats), all_of(o_feats))
  
  da_dft <- da_dft_m062x
  da <- merge(da_dft, da_fact)
  
  
  # helper functions
  normalize_yield <- function(yield) {
    y <- ifelse(yield > 97, 100 - runif(n(), 0.1, 3),
                ifelse(yield < 3, runif(n(), 0.1, 3), yield)) / 100
    return(y)
  }
  
  
  
  # build final dataframe for analysis
  data <- rxns %>%
    mutate(prob = normalize_yield(yield)) %>%
    left_join(da, by = "a") %>%
    mutate(a = as.factor(a))
  
  return(data)
  
  print("Data loaded successfuly")
  }
  



plot_yield_histo<-function(data)

{library(ggplot2)
  
  # Combine data
  y.all <- data$prob
  n.all <- length(y.all)
  tt <- 1:n.all/n.all
  
  # Find MLE for p
  ld<-function(p)
  {
    (n.all*p/(2*p-1)-sum(y.all)+n.all/log((1-p)/p))/(p*(p-1))
  }
  p.hat <- uniroot(ld,c(0.001,0.5))$root
  
  # Plot the histogram
  
  ll <- log((1-p.hat)/p.hat)*p.hat^tt*(1-p.hat)^(1-tt)/(1-2*p.hat)
  df <- data.frame(y.all,tt,ll)
  colnames(df) <- c("Yield","Probability","Density")
  
  
  
  
  
  
  element_textbox <- function(...) {
    el <- element_text(...)
    class(el) <- c("element_textbox", class(el))
    el
  }
  
  element_grob.element_textbox <- function(element, ...) {
    text_grob <- NextMethod()
    rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
    
    ggplot2:::absoluteGrob(
      grid::gList(
        element_grob(calc_element("strip.background", theme_bw())),
        text_grob
      ),
      height = grid::grobHeight(text_grob), 
      width = grid::unit(1, "npc")
    )
  }
  
  
  ggplot(df, aes(x=Yield*100)) + 
    geom_histogram(aes(y=..density..),color="black",fill="cornflowerblue",breaks=seq(0,100,length=11))+ 
    scale_y_continuous(breaks = c(0, 0.01, 0.02)) +
    theme(text = element_text(size = 25)) +labs(title = "Yield density") +
    theme(plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    ), aspect.ratio = 1)+
    geom_line(aes(x=tt*100,y =ll/100), color = "black",size=2)+ylab("Density")+xlab("Yield")
}



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



hamming_distance_sign<-function(beta, beta_hat,scale=TRUE)
{ beta<-c(beta)
  beta_hat<-c(beta_hat)
  total=length(beta)
 correct=sum(sign(beta)==sign(beta_hat))
 if (scale==TRUE)
   return(100-correct/total*100)
 return(total-corect)}

#print(hamming_distance_sign(c(1,2,3,4),c(0,0,0,1)))



l1_loss_beta<-function(beta, beta_hat, scale=TRUE){
  beta<-c(beta)
  beta_hat<-c(beta_hat)
  result<-sum(abs(beta-beta_hat))
  if (scale== TRUE)
  {result<-result/length(beta)}
  return(result)}


#l1_loss_beta(c(-1,2,3,5,0),c(-1,-0,0,-1,5))

MSE_beta<-function(beta, beta_hat){
  beta<-c(beta)
  beta_hat<-c(beta_hat)
  return(norm(beta-beta_hat, type="2")^2/length(beta))}

#print(MSE_beta(c(1,2,3,4),c(0,0,0,1)))

TPR_zeros<-function(beta,beta_hat) #TPR-sensitivity (TP /(TP+FN)) ##predicted as 0 from 0
  { beta<-c(beta)
  beta_hat<-c(beta_hat)
  idx_beta<-which(beta==0) #pred as positives
  idx_beta_hat<-which(beta_hat==0) ## True positives
  return( length(intersect(idx_beta, idx_beta_hat))/length(idx_beta)*100)}
  
#TPR_zeros(beta = c(0, 0, 0, 0, 1), beta_hat = c(0, 0, 0, 0, 0))

FPR_zeros<-function(beta,beta_hat) #FPR pred TRUE  ##predicted as 0 from non 0
{ beta<-c(beta)
beta_hat<-c(beta_hat)
  idx_beta<-which(beta!=0) #Negatives= FP+TN 
idx_beta_hat<-which(beta_hat==0) #pred positive
return(length(intersect(idx_beta, idx_beta_hat))/length(idx_beta)*100)}

#FPR_zeros(beta = c(1, 1, 0, 0, 1), beta_hat = c(1, 0, 0, 0, 0))

all_beta_functions<-function(beta, beta_hat, scale=TRUE)
{cat("TPR zeros: ", TPR_zeros(beta, beta_hat), "\n")
cat("FPR zeros: ", FPR_zeros(beta, beta_hat), "\n")
cat("MSE beta: ", MSE_beta(beta, beta_hat), "\n")
cat("L1 loss beta: ", l1_loss_beta(beta=beta, beta_hat = beta_hat, scale = scale ))
}
  

#all_beta_functions(beta = c(-1, 1, 0, 0, 1), beta_hat = c(1, 0, 0, 0, 0))
