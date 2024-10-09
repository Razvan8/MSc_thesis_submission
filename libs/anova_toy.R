######################################################################
## ANOVA FACTORIZATION FOR TOY DATA
######################################################################

# Clear all
rm(list=ls())
gc()
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
cat("\014")

# Memory
memory.limit(64000)



#######################
## Dummy Function
#######################

dummy.matrix <- function(NF=2, NL=rep(2,NF)) {
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
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    
    # Dummy design matrix 2-way
    n.2w <- L.12
    p.2w <- n.2w-1
    x.2w <- data.frame(matrix(0, nrow=n.2w, ncol=p.2w))
    rownames(x.2w) <- c(fac.12)
    colnames(x.2w) <- c(fac.1[-L.1], fac.2[-L.2],
                        sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":")))
    )
    for (col in 1:ncol(x.2w)) {
      col.tags <- unlist(strsplit(colnames(x.2w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.2w[grepl(col.tags,split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[1]]),col] <- 1
        x.2w[grepl(col.tags,split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[2]]),col] <- 1
        x.2w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[1]]),col] <- -1
        x.2w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[2]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.2w)==col.tags[1])
        col.2 <- which(colnames(x.2w)==col.tags[2])
        x.2w[,col] <- x.2w[,col.1]*x.2w[,col.2] 
      }
    }
    
    return(x.2w)
  }
  
  if (NF == 3) {
    
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    L.13 <- L.1*L.3
    L.23 <- L.2*L.3
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep=":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep=":")))
    
    # Three-ways
    L.123 <- L.1*L.2*L.3
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep=":")))
    
    # Dummy design matrix 3-way
    n.3w <- L.123
    p.3w <- n.3w-1
    x.3w <- data.frame(matrix(0, nrow=n.3w, ncol=p.3w))
    rownames(x.3w) <- c(fac.123)
    colnames(x.3w) <- c(fac.1[-L.1], fac.2[-L.2], fac.3[-L.3],
                        sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))), paste, sep=":")))
    )
    for (col in 1:ncol(x.3w)) {
      col.tags <- unlist(strsplit(colnames(x.3w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[1]]),col] <- 1
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[2]]),col] <- 1
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[3]]),col] <- 1
        x.3w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[1]]),col] <- -1
        x.3w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[2]]),col] <- -1
        x.3w[grepl(paste(fac.tag[1],L.3,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[3]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.3w)==col.tags[1])
        col.2 <- which(colnames(x.3w)==col.tags[2])
        x.3w[,col] <- x.3w[,col.1]*x.3w[,col.2]
      }
      if (length(col.tags)==3) {
        col.1 <- which(colnames(x.3w)==col.tags[1])
        col.2 <- which(colnames(x.3w)==col.tags[2])
        col.3 <- which(colnames(x.3w)==col.tags[3])
        x.3w[,col] <- x.3w[,col.1]*x.3w[,col.2]*x.3w[,col.3]
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
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep=".")
    fac.4 <- paste(fac.tags[4], fac.levs[1:L.4], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    L.13 <- L.1*L.3
    L.14 <- L.1*L.4
    L.23 <- L.2*L.3
    L.24 <- L.2*L.4
    L.34 <- L.3*L.4
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep=":")))
    fac.14 <- sort(as.vector(outer(fac.1, fac.4, paste, sep=":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep=":")))
    fac.24 <- sort(as.vector(outer(fac.2, fac.4, paste, sep=":")))
    fac.34 <- sort(as.vector(outer(fac.3, fac.4, paste, sep=":")))
    
    # Three-ways
    L.123 <- L.1*L.2*L.3
    L.124 <- L.1*L.2*L.4
    L.134 <- L.1*L.3*L.4
    L.234 <- L.2*L.3*L.4
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep=":")))
    fac.124 <- sort(as.vector(outer(fac.1, fac.24, paste, sep=":")))
    fac.134 <- sort(as.vector(outer(fac.1, fac.34, paste, sep=":")))
    fac.234 <- sort(as.vector(outer(fac.2, fac.34, paste, sep=":")))
    
    # Four-ways
    L.1234 <- L.1*L.2*L.3*L.4
    fac.1234 <- sort(as.vector(outer(fac.1, fac.234, paste, sep=":")))
    
    # Dummy design matrix 4-way
    n.4w <- L.1234
    p.4w <- n.4w-1
    x.4w <- data.frame(matrix(0, nrow=n.4w, ncol=p.4w))
    rownames(x.4w) <- c(fac.1234)
    colnames(x.4w) <- c(fac.1[-L.1], fac.2[-L.2], fac.3[-L.3], fac.4[-L.4],
                         sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], fac.3[-L.3], paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], fac.4[-L.4], paste, sep=":"))),
                         sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))),
                         sort(as.vector(outer(fac.2[-L.2], fac.4[-L.4], paste, sep=":"))),
                         sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))), paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                         sort(as.vector(outer(fac.2[-L.2], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                         sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))),
                                                                                paste, sep=":"))), paste, sep=":")))
    )
    for (col in 1:ncol(x.4w)) {
      col.tags <- unlist(strsplit(colnames(x.4w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[1]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[2]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[3]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[4]]),col] <- 1
        x.4w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[1]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[2]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.3,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[3]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.4,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[4]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]
      }
      if (length(col.tags)==3) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        col.3 <- which(colnames(x.4w)==col.tags[3])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]*x.4w[,col.3]
      }
      if (length(col.tags)==4) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        col.3 <- which(colnames(x.4w)==col.tags[3])
        col.4 <- which(colnames(x.4w)==col.tags[4])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]*x.4w[,col.3]*x.4w[,col.4]
      }
    }
    
    return(x.4w)
  }

}



#######################
## Toy Model
#######################

# Dummy Matrix
x.2w <- dummy.matrix(NF=2, NL=c(3,3))
x.3w <- dummy.matrix(NF=3, NL=c(9,6,4))
x.4w <- dummy.matrix(NF=4, NL=c(20,10,3,4))



Z=t(apply(x.2w[,c(1:4)], 1, function(x) outer(x, x, "*")))
minus_ones_count <- colSums(Z == -1)

# Identify columns with no occurrences of -1
no_minus_ones_columns <- which(minus_ones_count == 0)

# Set columns with no -1 to 0
Z[, no_minus_ones_columns] <- 0
Z

#x<-cbind(x.3w, rep(1, nrow(x.3w)))

#svd_result <- svd(x)

# Extract the singular values
#singular_values <- svd_result$d

# Find the largest and smallest singular values
#largest_singular_value <- max(singular_values)
#smallest_singular_value <- min(singular_values)
#largest_singular_value/smallest_singular_value

#libs_path<-file.path("..","libs")
#source(file.path(libs_path,'helper_functions.R'))
#data<-load_deoxi_flourination()

#colnames(data)

#X<-data[,c('a','b','s')]
#one_hot_encoding <- model.matrix(~ . - 1, X)

#dim(one_hot_encoding)
#data$y[data$y>100]

#unique(data$b)

# Hierarchical Coefficients (2way)
p.2w <- ncol(x.2w)
n.2w <- p.2w + 1
beta.min <- 1
beta.max <- 10
beta.true <- data.frame(rep(0, n.2w))
rownames(beta.true) <- c("interc", colnames(x.2w))
colnames(beta.true) <- c("coeffs")
beta.true$coeffs <- runif(n.2w, beta.min, beta.max)*sample(c(1,-1),size=n.2w,replace=TRUE)
levs.true <- c("interc","A.1", "A.2", "B.1", "A.1:B.1")
beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
beta.true

# Response vector (2way)
sigma.y <- 2.5
y.2w <- data.frame(row.names=rownames(x.2w))
y.2w$obs <- beta.true$coeffs[1] + as.matrix(x.2w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.2w), 0, sigma.y)
y.2w$true <- beta.true$coeffs[1] + as.matrix(x.2w)%*%as.vector(beta.true$coeffs)[-1]
y.2w



# Hierarchical Coefficients (2way)
p.3w <- ncol(x.3w)
n.3w <- p.3w + 1
beta.min <- 1
beta.max <- 10
beta.true <- data.frame(rep(0, n.3w))
rownames(beta.true) <- c("interc", colnames(x.3w))
colnames(beta.true) <- c("coeffs")
beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
levs.true <- c("interc","A.1", "A.2", "A.3","A.4", "A.5", "A.6", "A.7",  "B.1", "B.2","B.3", "B.4", "C.1", "C.2","C.3",  "A.1:B.1","A.1:B.2","A.1:B.3",
               "A.2:B.1","A.2:B.2","A.2:B.3", "A.3:B.1","A.3:B.2","A.3:B.3", "A.5:B.1","A.5:B.2","A.5:B.3", "B.1:C1","B.2:C2","B.3:C3", "B4:C1", "B.4:C3", "B4:C2", 
               "A.7:B.1","A.7:B.2","A.7:B.3","A.6:B.1","A.6:B.2","A.6:B.3",
               "A.2:B.1:C.3","A.2:B.2:C.3","A.2:B.3:C.3", "A.3:B.1:C.2","A.3:B.2:C.2","A.3:B.3:C.1","A.3:B.3:C.2","A.3:B.3:C.3", "A.5:B.1:C.3","A.5:B.2:C.3" )
beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
beta.true

# Response vector (2way)
sigma.y <- 4
y.3w <- data.frame(row.names=rownames(x.3w))
y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
y.3w




















