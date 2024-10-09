######################################################################
## COMPARISON PLSGLM AND LASSOGLM-HEREDITY ON FULL DATA
######################################################################


# Memory
memory.limit(64000)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
source(file.path(libs_path,'Shim3way_class.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))


library(glmnet)

data<-load_deoxi_flourination()

print(max(data$yield))
hist(data$prob)


data


plot_yield_histo(data)


levels(data$a)
levels(data$b)
levels(data$s)

y_true<-data$prob*100
X_data<-data[c('a','b','s')]
colnames(X_data)<-c("alcohol", "base", "sulfonyl_fluoride")

X_data <- lapply(X_data, as.factor)
X_data <- as.data.frame(X_data)  # Convert back to data frame
xs<-X_data

X_data[338,]

unique(X_data$a)
unique(X_data$b)
unique(X_data$s)

X_data[338,]
y_true[338]



# Create the dummy matrix without intercept
X_dummy <- model.matrix(~ . - 1, data = X_data, contrasts.arg = lapply(X_data, contrasts, contrasts = FALSE))

# Find rows where the 37th column has value 1
rows_to_replace1 <- X_dummy[,37] == 1
# Replace values in columns 1 to 36 with -1 for these rows
X_dummy[rows_to_replace1, 1:36] <- -1


rows_to_replace2 <- X_dummy[,41] == 1
X_dummy[rows_to_replace2, 38:40] <- -1


rows_to_replace3 <- X_dummy[,46] == 1
X_dummy[rows_to_replace3, 42:45] <- -1


print(dim(X_dummy))
print(colSums(X_dummy))

##only once
X_dummy<-X_dummy[,-c(37,41,46)]
print(dim(X_dummy))
print(colSums(X_dummy))
###

new_colnames <- c(paste0("alcohol ", 1:36), paste0("base ", 1:3), paste0("sulfonyl_fluoride ", 1:4))
colnames(X_dummy) <- new_colnames

X_dummy<-as.data.frame(X_dummy, col.names =TRUE)


for (i in 1:36) {
  for (jk in 37:43) {
    
    X_dummy[, paste0(colnames(X_dummy)[i], ':', colnames(X_dummy)[jk])] <-X_dummy[,i]*X_dummy[,jk]
  }
}

for (j in c(37:39))
{for (k in c(40:43))
{
  X_dummy[,paste0((colnames(X_dummy)[j]),':', colnames(X_dummy[k]) ) ] <- X_dummy[,j]*X_dummy[,k]
}
}




for (i in c(1:36))
{for (j in c(37:39))
{for (k in c(40:43))
{
  X_dummy[,paste0( (colnames(X_dummy)[i]), ':' , colnames(X_dummy[j]), ":",  colnames(X_dummy[k]) ) ] <- X_dummy[,i]*X_dummy[,j]*X_dummy[,k]
}
}}

X<-as.matrix(X_dummy)
y<-y_true
colnames(X)
l1=36
l2=3
l3=4
range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )



X[338,range_main]



#######################
## Data processing
#######################

x.all<-X[,range_main]
xxx.all<-X
y.all<- matrix(y_true/100, ncol=1)  
n.all <- nrow(y.all)
colnames(xxx.all)
# Identify label ijk for yield
rownames(y.all) <- as.character(1:n.all)

reorder<-order(y.all[,1],decreasing=TRUE)[1:n.all]


y.all<-as.data.frame(y.all[reorder,1], row.names=rownames(y.all)[reorder])
x.all<-x.all[reorder,]
xxx.all<-xxx.all[reorder,]

rownames(y.all)
rownames(x.all)
rownames(xxx.all)

cols123<-X_data[reorder,]



for (ijk in 1:n.all) {
  add.i <- colnames(x.all)[which(x.all[ijk,]!=0)][1]
  add.I <- sum(grepl("alcohol",colnames(x.all)[which(x.all[ijk,]!=0)]))
  if (add.I>1) add.i <- "alcohol 37"  # this is "additive 37"
  
  bas.j <- colnames(x.all)[which(x.all[ijk,]!=0)][1+add.I]
  bas.J <- sum(grepl("base",colnames(x.all)[which(x.all[ijk,]!=0)]))
  if (bas.J>1) bas.j <- "base 4" 
  
  # this is "base 1"
  suf.k <- colnames(x.all)[which(x.all[ijk,]!=0)][1+add.I+bas.J]
  suf.K <- sum(grepl("sulfonyl_fluoride",colnames(x.all)[which(x.all[ijk,]!=0)]))
  if (suf.K>1) suf.k <- "sulfonyl_fluoride 5"                                           
  rownames(y.all)[ijk] <- paste(add.i, bas.j, suf.k, sep=":")
}

rownames(xxx.all) <- rownames(y.all)


#print(xxx.all[500:504,range_main])


###############################
## Factorization PLSGLM   
###############################



# Initialize table
df.yields.colnames <- c("alcohol", "base", "sulfonyl_fluoride", 
                        "y.obs", "mu.hat", "g(mu.hat)", "intercept",
                        "a", "b", "sf",
                        "a:b", "a:sf", "b:sf", 
                        "a:b:sf")
df.yields <- data.frame(matrix(ncol=length(df.yields.colnames),nrow=0))
colnames(df.yields) <- df.yields.colnames

# Get coefficients of best model
coef.plsglm.all <- read.table(file="../Results/coefs_shim_deoxy_table1.txt") ##read coefs as here suitable form!!!!!!
beta.hat.plsglm.all <- coef.plsglm.all[,2]
eta.hat.plsglm.all <- as.vector(as.matrix(xxx.all)%*%beta.hat.plsglm.all[-1]+beta.hat.plsglm.all[1]) ### TAKE CARE TO MAKE xxx.all
mu.hat.plsglm.all <- kappa1(eta.hat.plsglm.all)      

coef_not_interc<-coef.plsglm.all[-1,2]

print(sum(coef_not_interc[range_main]==0)/length(coef_not_interc[range_main]) )
print(sum(coef_not_interc[range_theta]==0)/length(coef_not_interc[range_theta]) )
print(sum(coef_not_interc[range_psi]==0)/length(coef_not_interc[range_psi]) )

print(r2(y.all[,1], mu.hat.plsglm.all))

# Get top yields
n.yield.max <- length(y.all[,1])
IJK <- order(y.all[,1],decreasing=TRUE)[1:n.yield.max]
y.IJK <- as.data.frame(y.all[IJK,1], row.names=rownames(y.all)[IJK])

counter<-0

for (ijk in rownames(y.IJK)) {
  counter<-counter+1
  
  # For the given yield get: relevant features, non-standard features, yield value, yield position in data, yield estimated mean
  tags.ijk <- strsplit(ijk, split=":")[[1]]
  spec.ijk <- tags.ijk[(grepl("alcohol 37",tags.ijk))|
                           (grepl("base 4",tags.ijk))|
                           (grepl("sulfonyl_fluoride 5",tags.ijk))]        
  y.ijk <- format(y.IJK[ijk,], digits=8, nsmall=8, scientific=F)#get y
  ind.ijk <- which(rownames(y.all)==ijk)#get ind
  mu.ijk <- mu.hat.plsglm.all[ind.ijk] #get mu
  
  # Get coefficients with active features (we cross-multiply to get correct signs for positive and negative effects)
  active.ijk <- which(xxx.all[ind.ijk,]!=0)
  coef.plsglm.all.ijk <- rbind(coef.plsglm.all[1,], coef.plsglm.all[-1,][active.ijk,])
  coef.plsglm.all.ijk$V2 <- c(coef.plsglm.all.ijk$V2[1], as.numeric(xxx.all[ind.ijk,active.ijk])*coef.plsglm.all.ijk$V2[-1])
  
  # Deal with one-level combinations
  for (spec in spec.ijk) {
    spec.nonum <- strsplit(spec, split=" ")[[1]][1]
    cut.ijk <- (!grepl(":",coef.plsglm.all.ijk[,1]))&(grepl(spec.nonum,coef.plsglm.all.ijk[,1]))
    coef.plsglm.all.ijk[cut.ijk,]
    new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
    new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
    new.tag <- spec
    coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
    coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
  }
  
  # Deal with two-level combinations
  for (tags in setdiff(tags.ijk, spec.ijk)) {
    for (spec in spec.ijk) {
      spec.nonum <- strsplit(spec, split=" ")[[1]][1]
      cut.ijk <- (grepl("1",lengths(regmatches(coef.plsglm.all.ijk[,1], 
                                                gregexpr(":", coef.plsglm.all.ijk[,1])))))&(grepl(tags,coef.plsglm.all.ijk[,1]))&(grepl(spec.nonum,coef.plsglm.all.ijk[,1]))
      new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
      new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
      new.tag <- paste(sort(c(tags, spec)),collapse=":")
      coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
      coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
    }
  }
  
  for (spec1 in spec.ijk) {
    for (spec2 in setdiff(spec.ijk,c(spec1))) {
      spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
      spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
      cut.ijk <- (grepl("1",lengths(regmatches(coef.plsglm.all.ijk[,1], 
                                                gregexpr(":", coef.plsglm.all.ijk[,1])))))&(grepl(spec1.nonum,coef.plsglm.all.ijk[,1]))&(grepl(spec2.nonum,coef.plsglm.all.ijk[,1]))
      new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
      new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
      new.tag <- paste(sort(c(spec1, spec2)),collapse=":")
      coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
      coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
    }
  }
  
  # Deal with three-level combinations
  for (tags1 in setdiff(tags.ijk, spec.ijk)) {
    for (tags2 in setdiff(setdiff(tags.ijk, spec.ijk),c(tags1))) {
      for (spec in spec.ijk) {
        spec.nonum <- strsplit(spec, split=" ")[[1]][1]
        cut.ijk <- (grepl("2",lengths(regmatches(coef.plsglm.all.ijk[,1], gregexpr(":", coef.plsglm.all.ijk[,1])))))&
          (grepl(tags1,coef.plsglm.all.ijk[,1]))&(grepl(tags2,coef.plsglm.all.ijk[,1]))&(grepl(spec.nonum,coef.plsglm.all.ijk[,1]))
        new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
        new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
        new.tag <- paste(sort(c(tags1, tags2, spec)),collapse=":")
        coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
        coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
      }
    }
  }
  
  for (tags in setdiff(tags.ijk, spec.ijk)) {
    for (spec1 in spec.ijk) {
      for (spec2 in setdiff(spec.ijk,c(spec1))) {
        spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
        spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
        cut.ijk <- (grepl("2",lengths(regmatches(coef.plsglm.all.ijk[,1], gregexpr(":", coef.plsglm.all.ijk[,1])))))&
          (grepl(tags,coef.plsglm.all.ijk[,1]))&(grepl(spec1.nonum,coef.plsglm.all.ijk[,1]))&(grepl(spec2.nonum,coef.plsglm.all.ijk[,1]))
        new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
        new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
        new.tag <- paste(sort(c(tags, spec1, spec2)),collapse=":")
        coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
        coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
      }
    }
  }
  
  for (spec1 in spec.ijk) {
    for (spec2 in setdiff(spec.ijk,c(spec1))) {
      for (spec3 in setdiff(spec.ijk,c(spec1,spec2))) {
        spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
        spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
        spec3.nonum <- strsplit(spec3, split=" ")[[1]][1]
        cut.ijk <- (grepl("2",lengths(regmatches(coef.plsglm.all.ijk[,1], gregexpr(":", coef.plsglm.all.ijk[,1])))))&
          (grepl(spec1.nonum,coef.plsglm.all.ijk[,1]))&(grepl(spec2.nonum,coef.plsglm.all.ijk[,1]))&(grepl(spec3.nonum,coef.plsglm.all.ijk[,1]))
        new.ind <- (1:length(coef.plsglm.all.ijk[,1]))[cut.ijk]
        new.val <- sum(coef.plsglm.all.ijk[new.ind,2])
        new.tag <- paste(sort(c(spec1, spec2, spec3)),collapse=":")
        coef.plsglm.all.ijk <- coef.plsglm.all.ijk[-new.ind,]
        coef.plsglm.all.ijk <- rbind(coef.plsglm.all.ijk, list(new.tag,new.val))
      }
    }
  }
  

  # Recover original levels
  add.ijk <- as.numeric(strsplit(sort(tags.ijk)[1],split=" ")[[1]][2])
  bas.ijk <- as.numeric(strsplit(sort(tags.ijk)[2],split=" ")[[1]][2])
  suf.ijk <- as.numeric(strsplit(sort(tags.ijk)[3],split=" ")[[1]][2])

  print("start")
  add.true <- paste("alcohol ", apply(xs, 2, unique)$alcohol[which(apply(xs, 2, function(x) order(unique(x)))$alcohol==add.ijk)], sep="")
  print(add.true)
  bas.true <- paste("base ", apply(xs, 2, unique)$base[which(apply(xs, 2, function(x) order(unique(x)))$base==bas.ijk)], sep="")
  suf.true <- paste("sulfonyl_fluoride ", apply(xs, 2, unique)$sulfonyl_fluoride[which(apply(xs, 2, function(x) order(unique(x)))$sulfonyl_fluoride==suf.ijk)], sep="")
  coef.plsglm.all.ijk <- coef.plsglm.all.ijk[order(coef.plsglm.all.ijk$V1),]
  coef.plsglm.all.ijk$V2 <- format(as.numeric(coef.plsglm.all.ijk$V2),nsmall=7,digits=7, scientific=F)
  coef.plsglm.all.ijk$V1 <- sort(c("intercept", add.true, bas.true, suf.true, paste(add.true, bas.true, suf.true, sep=":"),
                              paste(add.true, bas.true, sep=":"), paste(add.true, suf.true, sep=":"), paste(bas.true, suf.true, sep=":")
                              ))
  coef.plsglm.all.ijk <- rbind( c("g.link(mu.hat)", format(eta.hat.plsglm.all[ind.ijk],nsmall=7,digits=7, scientific=F)), coef.plsglm.all.ijk)
  coef.plsglm.all.ijk <- rbind( c("mu.hat", format(mu.hat.plsglm.all[ind.ijk],nsmall=7,digits=7, scientific=F)), coef.plsglm.all.ijk)
  coef.plsglm.all.ijk <- rbind( c("y.obs", format(y.all[ind.ijk,],nsmall=7,digits=7, scientific=F)), coef.plsglm.all.ijk)
  #print(coef.plsglm.all.ijk)
  
  # Update global table
  df.yields <- rbind(df.yields, 
                     c(as.character(cols123[counter,1]), as.character(cols123[counter,2]), as.character(cols123[counter,3]),
                       coef.plsglm.all.ijk[1,2], coef.plsglm.all.ijk[2,2], coef.plsglm.all.ijk[3,2], coef.plsglm.all.ijk[10,2], #interc
                       coef.plsglm.all.ijk[4,2], coef.plsglm.all.ijk[8,2], coef.plsglm.all.ijk[11,2], #mains
                       coef.plsglm.all.ijk[5,2], coef.plsglm.all.ijk[7,2], coef.plsglm.all.ijk[9,2], #two-ways
                       coef.plsglm.all.ijk[6,2]# three ways 
                       ))
  colnames(df.yields) <- df.yields.colnames
  
}
gc()

# Save global table
write.table(df.yields, file=paste("../Results/yield_table_deoxy.txt",sep=""), sep=",", row.names=FALSE, col.names=TRUE)


X_data[which(y>97.3),]
