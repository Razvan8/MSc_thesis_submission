libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
data<-load_deoxi_flourination()

print(max(data$yield))
hist(data$prob)

which(data$yield_corrected==100)
which(data$yield>100)

data


plot_yield_histo(data)


levels(data$a)
levels(data$b)
levels(data$s)

y_true<-data$prob*100
X_data<-data[c('a','b','s')]

X_data <- lapply(X_data, as.factor)
X_data <- as.data.frame(X_data)  # Convert back to data frame


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

new_colnames <- c(paste0("a.", 1:36), paste0("b.", 1:3), paste0("s.", 1:4))
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

X<-X_dummy
y<-y_true


################################################################ USE LASSO AND MY CLASS ON DATA #################################





