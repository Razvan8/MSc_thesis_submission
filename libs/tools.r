# Functions can be accessed via the command source("tools.r")


# Functions for the textbox in ggplot

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




# Continuous Bernoulli likelihood

kappa0 <- function(x) {
  k0 <- rep(0,length(x))
  tt <- which((x<=700)&(x!=0))
  k0[tt] <- log((exp(x[tt])-1)/x[tt])
  tt <- which(x>700)
  k0[tt] <- x[tt] - log(x[tt])
  return(k0)
}

kappa1 <- function(x) {
  k1 <- rep(1/2,length(x))
  tt <- which(abs(x)<=0.0001)
  k1[tt] <- 1/2 + x[tt]/12 - x[tt]^3/720 + x[tt]^5/30240
  tt <- which(abs(x)>0.0001)
  k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
  return(k1)
}

kappa2 <- function(x) {
  k2 <- rep(1/12,length(x))
  tt <- which(abs(x)<=0.015)
  k2[tt] <- 1/12 - x[tt]^2/240 + x[tt]^4/6048
  tt <- which(abs(x)>0.015)
  k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
  return(k2)
}

g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN=function(v) min(max(v,0.001),0.999))
  g <- 3.5*tan(pi*(2*tt-1)/2) 
  return(g)
}

V.link <- function(x) {
  V <- kappa2(g.link(x))
  return(V)
}

A.link <- function(x) {
  A <- rep(0,length(x))
  for (tt in 1:length(x)) {
    A[tt] <- integrate(function(x) (V.link(x))^(-1/3), lower=0, upper=x[tt], subdivisions=1000)$value
  }
  return(A)
}




# Fisher scoring

fisher.scoring <- function(y,x,initial){
  
  # Initialization
  beta0 <- initial
  beta1 <- beta0
  
  # Main loop
  x <- cbind(rep(1,nrow(x)),x)
  counter <- 0
  repeat{
    counter <- counter+1
    eta <- as.vector(x%*%beta0)
    
    kp <- kappa1(eta)
    kpp <- kappa2(eta)
    kpp[which(kpp<0.01)] <- 0.01
    
    W <- kpp
    Z <- x%*%beta0+(y-kp)/W
    fit <- lm(Z~x-1,weights=W)
    beta1 <- fit$coef
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    print(paste("Epsilon: ", epsilon, sep=""))
    if(epsilon<=0.01) break
    if(counter==100) {print("no convergence"); break}
    beta0 <- beta1
  }
  list(beta.hat=as.vector(beta1),zstat=summary(fit)[["coefficients"]][, "t value"])
}




# The next function implements Partial Least Squares on
# uni/multivariate response Y and design X using KERNELPLS
# The function allows for centering/scaling (FALSE by default)
# The function allows for intercept (FALSE by default)

kernelpls2 <- function(X, Y, ncomp, 
                       centering=FALSE, scaling=FALSE, intercept=FALSE,
                       maxit=50, tol=1e-8,
                       verbose=FALSE){
  
  if (verbose) print("Performing KERNELPLS...")
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), tol)})
    sd.Y <- apply(Y, 2, function(v){max(sd(v), tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), tol)})
    Y <- apply(Y, 2, function(v){v/max(sd(v), tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X) # intercept is added as first column
  }    
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  ## Initialization
  R <- matrix(0, ncol=ncomp, nrow=p)        # projection
  P <- matrix(0, ncol=ncomp, nrow=p)        # X loadings
  tQ <- matrix(0, ncol=m, nrow=ncomp)       # Y loadings; transposed
  B <- array(0, c(p, m, ncomp))             # coefficients
  W <- P                                    # weights
  U <- matrix(0, ncol=ncomp, nrow=n)        # scores
  T <- matrix(0, ncol=ncomp, nrow=n)        # scores
  tsqs <- rep.int(1, ncomp)                 # t't
  
  ## 1.
  XtY <- crossprod(X, Y)
  
  for (a in 1:ncomp) {
    
    ## 2.
    if (m == 1) {
      
      w.a <- XtY / sqrt(c(crossprod(XtY)))
      
    } else {
      
      if (m < p) {
        
        q <- eigen(crossprod(XtY), symmetric = TRUE)$vectors[,1]
        w.a <- XtY %*% q
        w.a <- w.a / sqrt(c(crossprod(w.a)))
        
      } else {
        
        w.a <- eigen(XtY %*% t(XtY), symmetric = TRUE)$vectors[,1]
        
      }
    }
    
    ## 3.
    r.a <- w.a
    
    if (a > 5) {
      
      ## This is faster when a > 5:
      r.a <- r.a - colSums(crossprod(w.a, P[,1:(a-1), drop=FALSE])%*%t(R[,1:(a-1), drop=FALSE]))
      
    } else if (a > 1) {
      
      for (j in 1:(a - 1)) {
        r.a <- r.a - c(P[,j] %*% w.a) * R[,j]
      }
      
    }
    
    ## 4.
    t.a <- X %*% r.a
    tsq <- c(crossprod(t.a))
    p.a <- crossprod(X, t.a) / tsq
    q.a <- crossprod(XtY, r.a) / tsq
    
    ## 5.
    XtY <- XtY - (tsq * p.a) %*% t(q.a)
    
    ## 6. 7. 8.
    R[,a]  <- r.a
    P[,a]  <- p.a
    tQ[a,] <- q.a
    B[,,a] <- R[,1:a, drop=FALSE] %*% tQ[1:a,, drop=FALSE]
    tsqs[a] <- tsq
    
    ## Extra step to calculate Y scores:
    u.a <- Y %*% q.a / c(crossprod(q.a))
    
    ## Make u orth to previous X scores:
    if (a > 1) {
      u.a <- u.a - T %*% (crossprod(T, u.a) / tsqs)
    }
    U[,a] <- u.a
    T[,a] <- t.a
    W[,a] <- w.a
    
  }
  
  # Compute coeffs for original variables
  B.tilde <- array(0, c(p, m, ncomp))
  
  if (intercept) {
    
    for (nc in 1:ncomp) {
      for (mc in 1:m) {
        B.tilde[1,mc,nc] <- sd.Y[mc]*B[1,mc,nc] + mu.Y[mc] - sd.Y[mc]*(mu.X/sd.X)%*%B[2:p,mc,nc]  # intercept
        B.tilde[2:p,mc,nc] <- sd.Y[mc]*B[2:p,mc,nc]/sd.X                                          # core
      }
    }
    
  } else {
    
    for (nc in 1:ncomp) {
      for (mc in 1:m) {
        B.tilde[1:p,mc,nc] <- sd.Y[mc]*B[1:p,mc,nc]/sd.X  # no intercept, core
      }
    }
    
  }
  
  return(list(B=B,      # coeffs for centered/scaled variables
              T=T,      # scores for centered/scaled variables
              P=P,      # X loadings for centered/scaled variables
              W=W,      # weights for centered/scaled variables
              U=U,      # scores for centered/scaled variables
              Q=t(tQ),  # Y loadings for centered/scaled variables
              R=R,      # projection for centered/scaled variables
              muX=mu.X, muY=mu.Y,  # vectors of column-wise mean
              sdX=sd.X, sdY=sd.Y,  # vectors of column-wise sd
              Btilde=B.tilde       # coeffs for unscaled variables
              )
         )
  
}



# The next function implements Partial Least Squares in any case
# For the parameter method, "kernel" and "nipals" call the respective algorithms
# The function allows for centering/scaling (kernel only)
# The function allows for intercept (kernel only)

plslm <- function(X, Y, ncomp,
                  centering=FALSE, scaling=FALSE, intercept=FALSE,
                  maxit=50, tol=1e-8, 
                  method="kernel",
                  verbose=FALSE){
  
  if (method=="kernel") {
    if (verbose) print("Performing PLS with kernelpls...")
    return(kernelpls2(X=X, Y=Y, ncomp=ncomp, 
                      centering=centering, scaling=scaling, intercept=intercept, 
                      maxit=maxit, tol=tol,
                      verbose=verbose))
  }
  else {
    if (verbose) print("Invalid method selection, performing kernelpls as default...")
    return(kernelpls2(X=X, Y=Y, ncomp=ncomp, 
                      centering=centering, scaling=scaling, intercept=intercept, 
                      maxit=maxit, tol=tol,
                      verbose=verbose))
  }
  
}



# The next function implements PLSGLM with Continuous Bernoulli response
# The parameter beta0 is the initialization, choosing NULL initializes at zero 
# The function allows for centering/scaling (TRUE by default)
# The function allows for intercept (TRUE by default)

plsglm.cb <- function(X, Y, ncomp, beta0=NULL,
                      centering=TRUE, scaling=TRUE, intercept=TRUE,
                      maxit=10, tol=0.0545,
                      verbose=FALSE, clip=0.01){
  
  if (verbose) print("Performing PLSGLM-NEW")
  
  # CB Likelihood
  kappa1 <- function(x) {
    k1 <- rep(1/2,length(x))
    tt <- which(x!=0)
    k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
    return(k1)
  }
  
  kappa2 <- function(x) {
    k2 <- rep(1/12,length(x))
    tt <- which(x!=0)
    k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
    return(k2)
  }
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    #mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    #Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), tol)})
    #sd.Y <- apply(Y, 2, function(v){max(sd(v), tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), tol)})
    #Y <- apply(Y, 2, function(v){v/max(sd(v), tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X)
  }
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  # Initialization
  if (is.null(beta0)){ beta0 <- rep(0, p) }
  
  beta.old <- array(beta0, dim=c(p, 1, ncomp))
  Z.old <- array(0, dim=c(n, 1, ncomp))
  W.old <- array(0, dim=c(n, n, ncomp))
  
  beta <- array(beta0, dim=c(p, 1, ncomp))
  Z <- array(0, dim=c(n, 1, ncomp))
  W <- array(0, dim=c(n, n, ncomp))
  
  R <- list()
  
  nc <- ncomp
  it.stop <- rep(0, ncomp)
  
  # Run until convergence or stop
  counter <- 0
  
  repeat{
    if (verbose) print(paste("Performing PLSGLM iter", counter, sep=" "))
    
    # Keep old variables
    beta.old <- beta
    Z.old <- Z
    W.old <- W
    
    # Compute W and Z
    eta <- array(0, dim=c(n, 1, ncomp))
    k.p <- array(0, dim=c(n, 1, ncomp))
    k.pp <- array(0, dim=c(n, 1, ncomp))
    
    if (length(nc)==0) {
      if (verbose) print("No comps left...")
      break
    }
    
    if (verbose) print(paste("Comps left", paste(nc, collapse=" ")))
    for (m in nc) {
      
      eta[,,m] <- X%*%beta[,,m]
      
      k.p[,,m] <- kappa1(eta[,,m])
      
      k.pp[,,m] <- kappa2(eta[,,m])
      k.pp[,,m][which(k.pp[,,m]<clip)] <- clip
      
      W[,,m] <- diag(as.vector(k.pp[,,m]))
      
      Z[,,m] <- eta[,,m] + diag(1/as.vector(k.pp[,,m]))%*%(Y-as.vector(k.p[,,m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[,,m])%*%X)
      Z.W <- as.matrix(sqrt(W[,,m])%*%Z[,,m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept X.W
      fit.pls <- plslm(X=X.W, Y=Z.W, ncomp=m,
                       centering=FALSE, scaling=FALSE, intercept=FALSE)
      R[[m]] <- fit.pls$R
      beta[,,m] <- fit.pls$Btilde[,,m]
      
    }
    
    #print(paste("Norm of beta", paste(apply(beta^2, 3, sum), collapse=" "), sep=" "))
    
    epsilon <- sqrt(apply((beta-beta.old)^2, 3, sum)/apply((beta.old)^2, 3, sum) )
    #print(paste("Divergence", paste(epsilon, collapse=" "), sep=" "))
    print(paste("Min Divergence", min(epsilon[nc]), sep=" "))
    
    log.like <- apply(beta, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    #print(paste("Max Loglike", max(log.like[nc]), sep=" "))
    
    log.like.ratio <- log.like - apply(beta.old, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    print(paste("Min Loglike ratio", min(log.like.ratio[nc]), sep=" "))
    
    if (sum(is.nan(epsilon[nc]))>0) {
      nan.stop <- which(is.nan(epsilon))
      if (verbose) print(paste("Divergence NaN comps", paste(nan.stop, collapse=" ")))
      for (m in nc) {
        beta[,,m] <- beta.old[,,m]
        Z[,,m] <- Z.old[,,m]
        W[,,m] <- W.old[,,m]
      }
      nc <- setdiff(nc, nan.stop)
    }
    
    if ((min(epsilon[nc])<tol)|(min(log.like.ratio[nc])<tol)) { 
      nc.stop <- which((epsilon<tol)|(log.like.ratio<tol))
      it.stop[nc.stop] <- counter
      if (verbose) print(paste("Divergence/Loglike stop comps", paste(nc.stop, collapse=" ")))
      nc <- setdiff(nc, nc.stop)
    }
    
    if (counter==maxit) { 
      if (verbose) print("Maximum iterarion, no convergence...")
      it.stop[which(it.stop==0)] <- counter
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta*0
  
  for (m in 1:ncomp) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta[,,m][2:p]/sd.X                                         # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta[,,m][1:p]/sd.X                                       # no intercept, core
      }
    }
  }

  return(list(BETA=beta.tilde, 
              beta=beta, 
              Z=Z, W=W, 
              R=R,
              it=it.stop))
}



# Simplified implementation of PLSGLM-CB

plsglm.cb.simple <- function(x, y, m.plsglm, 
                             beta0=NULL,
                             centering=TRUE, scaling=TRUE, intercept=TRUE,
                             maxit=10, tol=0.0545,
                             verbose=FALSE){
  
  return(plsglm.cb(X=x, Y=y, ncomp=m.plsglm, 
                   beta0=beta0,
                   centering=centering, scaling=scaling, intercept=intercept,
                   maxit=maxit, tol=tol,
                   verbose=verbose)$BETA[,,m.plsglm])
}



# The next function computes DoF for PLS

dof.pls <- function(X, Y, m.max, Lambda){
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  A <- t(X)%*%X
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Initialize matrix of weights
  K <- matrix(0, nrow = p, ncol = m.max)
  
  # Run NIPALS algorithm
  X <- as.matrix(X)
  Y <- as.matrix(y.all)
  w <- matrix(0, nrow = p, ncol = 1)
  s <- matrix(0, nrow = n, ncol = 1)
  c <- 0
  
  for (i in 1:m.max){
    #print(paste("Checking",i))
    w <- t(X)%*%Y
    if (sqrt(sum(w^2))==0){
      break
    }
    w <- w/sqrt(sum(w^2))
    K[1:p, i] <- w
    s <- X%*%w
    P <- (t(X)%*%s)/drop(t(s)%*%s)
    X <- X - s%*%t(P)
    c <- drop(t(s)%*%Y)/drop(t(s)%*%s)
    Y <- Y - c*s
  }
  
  D <- t(K)%*%A%*%K
  
  # Compute Ritz values
  Theta <- matrix(0, nrow = m.max, ncol = m.max)
  for (i in 1:m.max){
    Theta[1:i,i] <- eigen(D[1:i,1:i])$values
  }
  
  # Compute shrinkage factors
  Z <- matrix(0, nrow = p, ncol = m.max)
  product <- 1
  for (m in 1:m.max){
    for (i in 1:p){
      for (j in 1:m){
        #product <- product*(1 - Lambda[i]/Theta[j,m])
        product <- product*(1 - Lambda[i]/Lambda[j])
      }
      Z[i,m] <- 1 - product
      product = 1
    }
  }
  
  # Compute DoF
  DOF <- apply(Z, 2, function(v) min(abs(sum(v)),p))
  
  return(list(DOF=DOF,Z=Z,Theta=Theta))
}



