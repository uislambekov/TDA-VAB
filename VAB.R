library(TDA) # load TDA package

computeVAB=function(D,homDim,scaleSeq){
  # D - N by 3 matrix (columns contain dimension, birth and death values respectively)
  # homDim - homological dimension (0 for H0, 1 for H1, etc.)
  # scaleSeq - a sequence of scale values for vectorization
  D <- D[D[,1]==homDim,2:3,drop=F] # extracting PD of dimension = homDim
  delta <- diff(scaleSeq) # deltas 
  l <- length(delta)
  vab <- numeric(length = l)
  for (k in 1:l){
    b <- pmin(scaleSeq[k+1],D[,2])-pmax(scaleSeq[k],D[,1])
    vab[k] <- sum(pmax(0,b))/delta[k] 
  }
  vab # returned object
}

# Example
# sample 100 points uniformly from unit circle and add Gaussian noise
X <- circleUnif(100) + rnorm(200,mean = 0,sd = 0.1)  
plot(X,asp = 1,pch=19,xlab = 'x',ylab = 'y') 

# compute PD using Rips filtration
D <- ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 

# compute VAB for homological dimensions H0 and H1
computeVAB(D,homDim = 0,scaleSeq = seq(0,2,length.out=11)) 
computeVAB(D,homDim = 1,scaleSeq = seq(0,2,length.out=11)) 

# compare with C++ version
library(Rcpp) # load Rcpp package
sourceCpp("VAB.cpp") # assuming VAB.cpp is in the current working directory

# compute VAB for homological dimensions H0 and H1
computeVABcpp(D,homDim = 0,scaleSeq = seq(0,2,length.out=11)) 
computeVABcpp(D,homDim = 1,scaleSeq = seq(0,2,length.out=11)) 
