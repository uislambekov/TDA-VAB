library(TDA) 
library(igraph)
library(Rcpp)

# Example 1
# sample 100 points uniformly from unit circle and add Gaussian noise
N <- 100
X <- circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.1)  
plot(X,asp = 1,pch=19,xlab = 'x',ylab = 'y') 

# compute PD using Rips filtration
source("VAB.R") # assuming VAB.R is in the current working directory
D <- ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 

# compute VAB for homological dimensions H0 and H1
computeVAB(D,homDim = 0,scaleSeq = seq(0,2,length.out=11)) 
computeVAB(D,homDim = 1,scaleSeq = seq(0,2,length.out=11)) 

# compare with C++ version
sourceCpp("VAB.cpp") # assuming VAB.cpp is in the current working directory

# compute VAB for homological dimensions H0 and H1
computeVABcpp(D,homDim = 0,scaleSeq = seq(0,2,length.out=11)) 
computeVABcpp(D,homDim = 1,scaleSeq = seq(0,2,length.out=11)) 

# Example 2
nNode <- 100
alpha <-c(1.5,1.5,1.5)
# sample from a Dirichlet distribution
lpvs <- sample_dirichlet(nNode,alpha) 
# generate a random graph according to the random dot product graph model
g <- sample_dot_product(lpvs) 

# construct a simplicial complex of dimension 2
zeroSimplx=as.list(V(g))
oneSimplx=data.frame(t(as_edgelist(g)))
twoSimplx=data.frame(matrix(triangles(g),nrow = 3))
cmplx=c(zeroSimplx,oneSimplx,twoSimplx)

# use cross-entropy to define node attributes
nodeAtr <- colSums(-lpvs*log2(lpvs))
nodeAtr <- nodeAtr/max(nodeAtr) # rescale node attributes to lie in [0,1]

# construct lower-star filtration using node attributes
fltrn <- funFiltration(FUNvalues = nodeAtr,cmplx = cmplx,sublevel = TRUE) 
# compute the PD of the filtration 
D <- filtrationDiag(filtration = fltrn,
                    maxdimension = 1,
                    library = 'Dionysus',
                    location = T,
                    diagLimit =1.01)$diagram

# compute VAB for homological dimensions H0 and H1
computeVAB(D,homDim = 0,scaleSeq = seq(0,1,length.out=11)) 
computeVAB(D,homDim = 1,scaleSeq = seq(0,1,length.out=11)) 
