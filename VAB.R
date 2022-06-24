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

