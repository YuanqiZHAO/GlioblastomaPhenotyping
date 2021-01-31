affinityMatrix <- function(Diff,K=20,sigma=0.5) {
###This function constructs similarity networks.
  N = nrow(Diff)
  
  Diff = (Diff + t(Diff)) / 2
  diag(Diff) = 0;

  sortedColumns = as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }

  means = apply(sortedColumns[,1:K+1,drop=F],1,finiteMean)+.Machine$double.eps;
  # print(c('means', class(means),dim(means)))
  # print(means)
  # print(c('Diff', class(Diff),dim(Diff)))
  # print(Diff)
  avg <- function(x,y) ((x+y)/2)
  Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
  densities = dnorm(Diff,0,sigma*Sig,log = FALSE)
  
  W = (densities + t(densities)) / 2
  
  return(W)
}
