nmds <- function(dis,k=2,y=cmdscale(dis,k),maxit=50)
{
  tmp <- isoMDS(dis,y=y,k=k,maxit=maxit)
  class(tmp) <- "nmds"
  return(tmp)
}
