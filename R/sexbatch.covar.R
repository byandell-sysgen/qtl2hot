sexbatch.covar <- function(cross, batch.effect, verbose = FALSE)
{
  ic <- qtl::getsex(cross)$sex
  ## Drop sex if only one present.
  if(length(unique(ic)) == 1)
    ic <- NULL
  
  if(!is.null(batch.effect)){
    batch <- cross$pheno[,batch.effect, drop = FALSE]
    tmp <- stats::formula(paste("~ factor(", batch.effect, ")"))
    if(verbose)
      cat("sexbatch.covar", names(tmp), unique(factor(batch[[1]])), "\n")
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    batch <- stats::model.matrix(tmp,batch)[,-1, drop = FALSE]
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    ac <- cbind(batch,ic)
  }
  else
    ac <- ic
  
  list(addcovar = ac, intcovar = ic)
}
