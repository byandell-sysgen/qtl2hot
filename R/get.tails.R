## Get the n.quant highest values above lod=4 (see scan.perm.R)
get.tails <- function(highs, n.quant = 2000, s.quant = seq(n.quant))
{
  ## Limit n.quant to range of data.
  n.quant <- min(n.quant, max(s.quant), max(table(highs[,"row"])))
  s.quant <- s.quant[s.quant <= n.quant]
  
  tmpfn <- function(x, sn) {
    x <- sort(x, decreasing = TRUE)
    x[sn]
  }
  
  out <- tapply(highs[,"lod"], highs[,"row"], tmpfn, s.quant)
  highrow.names <- names(out)
  
  ## Turn list of items of length n.quant into matrix. This step takes time!
  out <- matrix(unlist(out), n.quant)
  dimnames(out) <- list(s.quant, highrow.names)

  out
}
