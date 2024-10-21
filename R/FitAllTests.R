FitAllTests <- function(cross, pheno1, pheno2, Q.chr, Q.pos, verbose = TRUE)
{
  out <- CMSTtests(cross, pheno1, pheno2, Q.chr, Q.pos, 
                   NULL, NULL, NULL, NULL, "all", "both", verbose)
  
  nms <- pheno2
  ntests <- length(pheno2)
  out$pvals.cit <- matrix(NA, ntests, 2, dimnames = list(nms, c("pval.1", "pval.2")))
  
  for(k in 1 : ntests) {
    cit.mar <- qtl::find.marker(cross, Q.chr, Q.pos)
    LL <- qtl::pull.geno(cross)[, cit.mar]
    GG <- cross$pheno[, pheno1]
    TT <- cross$pheno[, pheno2[k]]
    aux2 <- try(CitTests(LL, GG, TT), silent = TRUE)
    if(class(aux2) != "try-error") {
      out$pvals.cit[k,] <- aux2
    }
    if(verbose)
      cat("CIT pheno2 = ", pheno2[k], "\n")   
  }
  out
}
