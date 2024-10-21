JoinTestOutputs <- function(comap, tests, file = NULL)
{
  if(!is.null(file) & missing(tests))
    tests <- CombineTests(comap, file)
  
  reg.nms <- names(comap)
  join.out <- tests[[1]]
  ## Add extra element to join.out: phenos.
  join.out$phenos <- cbind(rep(reg.nms[1], length(comap[[1]])), comap[[1]])
  
  for (k in 2 : length(comap)) {
    out <- tests[[k]]
    if (length(out) != 10) { ## CMSTtests if length(pheno2) == 1
      tmp <- t(out$Z.aic)
      AIC.stats <- c(out$AICs, tmp[!is.na(tmp)])
      tmp <- t(out$Z.bic)
      BIC.stats <- c(out$BICs, tmp[!is.na(tmp)])      
      out <- list(R2s = out$R2, AIC.stats = AIC.stats, BIC.stats = BIC.stats, 
                  pvals.j.BIC = out$pvals.j.BIC, pvals.p.BIC = out$pvals.p.BIC, 
                  pvals.np.BIC = out$pvals.np.BIC,  pvals.j.AIC = out$pvals.j.AIC, 
                  pvals.p.AIC = out$pvals.p.AIC,  pvals.np.AIC = out$pvals.np.AIC,
                  pvals.cit = out$pvals.cit)
    }
    
    for (i in 1 : 10) {
      join.out[[i]] <- rbind(join.out[[i]], out[[i]])
    }
    join.out[[11]] <- 
      rbind(join.out[[11]], cbind(rep(reg.nms[k], length(comap[[k]])), comap[[k]]))
  }
  join.out
}
