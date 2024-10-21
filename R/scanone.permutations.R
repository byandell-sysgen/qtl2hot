## Performs and saves qtl::nind on permuted dataset
scanone.permutations <- function(cross, pheno.col = seq(3, qtl::nphe(cross)),
                                 n.perm, seed=123456789, batch.effect = NULL,
                                 pheno.set = 1,
                                 lod.min, drop.lod = 1.5,
                                 addcovar = NULL, intcovar = NULL, ...)
{
  set.seed(seed[[1]])
  
  if(!is.null(batch.effect)) {
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))
    covars <- sexbatch.covar(cross, batch.effect)
  }
  else
    covars <- list(addcovar = addcovar, intcovar = intcovar)
  
  n.ind <- qtl::nind(cross)
  
  perms <- matrix(NA, n.ind, n.perm)
  
  for(i in 1:n.perm){
    perm.cross <- cross
    perms[,i] <- tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]
    
    per.scan <- qtl::scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                             addcovar=covars$addcovar, intcovar=covars$intcovar, ...)
    
    per.scan.hl <- highlod(per.scan, lod.thr = lod.min, drop.lod = drop.lod,
                           restrict.lod = TRUE)$highlod
    
    save(per.scan.hl, perms,
         file=paste("per.scan",pheno.set, i,"RData",sep="."))
  }
}
