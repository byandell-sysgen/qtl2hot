## without model C
get.power.type1.prec.matrix.2 <- function(out, models, alpha)
{
  n <- length(alpha)
  Power <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Type1 <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Prec <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                         "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                       as.character(alpha)))
  for(k in 1:n){
    outs <- array(NA, c(9,2,5), dimnames=list(c("aic","par.joint.aic","par.aic",
                                                "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                              c("TP","FP"),models))
    for(i in 1:5){
      outs[1,1:2,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="aic")[1:2])
      outs[2,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.aic")[1:2])
      outs[3,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.aic")[1:2])
      outs[4,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.aic")[1:2])
      outs[5,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="bic")[1:2])
      outs[6,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.bic")[1:2])
      outs[7,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.bic")[1:2])
      outs[8,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.bic")[1:2])
      outs[9,,i] <- as.numeric(performance.summaries.cit(out[[i]], model=models[i], alpha[k])[1:2])
    }
    all <- matrix(NA, 9, 2, dimnames=list(c("aic","par.joint.aic","par.aic",
                                            "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                          c("TP","FP")))
    for(i in 1:9)
      all[i,] <- apply(outs[i,1:2,-3],1,sum)
    Prec[,k] <- all[,1]/apply(all,1,sum)
    Power[,k] <- all[,1]/4000
    Type1[,k] <- all[,2]/4000
    print(k)
  }
  list(Power=Power, Type1=Type1, Prec=Prec)
}
