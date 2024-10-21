PerformanceSummariesKo <- function(alpha, nms, val.targets, all.orfs, 
                                   tests, cis.index)
{
  ## Unclear what this does. Part of KO data analysis.
  ## tests added.
  
  rnms <- c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", "p.aic", "np.aic", "cit")
  nt <- length(nms)
  TP <- FP <- TN <- FN <- NC <- data.frame(matrix(0, 9, nt))
  tar <- rep(NA, nt)
  names(TP) <- names(FP) <- names(TN) <- names(FN) <- names(NC) <- nms
  row.names(TP) <- row.names(FP) <- row.names(TN) <- row.names(FN) <- row.names(NC) <- rnms
  Causal <- NotCausal <- vector(mode = "list", length = 9)
  length.intersect <- function(a, b) {
    ints <- intersect(a, b)
    if(is.null(ints))
      0
    else
      length(!is.na(ints))
  }
  for (k in 1 : nt) {
    aux <- which(tests[[11]][,1] == nms[k])
    aux.nms <- tests[[11]][aux, 2]
    tar[k] <- length(aux)
    for (i in 2 : 3) {
      aux.rank <- apply(tests[[i]][aux, 1:4, drop = F], 1, rank)
      aux.best <- apply(aux.rank, 2, function(x) which.min(x))
      aux.index <- as.numeric(which(aux.best == 1))
      Causal[[i - 1]] <- aux.nms[aux.index] 
      NotCausal[[i - 1]] <- aux.nms[-aux.index] 
    }  
    for (i in 4 : 9) {
      Causal[[i - 1]] <- aux.nms[which(tests[[i]][aux, 1] <= alpha)]
      NotCausal[[i - 1]] <- 
        aux.nms[c(which(tests[[i]][aux, 2] <= alpha),
                  which(tests[[i]][aux, 3] <= alpha),
                  which(tests[[i]][aux, 4] <= alpha))] 
    }  
    Causal[[9]] <- 
      aux.nms[which(tests[[10]][aux, 1] <= alpha & tests[[10]][aux, 2] > alpha)]
    NotCausal[[9]] <- 
      aux.nms[c(which(tests[[10]][aux, 1] > alpha & tests[[10]][aux, 2] <= alpha),
                which(tests[[10]][aux, 1] >= alpha & tests[[10]][aux, 2] >= alpha))]
    val <- val.targets[[match(nms[k], names(val.targets))]]
    not.val <- all.orfs[-match(unique(c(as.character(nms[k]), val)), all.orfs)]
    for (i in 1 : 9) {
      TP[i, k] <- length.intersect(Causal[[i]], val)
      FP[i, k] <- length.intersect(Causal[[i]], not.val)
      TN[i, k] <- length.intersect(NotCausal[[i]], not.val)
      FN[i, k] <- length.intersect(NotCausal[[i]], val)
    }
    for (i in 4 : 9) {
      NC[i - 1, k] <- length(c(which(tests[[i]][aux, 1] > alpha),
                               which(tests[[i]][aux, 2] > alpha),
                               which(tests[[i]][aux, 3] > alpha),
                               which(tests[[i]][aux, 4] > alpha)))
    }
    NC[9, k] <- length(c(which(tests[[10]][aux, 1] < alpha),
                         which(tests[[10]][aux, 2] < alpha)))
  }
  tp <- apply(TP, 1, sum)
  fp <- apply(FP, 1, sum)
  tn <- apply(TN, 1, sum)
  fn <- apply(FN, 1, sum)
  nc <- apply(NC, 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.1 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  tp <- apply(TP[, cis.index, drop = FALSE], 1, sum)
  fp <- apply(FP[, cis.index, drop = FALSE], 1, sum)
  tn <- apply(TN[, cis.index, drop = FALSE], 1, sum)
  fn <- apply(FN[, cis.index, drop = FALSE], 1, sum)
  nc <- apply(NC[, cis.index, drop = FALSE], 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.2 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  list(overall.1, overall.2, tar)
}
