CitTests <- function(LL, GG, TT)
{
  no.bootstrap <- 50
  ### remove missing values ###
  sel <- (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  dat_f <- as.data.frame(cbind(LL, GG, TT), stringsAsFactors = FALSE)
  dat_f <- dat_f[sel,]
  names(dat_f) <- c("L", "G", "T")
  Lf <- as.factor(dat_f$L)
  dat_f$L <- as.integer(Lf) - 1
  llevels <- as.integer(unique(as.factor(dat_f$L)))
  dfL <- length(llevels) - 1
  pvec <- rep(NA, 4)
  
  if(dfL == 2){
    dat_f$L1 <- ifelse(dat_f$L == 1,1,0)
    dat_f$L2 <- ifelse(dat_f$L == 2,1,0)
    fit0 <- stats::lm(T ~ 1, data = dat_f)
    fit1 <- stats::lm(T ~ L1 + L2, data = dat_f)
    fit2 <- stats::lm(G ~ T, data = dat_f)
    fit3 <- stats::lm(T ~ G, data = dat_f)
    fit4 <- stats::lm(G ~ T + L1 + L2, data = dat_f)
    fit5 <- stats::lm(T ~ G + L1 + L2, data = dat_f)
    pvec[1] <- stats::anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- stats::anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- stats::anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- stats::anova(fit3, fit5)$F[2]
    fit1G <- stats::lm(G ~ L1 + L2, data = dat_f)
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    blg2 <- summary(fit1G)$coefficients["L2", 1]
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    blt2 <- summary(fit1)$coefficients["L2", 1]
    dat_f$eG <- stats::resid(fit1G)
    dat_f$eT <- stats::resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss * stats::runif(ss, 0, 1)) ;
      dat_f$G_ <- alg + blg1 * dat_f$L1 + blg2 * dat_f$L2 + dat_f$eG[nni]
      fit_0 <- stats::lm(T ~ G_, data = dat_f)
      fit_1 <- stats::lm(T ~ G_ + L1 + L2, data = dat_f)
      fvecr[i] <- stats::anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + blt2 * dat_f$L2 + dat_f$eT[nni]
      fit_0 <- stats::lm(G ~ T_, data = dat_f)
      fit_1 <- stats::lm(G ~ T_ + L1 + L2, data = dat_f)
      fvecr_r[i] <- stats::anova(fit_0, fit_1)$F[2]
    }
  }#End dfL == 2
  
  if(dfL == 1){
    dat_f$L1 <- ifelse(dat_f$L == 1, 1, 0)
    fit0 <- stats::lm(T ~ 1, data = dat_f)
    fit1 <- stats::lm(T ~ L1, data = dat_f)
    fit2 <- stats::lm(G ~ T, data = dat_f)
    fit3 <- stats::lm(T ~ G, data = dat_f)
    fit4 <- stats::lm(G ~ T + L1, data = dat_f)
    fit5 <- stats::lm(T ~ G + L1, data = dat_f)
    pvec[1] <- stats::anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- stats::anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- stats::anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- stats::anova(fit3, fit5)$F[2]
    fit1G <- stats::lm(G ~ L1, data = dat_f)
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    dat_f$eG <- stats::resid(fit1G)
    dat_f$eT <- stats::resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss * stats::runif(ss, 0, 1)) 
      dat_f$G_ <- alg + blg1 * dat_f$L1 + dat_f$eG[nni]
      fit_0 <- stats::lm(T ~ G_, data = dat_f)
      fit_1 <- stats::lm(T ~ G_ + L1, data = dat_f)
      fvecr[i] <- stats::anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + dat_f$eT[nni]
      fit_0 <- stats::lm(G ~ T_, data = dat_f)
      fit_1 <- stats::lm(G ~ T_ + L1, data = dat_f)
      fvecr_r[i] <- stats::anova(fit_0, fit_1)$F[2]
    }
  } #End dfL == 1
  
  #####F Method
  fvecr <- fvecr[!is.na(fvecr)]
  df1 <- stats::anova(fit3, fit5)$Df[2]
  df2 <- stats::anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- stats::pf(fvecr, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- stats::qnorm(npvals)
  npf <- stats::pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- stats::qnorm(npf)
  pvec[4] <- stats::pnorm(zf, mean = 0, sd = stats::sd(nfvecr))
  pvalc <- max(pvec)  ###Causal p-value
  
  #### Reactive p-value
  fit0G <- stats::lm(G ~ 1, data = dat_f)
  pvec1 <- rep(NA, 4)
  pvec1[1] <- stats::anova(fit0G, fit1G)$"Pr(>F)"[2]
  pvec1[2] <- stats::anova(fit3, fit5)$"Pr(>F)"[2]
  pvec1[3] <- stats::anova(fit1G, fit4)$"Pr(>F)"[2]
  f_ <- stats::anova(fit2, fit4)$F[2]
  #####F Method
  fvecr_r <- fvecr_r[!is.na(fvecr_r)]
  df1 <- stats::anova(fit3, fit5)$Df[2]
  df2 <- stats::anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr_r, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- stats::pf(fvecr_r, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- stats::qnorm(npvals)
  npf <- stats::pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- stats::qnorm(npf)
  pvec1[4] <- stats::pnorm(zf, mean = 0, sd = stats::sd(nfvecr))
  pvalr <- max(pvec1)  ###Reactive p-value
  ###
  c(pvalc, pvalr)
}
