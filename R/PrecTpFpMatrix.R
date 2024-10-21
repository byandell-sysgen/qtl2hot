#' Determine false positive and true positive rates for known targets.
#' 
#' Determine how well different tests do to predict candidates of regulation.
#' 
#' \code{FitAllTests} invokes 7 tests. The hidden routine \code{CitTests} is
#' invoked by call to \code{FitAllTests}; this is hidden because we do not
#' recommend its use.
#' 
#' \code{JoinTestOutputs} joins results of \code{\link{FitAllTests}}, either
#' from a list \code{tests} or from a collection of files prefixed by
#' \code{file}. The joined \code{tests} from \code{JoinTestOutputs} are
#' summarized with \code{PrecTpFpMatrix} using the biologically validated true
#' positives, false positives and precision, for the inferred causal relations.
#' We define a true positive as a statistically significant causal relation
#' between a gene and a putative target gene when the putative target gene
#' belongs to the known signature of the gene. Similarly, we define a false
#' positive as a statistically significant causal relation between a gene and a
#' putative target gene when the target gene does not belong to the signature.
#' (For the AIC and BIC methods that do not provide a p-value measuring the
#' significance of the causal call, we simply use the detected causal relations
#' in the computation of true and false positives). The validated precision is
#' computed as the ratio of true positives by the sum of true and false
#' positives. The \code{PrecTpFpMatrix} computes these measures to both all
#' genes, and to cis genes only. Simulations suggest only non-parametric tests
#' need to be adjusted using Benjamini-Hochberg via \code{p.adjust.np}.
#' 
#' @aliases FitAllTests JoinTestOutputs PrecTpFpMatrix CitTests p.adjust.np
#' @param cross object of class \code{cross}
#' @param pheno1 first phenotype column number or character string name
#' @param pheno2 second phenotype column number or character string name; if
#' more than one, then all phenotypes will be tested against \code{pheno1}
#' @param Q.chr QTL chromosome (number or label)
#' @param Q.pos QTL position in cM
#' @param verbose verbose printout if \code{TRUE}
#' @param comap list result of \code{GetComappingTraits}
#' @param alpha significance levels at which summaries are computed
#' @param val.targets validated targets of candidate regulators
#' @param all.orfs all trait names
#' @param tests list object as list of \code{FitAllTests} results, or of joined
#' output created by \code{JoinTestsOutputs}
#' @param file prefix for file names when running \code{FitAllTests} in
#' parallel and saving test results in separate files
#' @param cand.reg object from \code{\link{GetCandReg}}
#' @param cis.cand.reg object from \code{\link{GetCisCandReg}}
#' @param method method for p-value adjustment; see
#' \code{\link[stats]{p.adjust}}
#' @return List containing \item{Prec1,Prec2}{matrix of precision with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only} \item{Tp1,Tp2}{matrix of true positive rate with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only} \item{Fp1,Fp2}{matrix of false positive rate with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only}
#' @author Elias Chaibub Neto
#' @seealso \code{\link{GetCandReg}}, \code{\link{CMSTtests}},
#' \code{\link[stats]{p.adjust}}
#' @keywords utilities
#' @examples
#' 
#' example(GetCandReg)
#' ## Suppose y1 is causal with targets y2 and y3.
#' targets <- list(y1 = c("y2","y3"))
#' 
#' tests <- list()
#' for(k in seq(names(comap.targets))) {
#'   tests[[k]] <- FitAllTests(CMSTCross, pheno1 = names(comap.targets)[k],
#'                       pheno2 = comap.targets[[k]],
#'                       Q.chr = cand.reg[k, 4],
#'                       Q.pos = cand.reg[k, 5])
#' }
#' names(tests) <- names(comap.targets)
#' tests <- JoinTestOutputs(comap.targets, tests)
#' 
#' PrecTpFpMatrix(alpha = seq(0.01, 0.10, by = 0.01),
#'   val.targets = targets, all.orfs = CMSThigh$names, tests = tests,
#'   cand.reg = cand.reg, cis.cand.reg = cis.cand.reg)
#' 
#' @export PrecTpFpMatrix
PrecTpFpMatrix <- function(alpha, val.targets, all.orfs, tests, cand.reg, cis.cand.reg)
{
  nms <- as.character(cand.reg[,1])
  cis.index <- attr(cis.cand.reg, "cis.index")
  
  le <- length(alpha)
  Prec1 <- Tp1 <- Fp1 <- matrix(NA, 9, le, 
                                dimnames=list(c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", 
                                                "p.aic", "np.aic", "cit"), as.character(alpha)))
  Prec2 <- Tp2 <- Fp2 <- Prec1
  for(i in 1:le){
    aux <- PerformanceSummariesKo(alpha = alpha[i], nms,
                                  val.targets = val.targets, 
                                  all.orfs = all.orfs, 
                                  tests = tests,
                                  cis.index = cis.index)
    Prec1[,i] <- round(aux[[1]][,1], 2) 
    Prec2[,i] <- round(aux[[2]][,1], 2)
    Tp1[,i] <- aux[[1]][, 2]
    Tp2[,i] <- aux[[2]][, 2]
    Fp1[,i] <- aux[[1]][, 3]
    Fp2[,i] <- aux[[2]][, 3]
  }
  list(Prec1 = Prec1,
       Prec2 = Prec2,
       Tp1 = Tp1,
       Tp2 = Tp2,
       Fp1 = Fp1,
       Fp2 = Fp2)
}
