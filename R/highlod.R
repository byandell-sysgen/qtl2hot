######################################################################
# highlod.R
#
# Elias Chaibub Neto
# Brian S Yandell
# Aimee Teo Broman
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: highlod, sexbatch.covar,
#           cat.scanone, lod.quantile.permutation,
#           make.max.N, make.maxlod, smooth.neqtl
######################################################################


#' Pull high LOD values with chr and pos.
#' 
#' Pull high LOD values with chr and pos.
#' 
#' The \code{highlod} condenses a \code{scanone} object to the peaks above a
#' \code{lod.thr} and/or within \code{drop.lod} of such peaks. The
#' \code{pull.highlod} pulls out the entries at a particular genomic location
#' or interval of locations. Summary, print, plot, max and quantile methods
#' provide ways to examine a \code{highlod} object.
#' 
#' @aliases highlod print.highlod summary.highlod plot.highlod max.highlod
#' quantile.highlod pull.highlod
#' @param scans object of class \code{\link[qtl2]{scan1}}
#' @param map map object
#' @param lod.thr LOD threshold
#' @param drop.lod LOD drop from max to keep for support intervals
#' @param extend extend support interval just past \code{drop.lod}; matches
#' behavior of \code{\link[qtl]{lodint}} when \code{TRUE}
#' @param restrict.lod restrict to loci above LOD threshold if \code{TRUE};
#' matches behavior of \code{\link[qtl]{lodint}} when \code{FALSE} (default)
#' @param chr chromosome identifier
#' @param pos position, or range of positions, in cM
#' @param x,object object of class \code{highlod}
#' @param probs probability levels for quantiles (should be > 0.5)
#' @param n.quant maximum of \code{s.quant}
#' @param n.pheno number of phenotypes considered
#' @param max.quantile return only quantiles of max LOD across genome if
#' \code{TRUE}
#' @param window size of window for smoothing hotspot size
#' @param quant.level vector of LOD levels for 1 up to
#' \code{length(quant.level)} size hotspots
#' @param sliding plot as sliding hotspot if \code{TRUE}
#' @param \dots arguments passed along
#' @return Data frame with \item{row}{row number in \code{\link[qtl]{scanone}}
#' object} \item{phenos}{phenotype column number} \item{lod}{LOD score for
#' phenotype at locus indicated by \code{row}}
#' @author Brian S Yandell and Elias Chaibub Neto
#' @seealso \code{\link{highlod}}, \code{\link{hotperm}}
#' @keywords utilities
#' @examples
#' 
#' example(include.hotspots)
#' scan1 <- scanone(cross1, pheno.col = 1:1000, method = "hk")
#' high1 <- highlod(scan1, lod.thr = 2.11, drop.lod = 1.5)
#' pull.highlod(high1, chr = 2, pos = 24)
#' summary(high1, lod.thr = 2.44)
#' max(high1, lod.thr = seq(2.11, 3.11, by = .1))
#' 
#' @export highlod
highlod <- function(scans, map, lod.thr = 0, drop.lod = 1.5,
                    extend = TRUE, restrict.lod = FALSE, ...)
{
  pheno.col <- seq(ncol(scans) - 2)

  if(is.null(lod.thr))
    lod.thr <- 0
  
  ## Extract matrix of lod scores.
  x <- scans

  ## Keep only traits with some LOD above lod threshold.
  keep <- apply(x, 2, function(x, lod.thr) any(x >= lod.thr), lod.thr)
  x <- x[, keep, drop = FALSE]
  
  ## Find which values are at within drop.lod of maximum per chr and trait.
  if(restrict.lod) {
    ## Restrict to loci above LOD threshold.
    if(extend)
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod.thr) & (maxx <= x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          (g | (c(d,0) == 1) | (c(0,d) == -1)) & (x >= lod.thr)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod.thr, drop.lod) {
        (max(x) <= x + drop.lod) & (x >= lod.thr)
      }
  }
  else { ## Do not restrict support interval to be above lod.thr
    if(extend)
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod.thr) & (maxx <= x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          g | (c(d,0) == 1) | (c(0,d) == -1)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        (maxx >= lod.thr) & (maxx <= x + drop.lod)
      }
  }
  lodint.pos <- function(x, chr, lod.thr, drop.lod) {
    unlist(tapply(x, chr, tmpfn, lod.thr, drop.lod))
  }
  scans <- tbl_map(map)
  wh <- apply(x, 2, lodint.pos, scans$chr, lod.thr, drop.lod)
  
  ## Get row and column indices.
  rr <- row(x)[wh]
  cc <- seq(keep)[keep][col(x)[wh]]

  ## Find which are within drop.lod of max lod per chr.
  lod <- x[wh]
  
  ## return data frame with genome row, trait column and lod value.
  out <- list(highlod = cbind.data.frame(row = rr, phenos = pheno.col[cc], lod = lod),
              chr.pos = scans[,1:2],
              names = names(scans)[-(1:2)])
  class(out) <- c("highlod", "list")
  attr(out, "lod.thr") <- lod.thr
  attr(out, "drop.lod") <- drop.lod
  out
}
#' @method print highlod
#' @export
#' @rdname highlod
print.highlod <- function(x, ...) print(summary(x, ...))
#' @method summary highlod
#' @export
#' @rdname highlod
summary.highlod <- function(object, ...)
{
  summary(hotsize(object, ...))
}
#' @method quantile highlod
#' @export
#' @rdname highlod
quantile.highlod <- function(x, probs = NULL, lod.thr = NULL, n.quant, n.pheno,
                             max.quantile = TRUE, ...)
{
  highlod <- highlod.thr(x, lod.thr)
  lod.thr <- highlod$lod.thr
  highlod <- highlod$highlod
  if(!nrow(highlod))
    return(NULL)
  
  ## Probabilities if requested.
  if(!is.null(probs)) {
    if(missing(n.pheno))
      stop("need to supply n.pheno along with probs")
    s.quant <- ceiling(max(probs) * n.pheno)
    n.quant <- max(s.quant)
  }
  else {
    if(missing(n.quant))
      n.quant <- max(table(highlod[,"row"]))
    s.quant <- seq(n.quant)
  }
  
  ## Quantiles from 1 to n.quant.
  if(n.quant) {
    out <- get.tails(highlod, n.quant, s.quant)
    if(max.quantile) {
      s.quant <- dimnames(out)[[1]]
      out <- apply(out, 1, max, na.rm = TRUE)
      names(out) <- s.quant
    }
    out
  }
  else
    NULL
}
#' @method max highlod
#' @export
#' @rdname highlod
max.highlod <- function(x, lod.thr = NULL, window = NULL, quant.level = NULL,
                        ...)
{
  if(is.null(window))
    window <- attr(x, "window")
  if(is.null(quant.level))
    window <- attr(x, "quant.level")
  mymax <- function(x, window = NULL, quant.level = NULL) {
    if(is.null(x)) {
      out <- data.frame(chr = NA, pos = NA, max.N = 0)
      if(!is.null(window))
        out$max.N.window <- 0
      if(!is.null(quant.level))
        out$max.lod.quant <- 0
      out
    }
    else
      max.hotsize(x)
  }
  if(length(lod.thr) > 1) {
    out <- NULL
    out.thr <- NULL
    for(lod in lod.thr) {
      tmp <- hotsize(x, lod, ...)
      out <- dplyr::bind_rows(out, mymax(tmp, window, quant.level))
      out.thr <- c(out.thr, lod)
    }
    if(!is.null(out))
      out$lod.thr <- out.thr
    
    class(out) <- c("summary.scanone", "data.frame")
    out
  }
  else
    mymax(hotsize(x, lod.thr, window, quant.level, ...), window, quant.level)
}
#' @method plot highlod
#' @export
#' @rdname highlod
plot.highlod <- function(x, ..., quant.level = NULL, sliding = FALSE)
{
  if(sliding) {
    if(is.list(quant.level))
      quant.level <- quant.level$max.lod.quant
    
    ## Need to supply quant.level as second argument.
    slidingbar.plot(slidingbar.create(x, quant.level, ...), ...)
  }
  else
    graphics::plot(hotsize(x, ..., quant.level = quant.level), ...)
}
###################################################################################
highlod.thr <- function(highobj, lod.thr)
{
  if(is.null(lod.thr))
      lod.thr <- attr(highobj, "lod.thr")

  if(!is.null(lod.thr)) {
    highobj$highlod <- highobj$highlod[highobj$highlod$lod >= lod.thr,, drop = FALSE]
    attr(highobj, "lod.thr") <- lod.thr
  }
  
  highobj
}
###################################################################################
pull.highlod <- function(object, chr, pos, ...)
{
  ## Kludge to get names if not in object
  pheno.names <- object$names
  if(is.null(pheno.names)) {
    extra <- list(...)
    m <- match("names", names(extra))
    if(!is.na(m))
      pheno.names <- extra[[m]]
  }
  wh.chr <- which(object$chr.pos$chr == chr)
  ## Want to expand this to handle range of positions...
  wh.pos <- which(object$chr.pos$pos[wh.chr] - min(pos) >= 0 &
                  object$chr.pos$pos[wh.chr] - max(pos) <= 0)
  wh.high <- which(object$highlod$row %in% wh.chr[wh.pos])
  wh.row <- object$highlod[wh.high, "row"]

  out <- data.frame(object$chr.pos[wh.row,], object$highlod[wh.high, -1])
  ## Add phenotype names if available.
  if(!is.null(pheno.names))
    out$phenos <- pheno.names[out$phenos]
  out
}
