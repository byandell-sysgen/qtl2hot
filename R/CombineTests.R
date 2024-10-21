CombineTests <- function(comap, file)
{
  reg.nms <- names(comap)
  out <- NULL
  join.out <- list()
  for (k in 1 : length(comap)) {
    load(paste(file, reg.nms[k], "Rdata", sep="."))
    join.out[[k]] <- out
  }
  names(join.out) <- reg.nms
  join.out
}
