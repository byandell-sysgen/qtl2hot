tbl_map <- function(map) {
  tibble::as_tibble(
    data.frame(
      chr = factor(rep(names(map), sapply(map, length)), names(map)),
      pos = as.vector(unlist(map)),
      mar = c(sapply(map, names))))
}