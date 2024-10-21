cat.scanone <- function(dirpath = ".", filenames = permfiles, chr.pos)
{
  ## Folder should contain qtl::scanone highlods data across all traits for ONE permutation
  permfiles <- list.files(dirpath, paste("per.scan", "*", "RData", sep = "."))
  
  ## Make and remove per.scan.hl. Below use version from files.
  per.scan.hl <- NULL
  rm(per.scan.hl)
  
  for(i in 1:length(filenames)) {
    highobj <- with(filenames[i], {
      if(i==1) per.scan.hl else rbind.data.frame(highobj, per.scan.hl)
    })
  }
  cbind.data.frame(chr.pos[highobj$row,],highobj)
}
