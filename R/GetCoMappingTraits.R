GetCoMappingTraits <- function(highobj, cand.reg)
{
  chr.pos <- highobj$chr.pos
  chrs <- unique(chr.pos$chr)
  chr <- ordered(cand.reg$peak.chr, chrs)
  phys <- cbind(phenos = match(as.character(cand.reg[,1]), highobj$names),
                chr = unclass(chr), pos = cand.reg$peak.pos)
  
  in.range <- function(pos, x.pos) {
    r <- range(pos)
    r[1] <= x.pos & r[2] >= x.pos
  }
  ## Find traits that are co-mapping.
  find.comap <- function(x, highlod, chr.pos, chrs, all.traits) {
    ## Traits must map to same chromosome.
    h <- highlod[chrs[x[2]] == chr.pos$chr[highlod$row],, drop = FALSE]
    ## And peaks must be in range.
    h <- tapply(chr.pos$pos[h$row], h$phenos, in.range, x[3])
    ## But have to remove trait from its list.
    h <- as.numeric(names(h[h]))
    h <- h[-match(x[1], h)]
    all.traits[h]
  }
  
  ## This list is too restrictive compared with earlier list of Elias.
  ## Try re-running deprecated code using scan.orf to compare.
  
  out <- apply(phys, 1, find.comap, highobj$highlod, chr.pos, chrs, highobj$names)
  if(is.matrix(out))
    out <- as.data.frame(out)
  names(out) <- as.character(cand.reg[,1])
  
  lapply(out, as.character)
}
