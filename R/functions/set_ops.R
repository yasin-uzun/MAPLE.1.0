library("GenomicRanges")

#a=bed_gr_region
#b=bed_gr_cov

intersect_bed <- function(a, b){
  library(GenomicRanges)
  my_hit <- findOverlaps(a, b)
  
  da = as.data.frame(a[queryHits(my_hit)])
  db = as.data.frame(b[subjectHits(my_hit)])
  
  colnames(db) = paste0('b_', colnames(db))
  
  inter  <- cbind(da, db)
  return(inter)
}


union_bed <- function(a, b){
  library(GenomicRanges)
  my_hit <- findOverlaps(a, b)
  
  da = as.data.frame(a[queryHits(my_hit)])
  db = as.data.frame(b[subjectHits(my_hit)])
  
  colnames(db) = paste0('b_', colnames(db))
  
  inter  <- cbind(da, db)
  return(inter)
}
