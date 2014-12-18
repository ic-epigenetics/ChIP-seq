require(ShortRead)

removeSRduplicates <- function(inputFile){
  RA <- readAligned(inputFile,type="BAM",filter=nFilter())
  RA <- RA[!is.na(position(RA))]
  RA[!duplicated(position(RA))]
}
