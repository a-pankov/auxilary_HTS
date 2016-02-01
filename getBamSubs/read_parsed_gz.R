read_parsed_gz <- function(file_location){
  tst <- readLines(gzfile(file_location))
  split_df <- data.frame(matrix(unlist(strsplit(tst[-(1:4)], '\t')), ncol = 11, byrow = T), stringsAsFactors = F)
  
  split_df[,1] <- factor(split_df[,1], levels = c("0", "1"))  
  split_df[,2] <- factor(split_df[,2], levels = c("0", "1"))
  split_df[,4] <- as.integer(split_df[,4])
  split_df[,6] <- as.integer(split_df[,6])
  split_df[,5] <- factor(split_df[,5])
  split_df[,8] <- factor(split_df[,8])
  split_df[,9] <- factor(split_df[,9])
  split_df[,10] <- factor(split_df[,10])
  
  names(split_df) <- c("isFirstMate", "isForwardStrand", "cigar", "readLength", "chrom", "md", "queryBases", "substitution", "queryQualities", "readInd")
  
  return(split_df)
}

