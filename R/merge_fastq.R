## Script to concatenate fastq_pass data by barcode for each sample.

# libraries
library(Biostrings)
library(stringr)

options(timeout=9999)

pool = list.files('./SURE2021', full.names = T)
for(p in 1:length(pool)){
  target = pool[p]
  target_str = list.files(target, recursive=T, full.names = T)
  fastq_pass = target_str[grepl('barcode', target_str) & grepl('fastq_pass', target_str)]
  
  outname = word(target, start=-1, sep = "/")
  
  # separate fastq_pass by barcode and write new single fastq files for each folder
  barcode_id = word(fastq_pass, start=-2, sep="/")
  dir.create('data', showWarnings = F)
  dir.create(paste('data/', outname, sep=''), showWarnings = F)
  for(i in unique(barcode_id)){
    filepath = paste('data/', outname, "/", i, ".fastq", sep='')
    agg_read = readDNAStringSet(fastq_pass[grepl(i, fastq_pass)], format = 'fastq')
    # filter longer reads out to prevent write errors
    agg_read = agg_read[which(agg_read@ranges@width<10000)]
    
    ## Write an XStringSet object to a FASTA (or FASTQ) file:
    writeXStringSet(agg_read, 
                    filepath, 
                    append=FALSE,
                    compress=FALSE,  
                    format="fastq")
    
  }
}
