#### R Script to analyze DNA metabarcode data obtained from ancient packrat middens
#load libraries
library(stringr)

# get data from zenodo repository and unzip
download.file("https://zenodo.org/record/5068129/files/SURE2021.tar.gz?download=1", 'SURE2021.tar.gz')
system('tar -xvf SURE2021.tar.gz')


# move data from pool to data folders
source('./R/merge_fastq.R')

# setup kraken databases
#source('./R/get_amplicon_db.R') # this takes a long time, run only once!

# run kraken search
data = list.files('data', full.names = T)
dir.create('results', showWarnings = F)
for(d in data){
  barcodes = list.files(d, full.names = T)
  sambase = word(d, start=-1, sep='/')
  results = paste('results/', sambase, sep='')
  dir.create(results, showWarnings = F)
  print(d)
  print(results)
  
  for(b in barcodes){
    # run kraken query
    print(b)
    base = word(b, start=-1, sep ='/')
    base = word(base, start=1, sep='[.]')
    kdb = 'kraken_db'
    nclus = 4
    kreport = paste(results, "/", base, '.kreport', sep = '')
    readsout = paste(results, "/", base, '.fastq', sep = '')
    kout = paste(results, "/", base, '.kout', sep = '')
    run_kraken = paste('kraken2 --db',
                       kdb,
                       '--threads',
                       nclus, 
                       '--use-names --report', 
                       kreport, 
                       '--classified-out',
                       readsout,
                       b,
                       ">", 
                       kout)
    system(run_kraken)
  }
}

# visualize results


