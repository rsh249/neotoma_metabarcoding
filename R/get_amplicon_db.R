# script to get gene specific databases:
require(rentrez)
set_entrez_key("55549707f8ed138bcb423d98d739c7789408")
Sys.getenv("ENTREZ_KEY")


get_amplicon_db <- function(target, db.dir = 'kraken_db', key = NULL, threads=2) {
  #get fasta file with gene sequences from nuccore
  #write to blast database
  #works for plants only
  #no fuzzy matching of gene names? No handling of gene IDs. You get what you get.
  
  require(rentrez)

  by = 1000
  query = paste(target, '[All Fields] AND (plants[filter] AND ("0"[SLEN] : "100000"[SLEN]))',sep='') # allow more flexibility in searches
  genesearch <-
    entrez_search(db = "nuccore",
                  term = query,
                  use_history = TRUE, key=key)
  db_file = paste(db.dir, '/', target, ".fasta", sep='')
  if(dir.exists(db.dir)){} else {dir.create(db.dir)}
  get_gene = function(x) {
    attempt <- 1
    recs = NULL
    while( is.null(recs) && attempt <= 10 ) {
      attempt <- attempt + 1
      recs <-
        try(
            entrez_fetch(
            db = "nuccore",
            web_history = genesearch$web_history,
            rettype = "fasta",
            retmax = by,
            retstart = x, 
            key=key
          ),
          silent=TRUE
        )
    } 
    
    cat(recs, file = db_file, append = TRUE)
    cat(x + by, "sequences downloaded out of", genesearch$count, "\r")
  }
  sapp = lapply(seq(1, genesearch$count, by), get_gene)
  
  if(dir.exists(paste(db.dir, "/taxonomy"))){
    
  } else{
    getKDBtaxonomy = paste('kraken2-build --download-taxonomy --db ', db.dir)
    system(getKDBtaxonomy) # get taxonomy files
  }  
  
  fasta = list.files(db.dir, pattern='fasta$', full.names = T)
  print(fasta)
  
  for(f in fasta) {
    kdb_add = paste('kraken2-build --add-to-library', f, '--db', db.dir)
    system(kdb_add)
  }
  kdb_build = paste('kraken2-build --build --threads', threads, '--db', db.dir)
  system(kdb_build)
 
  return()
}

options(timeout=999999)

get_amplicon_db('rbcL', key =Sys.getenv("ENTREZ_KEY"))
get_amplicon_db('trnL', key =Sys.getenv("ENTREZ_KEY"))
get_amplicon_db('psbA', key =Sys.getenv("ENTREZ_KEY"))
get_amplicon_db('MATK', key =Sys.getenv("ENTREZ_KEY"))
get_amplicon_db('ITS2', key =Sys.getenv("ENTREZ_KEY"))
get_amplicon_db('18S', key =Sys.getenv("ENTREZ_KEY"))



