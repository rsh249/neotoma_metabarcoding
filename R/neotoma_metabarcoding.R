#### R Script to analyze DNA metabarcode data obtained from ancient packrat middens
#load libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(rentrez)
pal = brewer.pal(name='Paired', n = 12)

# get data from zenodo repository and 
options(timeout=9999)
download.file("https://zenodo.org/record/5110824/files/SURE2021.tar.gz?download=1", 'SURE2021.tar.gz')
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
krep.files = list.files('results', full.names = T, recursive = T, pattern='kreport')
l <- lapply(krep.files, data.table::fread)
for(i in 1:length(l)){
  filename = krep.files[[i]]
  #samname = tools::file_path_sans_ext(basename(filename)) 
  samname = paste(word(filename, start=-2, sep='/'), word(filename, start=-1, sep='/'), sep='_')
  samname = word(samname, start = 1, sep = '[.]')
  l[[i]]$samname = samname
}

krep <- bind_rows(l)

kproc.sample = krep %>%
  filter(V4=='F') %>% 
  filter(V1>=0.05) %>%
  group_by(samname, V6) %>%
  summarise(readcount = sum(V2)) %>%
  mutate(freq = readcount / sum(readcount)) %>% 
  filter(freq>0.05)



ktileplot_scaled = ggplot(data = kproc.sample ) +
  geom_tile(aes(x=samname, y = V6, fill = readcount)) +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.y = element_text(angle=45),
        axis.text.x = element_text(angle=45, vjust =1, hjust = 1)) +
  scale_fill_steps(low=pal[1], high=pal[2], n.breaks = 10) +
  xlab('') + ylab('')
ktileplot_scaled

ktileplot = ggplot(data = krep %>% filter(V4=='G') %>% filter(V2>5) ) +
  geom_tile(aes(x=samname, y = V6, fill = V2)) +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.y = element_text(angle=45),
        axis.text.x = element_text(angle=45, vjust =1, hjust = 1)) +
  scale_fill_steps(low=pal[1], high=pal[2], n.breaks = 10) +
  xlab('') + ylab('')
ktileplot

