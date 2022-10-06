library(dplyr)
nextclade <- read.csv(snakemake@input[[1]],sep = '\t')
col.namesf2 = strsplit("query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",split=',')[[1]] 
mmseq <- read.csv(snakemake@input[[2]],sep = '\t',
                      row.names=NULL,
                      header = FALSE,
                      col.names = col.namesf2)
mergd <- left_join(nextclade,mmseq,by=c("seqName"="target"))
mergd <- mergd %>%
  filter(totalMissing <= 100) %>% 
  filter(fident > 0.95) %>% 
  filter( 
          grepl("S:E156",aaDeletions) | 
          grepl( "S:F157",aaDeletions) | 
          grepl( "S:R158",aaDeletions) | 
           grepl("S:S494P",aaSubstitutions) | 
           grepl("S:E484K",aaSubstitutions) 
          )
out <- mergd %>%
  select(seqName)
write.table(out,snakemake@output[[1]],sep="", row.names = FALSE, quote=FALSE, col.names = FALSE)
