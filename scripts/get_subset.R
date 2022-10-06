library(dplyr)
library(data.table)
allids <- snakemake@input[["lineage_ids"]]
meta <- snakemake@input[["meta"]]
output <- snakemake@output[["lineage_seqs_lt"]]
maxout <- snakemake@params[["subset_max"]]

allidsdf <- fread(allids,col.names=c("seqName"),showProgress=TRUE,colClasses=c(seqName="character"))
metadf <- fread(meta,select=c("strain","date","region","country","division"),showProgress=TRUE)



allidsdf <- allidsdf %>%
    left_join(metadf,by=c("seqName"="strain")) %>%
    filter(country=="Lithuania") %>%
    select(seqName,date) %>%
    slice_min(order_by=date,n = maxout,with_ties = T) %>%
    select(seqName)


fwrite(allidsdf,output,showProgress=TRUE,col.names=FALSE)