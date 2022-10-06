library(dplyr)
library(data.table)
pangolin <- snakemake@input[["pangolin"]]
nextclade <- snakemake@input[["nextclade"]]
output <- snakemake@output[["lineage_ids"]]

lg <- snakemake@params[["lineage"]]
print(paste("looking for", lg, "lineage"))
nextcladedf <- fread(nextclade,select=c("seqName","totalMissing"),showProgress=TRUE,colClasses=c(seqName="character",totalMissing="integer")) 
pangolindf <- fread(pangolin,select=c("taxon","lineage"),showProgress=TRUE,colClasses=c(taxon="character",lineage="character")) 

nextcladedf <- nextcladedf %>%
    filter(totalMissing <= 1000) %>%
    left_join(pangolindf, by=c("seqName"="taxon")) %>%
    filter(lineage==lg) %>%
    select(seqName)
fwrite(nextcladedf,output,showProgress=TRUE,col.names=FALSE)

