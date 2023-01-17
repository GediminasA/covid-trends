library(dplyr)
library(data.table)
library(lubridate)
pangolin <- snakemake@input[["pangolin"]]
nextclade <- snakemake@input[["nextclade"]]
lineages <- snakemake@params[["lineages"]] 
meta <- snakemake@input[["meta"]]
break_date <- snakemake@params[["break_date"]]
start_date <- snakemake@params[["start_date"]]
out_chosenini_ids <- snakemake@output[["chosen_ids"]]
out_all_ids <- snakemake@output[["all_ids"]]
out_chosenini_data <- snakemake@output[["chosen_data"]]
out_all_data <- snakemake@output[["all_data"]]


lg <- strsplit(lineages,split=" ")[[1]]

print(paste("looking for", lg, "lineage"))
nextcladedf <- fread(nextclade,select=c("seqName","totalMissing"),showProgress=TRUE,colClasses=c(seqName="character",totalMissing="integer")) 
pangolindf <- fread(pangolin,select=c("taxon","lineage"),showProgress=TRUE,colClasses=c(taxon="character",lineage="character")) 
metadf <- fread(meta,select=c("strain","date","country"),showProgress=TRUE,colClasses=c(strain="character",country="character",date="character")) 

nextcladedf <- nextcladedf %>%
  filter(totalMissing <= 1000) %>%
  left_join(pangolindf, by=c("seqName"="taxon")) %>%
  #filter(lineage %in% lg) %>%
  select(seqName,lineage)
nextcladedf2 <- nextcladedf %>%
  left_join(metadf, by=c("seqName"="strain")) %>%
  filter(country=="Lithuania" & nchar(date)==10)
nextcladedf2$dataparsed <-lubridate::as_date(nextcladedf2$date)

alldata <- nextcladedf %>%
  left_join(metadf, by=c("seqName"="strain")) %>%
  filter( nchar(date)==10) %>%
  mutate(dataparsed=lubridate::as_date(date)) %>%
  filter(!country=="Lithuania")


part1 <- nextcladedf2 %>%
  filter(dataparsed < lubridate::as_date(break_date) &  dataparsed >= lubridate::as_date(start_date))
npart1 <- nrow(part1)
part2ini <- nextcladedf2 %>%
  arrange(dataparsed) %>%
  filter(dataparsed > lubridate::as_date(break_date) ) %>%
  mutate(N = 1) %>%
  mutate(cum_count = cumsum(N)) %>%
  filter(cum_count <= npart1)
end_break <- max(part2ini$dataparsed)
start_break <- lubridate::as_date(start_date)

alldata <- alldata %>%
  filter(dataparsed >= start_break & dataparsed <=  end_break )

outdata <- bind_rows(part1,part2ini)
outdata_ids <- outdata %>%
  select(seqName)
alldata_ids <- alldata %>%
  select(seqName)
fwrite(outdata,out_chosenini_data)
fwrite(outdata_ids,out_chosenini_ids,col.names = F)
fwrite(alldata,out_all_data)
fwrite(alldata_ids,out_all_ids,col.names = F)