library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(tidyr)

meta <- fread(snakemake@input$data)
head(meta)
nrow(meta)

ids <- fread(snakemake@input$ids, header = F)
ids <- unlist(ids$V1)
meta <- meta %>%
    filter(taxon %in% ids) %>%
    as.data.frame()
nrow(meta)

lineage4focus=snakemake@params$id

df_focus = meta %>%
    lazy_dt() %>%
    mutate(date = date(date)) %>%
    mutate(week = lubridate::floor_date(date, "week")) %>%
    as.data.table() 
head(df_focus)


df_focus_per_week <- df_focus %>%
    lazy_dt() %>%
    group_by(week) %>%
    summarise(N=n(), IDs=list(taxon)) %>%
    as.data.table()

set.seed(snakemake@wildcards$seed)
downsampled = list()
labels = list()
ct = 0
for (lini in df_focus_per_week$IDs) {
    ct = ct + 1
    n_ini = length(lini)
    n_after = round(n_ini*as.double(snakemake@wildcards$perc)/100,digits = 0)
    if (n_after == 0) n_after = 1 
    downsampled[[ct]]=sample(lini,n_after,replace=FALSE)
    labels[[ct]] = paste("week",ct,sep="")
}

df_focus_per_week$ID_sampled <- downsampled
df_focus_per_week$Label <- labels

#output info
outdf <- df_focus_per_week %>%
    select(Label,ID_sampled) %>%
    as.data.frame() %>%
    unnest(ID_sampled)
fwrite(x = outdf, file = snakemake@output$partitions_info)

outdf2 <- outdf %>%
    select(ID_sampled) %>%
    as.data.frame()
fwrite(x = outdf2, file = snakemake@output$cluster, row.names = F, col.names = F)
