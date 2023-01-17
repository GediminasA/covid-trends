library(dplyr)
library(data.table)
inf = snakemake@input[[1]]
targetamount = snakemake@params[["target_amount"]]
out_all_ids = snakemake@output[["ids"]]
out_target_ids = snakemake@output[["ids_target"]]
out_addiotan_ids = snakemake@output[["ids_additional"]]  
additional_id_cut_off = 0.996
additional_id_amount = 1000
df <- fread(inf,sep="\t",header = FALSE,fill=T)
df2 <- df %>%
  select(Sequence=V1,Target=V6,Length=V2,LengthRef=V7,Start=V3,End=V4,Matches=V10,AlignmentLength=V11)%>%
  group_by(Sequence) %>%
  arrange(-Matches, .by_group = TRUE) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  mutate(Identity4Seq=Matches/abs(End-Start),Coverage=abs(End-Start)/Length,CoverageRef=abs(End-Start)/30000,Identity4Aln=Matches/AlignmentLength,Identity4Ref=Matches/LengthRef) %>%
  mutate(Identity4RefRound=round(Identity4Ref,digits=5)) %>%
  mutate(Identity4SeqRound=round(Identity4Seq,digits=5)) %>%
  filter(Identity4RefRound >= 0.9000)

dfadditional <- df2 %>%
  filter(Identity4SeqRound >= additional_id_cut_off)
additional_ids = as.vector(unique(dfadditional$Sequence))

range_c = seq(0.9997,1.0,0.00001)
range_n  = lapply(range_c,function(s){
  df = df2 %>%
    filter(Identity4SeqRound >= as.double(s))
  return(nrow(df))
})

dfw = data.frame(limit=range_c,count = unlist(range_n))
print(dfw)
dfw2 <- dfw %>%
  filter(count <= targetamount) %>%
  arrange(-count)
limit_cnt = dfw2$limit[1]
df2 <- df2 %>%
  mutate(Identity4SeqRound = as.double(Identity4SeqRound)) %>%
  filter(Identity4SeqRound >= as.double(limit_cnt))
print(paste("Chosen",nrow(df2),"sequences with identity cut-off ",limit_cnt))


target_ids =  as.vector(unique(df2$Sequence))
additional_ids <- setdiff(additional_ids,target_ids)
if (length(additional_ids) >= additional_id_amount) {
  additional_ids <- sample(additional_ids,additional_id_amount,replace = F)
  print(paste("Sampled ",additional_id_amount,"of additional sequences"))
  } else {
    lad = length(additional_ids)
    print(paste("Aded ",lad,"Additional sequences"))
  }

out_additional_ds = data.table(IDS = additional_ids)
out_target_ds = data.table(IDS = target_ids)
out_all_ds = bind_rows(out_target_ds,out_additional_ds)

fwrite(out_additional_ds,out_addiotan_ids,col.names = F)
fwrite(out_target_ds,out_target_ids,col.names = F)
fwrite(out_all_ds,out_all_ids,col.names = F)

