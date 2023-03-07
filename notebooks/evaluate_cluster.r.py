#!/usr/bin/env python
# coding: utf-8

# In[2]:


library(dplyr)
library(data.table)
library(tidyr)
library(dtplyr)
library(stringr)


# In[3]:


cldf <- fread(snakemake@input$merged)


# In[5]:


head(cldf)


# In[6]:


df <- fread(snakemake@input$cluster)


# In[8]:


data <- fread(snakemake@input$data, select = c("taxon","month"))


# In[9]:


cldf2 <- cldf %>%
    left_join(data, by = c("Member"="taxon")) %>%
    mutate(Cluster = str_replace_all(Cluster,pattern = "'",replacement = "")) %>%
    left_join(df, by=c("Cluster"="SequenceName")) %>%
    as.data.frame()


# In[10]:


head(cldf2)


# In[11]:


mnths = unique(sort(data$month))
m1 = mnths[[1]]
m2 = mnths[[length(mnths)]]


# In[15]:


df2 <- cldf2 %>%
    filter(month == m1 | month ==m2)
df3 <- df2  %>%
    group_by(ClusterNumber, month) %>%
    summarise(N=n()) %>%
    tidyr::pivot_wider(names_from = month, values_from = N) %>%
    na.omit() %>%
    rowwise() %>%
    mutate(total_nmb = sum(c_across(starts_with("20")))) %>%
    as.data.frame()
df3$ratio <- df3[,2] / df3[,3]
good <- df3 %>%
    filter(total_nmb >= 200) %>%
    filter(ratio < 1, ratio > 0.1, ClusterNumber != -1)
# good <- df3 %>%
#     filter(total_nmb >= 200) %>%
#     filter(ratio < 1, ratio > 0.1, ClusterNumber != -1)
sd = sd(good$total_nmb)

out=data.frame(File = snakemake@input$cluster ,N=nrow(good), Nall=nrow(df3), snakemake@wildcards$p,
              meanv = mean(good$total_nmb), minv = min(good$total_nmb), maxv = max(good$total_nmb), sdv = sd(good$total_nmb))
fwrite(x = out, row.names = F,col.names = F, file = snakemake@output[[1]],na='NA')

