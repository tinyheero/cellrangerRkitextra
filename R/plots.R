source("scripts/utils.R")
library(tidyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(ggfortify)
library(cowplot)
library(purrr)
library(doMC)
library(foreach)

library("dplyr")
library("readr")
library("cellrangerRkit")

theme_set(theme_grey())

#Single Cell Data
matrix.df <- read_delim("10X_demo/sample_results/293t/filtered_matrices_mex/hg19/matrix.tsv", 
                        col_names = FALSE, delim = " ")

matrix.df <- rename(matrix.df, gene_ind = X1, barcode_ind = X2)

#----------
# Loading Gene List Data
#----------
gene.list.df <- read_tsv("10X_demo/sample_results/293t/filtered_matrices_mex/hg19/genes.tsv", 
                         col_names = c("gene_id", "gene_name"))

gene.list.df <- 
  gene.list.df %>%
  mutate(gene_ind = row_number())

#----------
# Loading Barcode Data
#----------
barcode.df = read_tsv("10X_demo/sample_results/293t/filtered_matrices_mex/hg19/barcodes.tsv",
                      col_names = c("barcode"))

barcode.df <- 
  barcode.df %>%
  mutate(barcode_ind = row_number())

matrix.gene.barcode.df <- 
  matrix.df %>%
  left_join(gene.list.df) %>%
  left_join(barcode.df)


full_gene_list = read.table("10X_demo/sample_results/293t/matrices_mex/hg19/genes.tsv",header=F,sep="\t")
full_barcode_list = read.table("10X_demo/sample_results/293t/matrices_mex/hg19/barcodes.tsv",header=F,sep="\t")

93t.df

#293t.df = 293t.df[-1,]
#
#293t.df$ens_id = gene_list[as.numeric(as.character(293t.df$V1)),1]
#293t.df$barcode = barcode_list[as.numeric(as.character(293t.df$V2)),1]
293t.df = 293t.df[,c(4,5,3)]

df_293_spread = spread(293t.df,"barcode","V3",fill=NA)
df_293_spread = firstColAsRowNames(df_293_spread)

exp_gene_count = apply(df_293_spread,1,function(x){length(x[!is.na(x)])})
qplot(exp_gene_count,
      geom="histogram", 
      main = "Histogram For Number of Expressed Genes", 
      xlab = "Count")


#Calculating gene coverage saturation plot
registerDoMC(30)
gene_counts_df = foreach(i = c(1:ncol(df_293_spread)))%dopar%{
  x = as.numeric(df_293_spread[,i])
  x[is.na(x)] = 0
  x = sample(x,length(x))
  csum = unlist(lapply(c(1:length(x)),function(y){
    sum(x[1:y])
  }))
  csum
}
gene_counts_df = do.call(rbind,gene_counts_df)
gene_counts_df = t(gene_counts_df)
rownames(gene_counts_df) = c(1:nrow(gene_counts_df))
colnames(gene_counts_df) = colnames(df_293_spread)

count_int = seq(0,52302,1000)
gene_counts_df_summary = foreach(i = c(1:ncol(gene_counts_df)))%dopar%{
  x = as.numeric(gene_counts_df[,i])
  gene_read_counts = df_293_spread[,i]
  unlist(lapply(count_int,function(z){
    if(max(x)>=z){
      y = abs(x-z)
      ywhich = which(y==min(y))[1]
      num_genes = gene_read_counts[1:ywhich]
      num_genes = num_genes[!is.na(num_genes)]
      rr = length(num_genes)
    }else{
      rr = NA
    }
    rr
  }))
}
gene_counts_df_summary = do.call(rbind,gene_counts_df_summary)
rownames(gene_counts_df_summary) = colnames(gene_counts_df)
colnames(gene_counts_df_summary) = count_int

save(gene_counts_df_summary,file="../10X_demo/sample_results/293t/filtered_matrices_mex/hg19/gene_counts_df_summary.RData")
head(gene_counts_df_summary)

gene_counts_df_summary.melt <- 
  melt(gene_counts_df_summary, value.name = "num_genes") %>%
  rename(barcode = Var1, num_mapped_reads = Var2) %>%
  arrange(num_mapped_reads) %>%
  mutate(num_mapped_reads = factor(as.character(num_mapped_reads), 
                                   levels = unique(num_mapped_reads)))

gene_counts_df_summary.melt.num.na <- 
  gene_counts_df_summary.melt %>%
  group_by(num_mapped_reads) %>%
  filter(!is.na(num_genes)) %>%
  summarize(num_non_na = n())

gene_counts_df_summary.melt %>%
  ggplot(aes(x = barcode, y = num_genes)) +
  geom_boxplot()

  group_by(barcode) %>%
  summarize(mean_num_genes = mean(num_genes, na.rm = TRUE)) %>%
  ))

p1 <- 
  gene_counts_df_summary.melt %>%
  ggplot(aes(x = num_mapped_reads, y = num_genes)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Number of Genes Represented") + 
  xlab("Number of Mapped Reads In Genic Regions")

p2 <- 
  gene_counts_df_summary.melt.num.na %>%
  ggplot(aes(x = num_mapped_reads, y = num_non_na)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_text(aes(label = num_non_na), vjust = -0.5) +
  ylab("Number of Cells Represented By The Box") + 
  xlab("Number of Mapped Reads In Genic Regions")

coverage_saturation_plot = plot_grid(p1, p2, nrow = 2)
coverage_saturation_plot
save.image("coverage_saturation.RData")
