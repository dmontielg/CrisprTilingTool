
setwd("~/Documents/github/CrisprTilingTool")
library(dplyr)
source("get-phase-trancript.R")
source("get-sgRNA.R")
source("locate-base-changes.R")
source("run_locate_modifications.R")

CDS_Human   <- readr::read_csv("CDS-human.csv")
Genome      <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
#gene.search <- 'TP53'
trans <- 'ENST00000269305.9'
#gene.search <- 'CFTR'
#trans <- 'ENST00000003084.11'
flank <- 20
sgrna_len <- 20
window_start <- 4
window_end <- 8
#base_change <- c("C->T")
base_change <- c("A->G")
# Append vector on sgRNA
seq_before <- "GAGCCTCGTCTCCCACCG"
seq_after <- "GTTTTGAGACGCATGCTGCA"

#------------------------------------------------------------------------------#

#df.target.trans <- CDS_Human %>% dplyr::filter(gene == gene.search, tx == trans)
df.target.trans <- CDS_Human %>% dplyr::filter(tx == trans)

gene.search <- unique(df.target.trans$gene)

df_table_out <- run_locate_modifications(df.target.trans,
                                         gene.search,
                                         trans,
                                         flank,
                                         sgrna_len,
                                         window_start,
                                         window_end,
                                         base_change,
                                         seq_before,
                                         seq_after, 
                                         Genome)

gc()

#print(filename_out)
#write.table(x = df_table_all, file = filename_out, sep = "\t", quote = F)