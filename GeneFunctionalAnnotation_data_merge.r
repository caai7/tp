install.packages("readr") # 下載 readr 套件
install.packages("dplyr") # 下載 dplyr 套件
install.packages("tidyr") # 下載 tidyr 套件
library(readr)
library(dplyr)
library(tidyr)
setwd("D:/genome_assembly/no_filter/annotation/")
#D:\download\anno
interpro <- read.delim("tala2_PNC_iprs.interpro.anno.txt",header=T,sep='\t')
panther <- read.delim("tala2_PNC_iprs.panther.anno.txt",header=T,sep='\t')
pfam <- read.delim("tala2_PNC_iprs.pfam.anno.txt",header=T,sep='\t')
pannzer.go <- read.delim("tala2_PNC_pannzer2.go.anno.txt",header=T,sep='\t')
pannzer.uniport <- read.delim("tala2_PNC_pannzer2.uniprot.anno.txt",header=T,sep='\t')
nr <- read.delim("tala2.nr.protein_annoation.txt",header=T,sep='\t')
ko <- read.delim("tala2_query.ko.txt",header=T,sep='\t') 

merged_data <- full_join(ko, interpro, by = "GeneID")
merged_data <- full_join(merged_data, panther, by = "GeneID")
merged_data <- full_join(merged_data, pfam, by = "GeneID")
merged_data <- full_join(merged_data, pannzer.go, by = "GeneID")
merged_data <- full_join(merged_data, pannzer.uniport, by = "GeneID")
merged_data <- full_join(merged_data, nr, by = "GeneID")

write.csv(merged_data, "tala_merged_functional_annotations.csv", row.names = FALSE) 

####Done######



#(optional)把GO 重複的地方整合
#因本次無採用eggNOG-mapper2的data故不採用此步驟
consolidate_go_terms <- function(row) {
  go_terms <- unique(c(
    if (!is.na(row$pannzer2.go)) unlist(strsplit(row$pannzer2.go, ",")) else NULL,
    if (!is.na(row$eggnog.go)) unlist(strsplit(row$eggnog.go, ",")) else NULL
  ))
  return(paste(go_terms, collapse = ","))
}

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(Consolidated_GO = consolidate_go_terms(cur_data_all()))


