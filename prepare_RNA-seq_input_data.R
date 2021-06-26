library(openxlsx)
library(tidyverse)

rnaseq_input <- read.xlsx("RNA-seq_data/RNAseq based CDR3s for aromaticity and UH with alpha beta TCR indicated.xlsx")

cancers <- c("LUSC", "KIRC", "CESC",
             "LUAD", "PAAD", "COAD",
             "BLCA", "BRCA", "PRAD",
             "SKCM")

for (cancer in cancers){
  
  cancer_cdr3_properties <- rnaseq_input[rnaseq_input$Study==cancer,]
  
  # rename the columns
  cancer_cdr3_properties <- cancer_cdr3_properties %>% rename_at("SampleTypeLetterCode", ~"Sample.Type")
  cancer_cdr3_properties <- cancer_cdr3_properties %>% rename_at("ParticipantBarcode", ~"Filename")
  cancer_cdr3_properties <- cancer_cdr3_properties %>% mutate(Sample.Type = ifelse(Sample.Type == "TP", "Primary Tumor", Sample.Type),
                                                              Sample.Type = ifelse(Sample.Type == "TM", "Metastatic", Sample.Type))
  cancer_cdr3_properties <- cancer_cdr3_properties %>% mutate(Receptor = ifelse(chain == "alpha", "TRA", "TRB"))
  
  #exclude NAs
  cancer_cdr3_properties <- cancer_cdr3_properties[complete.cases(cancer_cdr3_properties),]
  
  # write to excel file
  write.xlsx(cancer_cdr3_properties, file = paste0("RNA-seq_data/VDJ_Recoveries Boris ", cancer,".xlsx"), sheetName = "physicochem")
  
  
}


