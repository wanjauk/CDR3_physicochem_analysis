library(openxlsx)
library(tidyverse)

rnaseq_input <- read.xlsx("RNAseq based CDR3s for aromaticity and UH.xlsx")

cancers <- c("LUSC", "KIRC", "CESC",
             "LUAD", "PAAD", "COAD",
             "BLCA", "BRCA", "PRAD",
             "SKCM")

for (cancer in cancers){
  
  cancer_cdr3_properties <- rnaseq_input[rnaseq_input$Study==cancer,]
  
  # rename the columns
  cancer_cdr3_properties <- cancer_cdr3_properties %>% rename_at("SampleTypeLetterCode", ~"Sample.Type")
  cancer_cdr3_properties <- cancer_cdr3_properties %>% rename_at("ParticipantBarcode", ~"Filename")
  
  #exclude NAs
  cancer_cdr3_properties <- cancer_cdr3_properties[complete.cases(cancer_cdr3_properties),]
  
  # write to excel file
  write.xlsx(cancer_cdr3_properties, file = paste0("VDJ_Recoveries Boris ", cancer,".xlsx"), sheetName = "physicochem")
  
  
}


