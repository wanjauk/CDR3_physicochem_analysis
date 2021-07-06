library(openxlsx)
library(tidyverse)


for (cancer in cancers) {
 
  exp_output_dir <- paste("Expression_profiles", cancer, sample ,"", sep = "/")
  
  if (!dir.exists(exp_output_dir)){
    dir.create(exp_output_dir, recursive = TRUE)
  }
  
  expression_profiles <- read.csv(paste0("Expression_profiles/",cancer, " mRNA expression for immune markers (RNA Seq V2 RSEM).txt"), sep="\t")
  
  expression_profiles$samples_id <- str_sub(expression_profiles$SAMPLE_ID, 1, 12)
  
  # path from survival analysis script
  cdr3_scores <- read.csv(paste0(test_dir, cancer,"_",physicochemical_property,"_",analysis_category,"_survival.csv"))
  
  expression_profiles_subset <- expression_profiles[expression_profiles$samples_id %in% cdr3_scores$sample_id,]
  
  expression_profiles_subset <- expression_profiles_subset %>% merge(., cdr3_scores, by.x = "samples_id", by.y = "sample_id") %>% select(-STUDY_ID, -SAMPLE_ID)
  
  for (gene in genes) {
    
    p <- ggplot(expression_profiles_subset , aes(cdr3_score, expression_profiles_subset[,gene]))
    p + geom_boxplot() + ylab("mRNA expression") + labs(title = gene) + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(exp_output_dir, gene,"_",cancer,"_",physicochemical_property,"_",analysis_category,".png"))
    
  }
  
}




