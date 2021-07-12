library(openxlsx)
library(tidyverse)


for (cancer in cancers) {
  
  # input data
  expr_input_data <- paste(paste0(analysis_data, "_results"), cancer, sample, physicochemical_property,"", sep = "/")
  test_dir <-  paste0(expr_input_data,"test_output/")
  
  # if the genes are not provided, skip performing expression profiles analysis
  if (!exists("genes")){ 
    next
  }
 
  exp_output_dir <- paste(expression_profiles_dir, cancer, analysis_data, sample ,"", sep = "/")
  
  if (!dir.exists(exp_output_dir)){
    dir.create(exp_output_dir, recursive = TRUE)
  }
  
  expr_file <- paste0(expression_profiles_dir ,"/", cancer, " mRNA expression (RNA Seq V2 RSEM).txt")
  
  if (!file.exists(expr_file)){
    next
  }
  
  expression_profiles <- read.csv(expr_file, sep="\t")
  
  expression_profiles$samples_id <- str_sub(expression_profiles$SAMPLE_ID, 1, 12)
  
  # path from survival analysis script
  cdr3_scores <- read.csv(paste0(test_dir, cancer,"_",physicochemical_property,"_",analysis_category,"_survival.csv"))
  
  expression_profiles_subset <- expression_profiles[expression_profiles$samples_id %in% cdr3_scores$sample_id,]
  
  expression_profiles_subset <- expression_profiles_subset %>% merge(., cdr3_scores, by.x = "samples_id", by.y = "sample_id") %>% select(-STUDY_ID, -SAMPLE_ID)
  
  # for (gene in genes) {
  #   
  #   p <- ggplot(expression_profiles_subset , aes(cdr3_score, expression_profiles_subset[,gene]))
  #   p + geom_boxplot() + ylab("mRNA expression") + labs(title = gene) + theme(plot.title = element_text(hjust = 0.5))
  #   ggsave(paste0(exp_output_dir, gene,"_",cancer,"_",physicochemical_property,"_",analysis_category,".png"))
  #   
  # }
  
  cancer_ <- cancer
  gene_ <- c()
  t_statistic <- c()
  p_value <- c()
  conf.int.1 <- c()
  conf.int.2 <- c()
  
  
  
  for (i in 1:length(genes)){
    
    top_RSEM <- expression_profiles_subset[,c(genes[i], "cdr3_score")][expression_profiles_subset[,c(genes[i], "cdr3_score")]$cdr3_score == "top",][,genes[i]]
    
    bottom_RSEM <- expression_profiles_subset[,c(genes[i], "cdr3_score")][expression_profiles_subset[,c(genes[i], "cdr3_score")]$cdr3_score == "bottom",][,genes[i]]
    
    ttest <- t.test(top_RSEM, bottom_RSEM, var.equal = TRUE)
    
    gene_[i] <- genes[i]
    t_statistic[i] <- ttest$statistic[[1]]
    p_value[i] <- ttest$p.value
    conf.int.1[i] <- ttest$conf.int[1]
    conf.int.2[i] <- ttest$conf.int[2]
  }
  
  df <- data.frame(cancer_, gene_, t_statistic, p_value, conf.int.1, conf.int.2)
  
  write.csv(df, paste0(exp_output_dir, cancer,"_",physicochemical_property,"_",analysis_category,".csv"), row.names = FALSE)
  
}




