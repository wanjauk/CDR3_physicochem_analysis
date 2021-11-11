library(tidyverse)
library(Rmisc)

genes <- c("CTLA4", "PDCD1", "PDCD10","PDCD1LG2", "LAG3")

cancers <- c("LUSC", "BLCA", "SKCM")
input_data_dir <- c("Exome_data")
analysis_data <- sub("\\_.*", "", input_data_dir)

expression_profiles_dir <- "Expression_profiles"

sample <- "tumor"

physicochemical_property <- "aromaticity"
analysis_category <- "average"

t_test_scores <- data.frame()

for (cancer in cancers){

        # input data
        expr_input_data <- paste(paste0(analysis_data, "_results"), cancer, sample, physicochemical_property,"", sep = "/")
        test_dir <-  paste0(expr_input_data,"test_output/")
        
        # exp_output_dir <- paste(expression_profiles_dir, cancer, analysis_data, sample ,"", sep = "/")
        
        expr_file <- paste0(expression_profiles_dir ,"/", cancer, " mRNA expression (RNA Seq V2 RSEM).txt")
        
        expression_profiles <- read.csv(expr_file, sep="\t")
        
        expression_profiles$samples_id <- str_sub(expression_profiles$SAMPLE_ID, 1, 12)
        
        
        # path from survival analysis script
        cdr3_scores <- read.csv(paste0(test_dir, cancer,"_",physicochemical_property,"_",analysis_category,"_survival.csv"))
        
        expression_profiles_subset <- expression_profiles[expression_profiles$samples_id %in% cdr3_scores$sample_id,]
        
        expression_profiles_subset <- expression_profiles_subset %>% merge(., cdr3_scores, by.x = "samples_id", by.y = "sample_id") %>% select(-STUDY_ID, -SAMPLE_ID)
        
        
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
        # t_test_scores <- rbind(t_test_scores, df)
}






skcm <- t_test_scores[t_test_scores$cancer_ == "SKCM",] %>% select(-cancer_, -t_statistic,-p_value)

keycol <- "conf.int"
valuecol <- "confidence_interval"
gathercols <- c("conf.int.1", "conf.int.2")

data_long <- gather_(skcm, keycol, valuecol, gathercols)
# data_long$conf.int <- as.factor(data_long$conf.int)

ggplot(data_long ,aes(x=gene_,y=confidence_interval,group=gene_)) +
  geom_boxplot(position = "identity")


ggplot(data=data_long, aes(x=gene_, y=confidence_interval, fill=conf.int)) +
  stat_boxplot( aes(gene_, confidence_interval), 
                geom='errorbar', linetype=1, width=0.5)+
  geom_boxplot(position = "identity") + 
  scale_fill_manual(values=c("black", "blue")) +
  geom_hline(yintercept=0, color="darkred", linetype="dotted") + stat_summary(fun.data = mean_se, geom = "errorbar")# +
  # theme(text=element_text(family="serif"), panel.background=element_blank(),
  #       axis.text.x=element_text(angle=90,hjust=1,vjust=.3))

ggsave("plot.png")






genes_of_interest <- expression_profiles_subset %>% select(cdr3_score, `genes`)

keycol_ <- "genes_"
valuecol_ <- "rsem_values"
gathercols_ <- genes

tg <- gather_(genes_of_interest, keycol_, valuecol_, gathercols_)

tgc <- summarySE(tg, measurevar="rsem_values", groupvars=c("cdr3_score","genes_"))

pd <- position_dodge(0.1)

ggplot(tgc, aes(x=genes_, y=rsem_values, group=genes_)) + 
  geom_errorbar(aes(ymin=rsem_values-ci, ymax=rsem_values+ci), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)












# z_stat <- (mean(top_RSEM, na.rm = TRUE) - mean(bottom_RSEM) - delta_0) / 
#   sqrt(var(top_RSEM, na.rm = TRUE) / n_1 + var(bottom_RSEM) / n_2)
# 
# 
# biology <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2, 4, 7, 2, 9)
# english <- c(8, 5, 4, 10, 4, 5, 7, 2, 6, 1, 2, 7, 0, 6, 4, 12, 5, 2)
# 
# delta_0 <- 0 # hypothesized mean of the population
# 
# # by assumption --> population variance
# sigma_sq_1 <- 3
# sigma_sq_2 <- 2
# 
# n_1 <- length(top_RSEM)
# n_2 <- length(bottom_RSEM)
# 
# # calculate the z-statistic
# z_stat <- (mean(biology) - mean(english) - delta_0) / 
#   sqrt(sigma_sq_1 / n_1 + sigma_sq_2 / n_2)
# 
# z_stat
# 
# 
# library(BSDA)
# ztest <- z.test(x = biology, y = english, alternative = "two.sided", mu = delta_0, sigma.x = sigma_sq_1,
#        sigma.y = sigma_sq_2, conf.level = 0.95)
# 
# ztest
# 
# x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
# 
# # calculating sample variance
# sample_variance <- var(x)
# sample_mean <- mean(x)
# 
# z.test(x = x, y = NULL, alternative = "two.sided", mu = 0, sigma.x = sample_variance,
#        sigma.y = NULL, conf.level = 0.95)
# 
# ttest <- t.test(biology, english, var.equal = TRUE)
# 
# 
# 
# 
# 
