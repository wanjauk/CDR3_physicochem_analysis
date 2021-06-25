library(readr)
library(tidyverse)
library(survminer)
library(survival)
options(bitmapType='cairo')

source("get_clinical_data.R")
source("utils.R")

# input data directory
input_data <- paste("analysis_input","exome","physicochem_analysis_output/", sep = "/")
output_dir <- paste("analysis_input","exome","survival_analysis_figures/", sep = "/")
test_dir <-  paste("analysis_input","exome","test_output/", sep = "/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

if (!dir.exists(test_dir)){
  dir.create(test_dir)
}

# physicochemical_property <- "aromaticity"
# analysis_category <- "average"

physicochemical_property <- "aromaticity"
analysis_category <- "jacki_factor"

# physicochemical_property <- "uversky_hydropathy"
# analysis_category <- "average"

# physicochemical_property <- "uversky_hydropathy"
# analysis_category <- "jacki_factor"


cancers <- c("LUSC", "KIRC", "CESC", # no clinical data found for "COAD"
            "LUAD", "PAAD", "SKCM",
            "BLCA", "BRCA", "PRAD")



for (cancer in cancers){
  
  top50 <- read.csv(paste0(input_data, cancer,"_",physicochemical_property,"_",analysis_category,"_top50.csv"))
  cdr3_input <- read.csv(paste0(input_data, cancer,"_",physicochemical_property,"_",analysis_category,".csv"))
  
  clinical_data <- get_clinical_data(cancer = cancer)
  
  # names to lower case
  names(clinical_data) <- tolower(names(clinical_data))
  
  # convert blanks to NA
  clin_data <- clinical_data %>%
    dplyr::mutate_each(funs = funs(convert_blank_to_na), everything())
  
  clin_data$sample_id <- str_sub(clin_data$case_id, 1, 12)
  
  
  # subset the data to get only ids in jacki factor analysis
  clin_data <- clin_data[clin_data$sample_id %in% cdr3_input[,1],]
  
  clin_data <- clin_data %>% mutate(cdr3_score = ifelse(clin_data$sample_id %in% top50[,1], "top", "bottom"))
  
  clin_data <- clin_data %>% separate(os_status, c("vital_status_code", "vital_status"))
  clin_data$vital_status_code <- as.numeric(clin_data$vital_status_code)
  
  # Tabulate by outcome
  xtabs(~cdr3_score+vital_status_code, data=clin_data) %>% addmargins()
  
  sfit <- survfit(Surv(os_months, vital_status_code)~cdr3_score, data=clin_data)
  
  ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, title = paste0(cancer," ",physicochemical_property," ",analysis_category), ggtheme=custom_theme()) +
    labs(x = "Overall survival (months)") 

  ggsave(paste0(output_dir, cancer,"_",physicochemical_property,"_",analysis_category,"_overall_survival.png"))
  
  # # write out the data for inspection
  to_csv_file <- clin_data %>% select(sample_id, cdr3_score)
  write.csv(to_csv_file, paste0(test_dir, cancer,"_",physicochemical_property,"_",analysis_category,"_survival.csv"), row.names = FALSE)
  
}

