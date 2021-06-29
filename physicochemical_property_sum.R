library(tidyverse)

for (cancer in cancers){
  
  physicochemical_property_1 <- "aromaticity"
  input_dir <- paste(paste0(analysis_data, "_results"), cancer, sample, physicochemical_property_1,"", sep = "/")
  physicochemical_property_1 <- read.csv(file = paste0(input_dir, cancer,"_",physicochemical_property_1,"_jacki_factor.csv"))
  
  
  physicochemical_property_2 <- "uversky_hydropathy"
  input_dir <- paste(paste0(analysis_data, "_results"), cancer, sample, physicochemical_property_2,"", sep = "/")
  physicochemical_property_2 <- read.csv(file = paste0(input_dir, cancer,"_",physicochemical_property_2,"_jacki_factor.csv"))
  
  
  physicochemical_property_sum <- physicochemical_property_1 %>% 
    left_join(physicochemical_property_2, by = 'Tumor_Sample_ID') %>%
    mutate(physicochem_average_sum = aromaticity_average + uversky_hydropathy_average,
           physicochem_jacki_factor_sum = jacki_factor.x + jacki_factor.y) %>% 
    select(Tumor_Sample_ID, physicochem_average_sum, physicochem_jacki_factor_sum)
  
  output_dir <- paste(paste0(analysis_data, "_results"), cancer, sample, "physicochem","", sep = "/")
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  #-------------------------------------------------------------------------------------------------
  physicochem_average_sum <- physicochemical_property_sum %>% select(Tumor_Sample_ID, physicochem_average_sum)
  
  # For survival analysis:
  # Number of rows
  number_of_rows <- nrow(physicochem_average_sum)
  print(paste0(cancer," physicochem_average_sum:"))
  print(paste0("Number of samples for survival analysis: ", number_of_rows))
  
  # write out the  analysis output
  write.csv(physicochem_average_sum %>% arrange(desc(physicochem_average_sum)), file = paste0(output_dir, cancer,"_physicochem_average_sum.csv"), row.names = FALSE)
  
  # 1. Select the top n percent tumor samples
  top_physicochem_average_sum <- physicochem_average_sum %>% top_n(round(number_of_rows * select_n_percent))
  print(paste0("Number of samples in the top 50%: ", nrow(top_physicochem_average_sum)))
  
  # 2. Select the bottom n percent tumor samples
  bottom_physicochem_average_sum <- physicochem_average_sum  %>% top_n(-(number_of_rows - round(number_of_rows * select_n_percent)))
  
  print(paste0("Number of samples in the bottom 50%: ", nrow(bottom_physicochem_average_sum)))
  
  # write out the top and bottom n percent tumor samples.
  write.csv(top_physicochem_average_sum, file = paste0(output_dir, cancer,"_physicochem_average_sum_top50.csv"), row.names = FALSE)
  write.csv(bottom_physicochem_average_sum, file = paste0(output_dir, cancer,"_physicochem_average_sum_bottom50.csv"), row.names = FALSE)
  
  
  #-------------------------------------------------------------------------------------------------
  physicochem_jacki_factor_sum <- physicochemical_property_sum %>% select(Tumor_Sample_ID, physicochem_jacki_factor_sum)
  
  # For survival analysis:
  # Number of rows
  number_of_rows <- nrow(physicochem_jacki_factor_sum)
  print(paste0(cancer," physicochem_jacki_factor_sum:"))
  print(paste0("Number of samples for survival analysis: ", number_of_rows))
  
  # write out the  analysis output
  write.csv(physicochem_jacki_factor_sum %>% arrange(desc(physicochem_jacki_factor_sum)), file = paste0(output_dir, cancer,"_physicochem_jacki_factor_sum.csv"), row.names = FALSE)
  
  # 1. Select the top n percent tumor samples
  top_physicochem_jacki_factor_sum <- physicochem_jacki_factor_sum %>% top_n(round(number_of_rows * select_n_percent))
  print(paste0("Number of samples in the top 50%: ", nrow(top_physicochem_jacki_factor_sum)))
  
  # 2. Select the bottom n percent tumor samples
  bottom_physicochem_jacki_factor_sum <- physicochem_jacki_factor_sum  %>% top_n(-(number_of_rows - round(number_of_rows * select_n_percent)))
  
  print(paste0("Number of samples in the bottom 50%: ", nrow(bottom_physicochem_jacki_factor_sum)))
  
  # write out the top and bottom n percent tumor samples.
  write.csv(top_physicochem_jacki_factor_sum, file = paste0(output_dir, cancer,"_physicochem_jacki_factor_sum_top50.csv"), row.names = FALSE)
  write.csv(bottom_physicochem_jacki_factor_sum, file = paste0(output_dir, cancer,"_physicochem_jacki_factor_sum_bottom50.csv"), row.names = FALSE)
  
  
}