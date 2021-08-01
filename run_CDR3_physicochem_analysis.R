# required packages
library(openxlsx)
library(tidyverse)

# load functions that calculates the jacki factor and calculate physicochemical average
source("calculate_jacki_factor.R")
source("physicochemical_property_average.R")

# # select n percent of tumor samples for survival analysis
# select_n_percent <- 0.5


for (cancer in cancers) {
  # Load and prepare the required data----------------------------------------
  
  # create output directory
  output_dir <- paste(paste0(analysis_data, "_results"), cancer, sample, physicochemical_property,"", sep = "/")
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  #CDR3 physicochemical properties
  cdr3_physicochem <- read.xlsx(paste0(input_data_dir, "/","VDJ_Recoveries Boris ",cancer,".xlsx"), sheet = "physicochem")
  
  if (physicochemical_property == "aromaticity"){
    
    # remove case ids where physicochemical property is zero e.g. aromaticity is 0
    cdr3_physicochem <- cdr3_physicochem[!cdr3_physicochem$aromaticity == 0,]
  }
  
  # load cancer mutect file
  mutect <- read.xlsx(paste0(input_data_dir, "/", "mutect Boris ",cancer,".xlsx"))
  
  # exclude row with missing NAs in the cdr3 physicochemical properties data
  cdr3_physicochem <- cdr3_physicochem[complete.cases(cdr3_physicochem),]

  # filter to obtain missense mutations in the mutect file
  mutect_missense_mutations <- mutect[mutect$Variant_Classification == mutation_classification,]
  
  # separate Amino acid into wild type and mutant amino acids in the mutect file
  mutect_missense_mutations <- mutect_missense_mutations %>% separate(Amino_acids, c("Wild_Type_AA", "Mutant_AA"))
  
  # obtain tumor sample IDs from Tumor_Sample_Barcode column in the mutect file
  mutect_missense_mutations$Tumor_Sample_ID <- str_sub(mutect_missense_mutations$Tumor_Sample_Barcode, 1, 12) # find better code to do this
  
  #----------------------------------------------------------------------------
  # calculate the jacki factor
  jacki_factor <- calculate_jacki_factor(mutect_missense_mutations = mutect_missense_mutations,
                                              cdr3_physicochem = cdr3_physicochem,
                                              sample_type = sample_type,
                                              receptor = receptors,
                                              aromatic_aa = aromatic_aa,
                                              physicochemical_property = physicochemical_property,
                                              analysis_strategy = analysis_strategy)

  
  # For survival analysis:
  # Number of rows
  number_of_rows <- nrow(jacki_factor)
  print(paste0(cancer," jacki factor:"))
  print(paste0("Number of samples for survival analysis: ", number_of_rows))
  
  # write out the jacki factor analysis output
  write.csv(jacki_factor %>% arrange(desc(jacki_factor)), file = paste0(output_dir, cancer,"_",physicochemical_property,"_jacki_factor_",select_n_percent * 100,"_",(1 - select_n_percent) * 100,".csv"), row.names = FALSE)
  
  # 1. Select the top n percent tumor samples
  top_jacki_factor <- jacki_factor %>% top_n(round(number_of_rows * select_n_percent))
  print(paste0("Number of samples in the top ",select_n_percent * 100,"%: ", nrow(top_jacki_factor)))
  
  # 2. Select the bottom n percent tumor samples
  bottom_jacki_factor <- jacki_factor  %>% top_n(-(number_of_rows - round(number_of_rows * select_n_percent)))
  
  print(paste0("Number of samples in the bottom ",(1 - select_n_percent) * 100,"%: ", nrow(bottom_jacki_factor)))
  
  # write out the top and bottom n percent tumor samples.
  write.csv(top_jacki_factor, file = paste0(output_dir, cancer,"_",physicochemical_property,"_jacki_factor_top",select_n_percent * 100,".csv"), row.names = FALSE)
  write.csv(bottom_jacki_factor, file = paste0(output_dir, cancer,"_",physicochemical_property,"_jacki_factor_bottom",(1 - select_n_percent)*100,".csv"), row.names = FALSE)
  
  
  ##----------------------------------------------------------------------------------------------
  # Calculate the physicochemical property average alone
  
  # physicochem_average <- physicochemical_property_average(cdr3_physicochem = cdr3_physicochem,
  #                                                         receptor = receptors,
  #                                                         sample_type = sample_type,
  #                                                         physicochemical_property = physicochemical_property,
  #                                                         analysis_strategy = analysis_strategy)
  physicochem_average <- jacki_factor %>% select(1,3)
  physicochem_average <- physicochem_average %>% arrange(desc(physicochem_average[,2]))

  number_of_rows <- nrow(physicochem_average)
  print(paste0(cancer," ",analysis_strategy,":"))
  print(paste0("Number of samples for survival analysis: ", number_of_rows))


  # write out the analysis output
  write.csv(physicochem_average, file = paste0(output_dir,cancer,"_",analysis_strategy,"_",select_n_percent * 100,"_",(1 - select_n_percent) * 100,".csv"), row.names = FALSE)

  # 1. Select the top n percent tumor samples
  top_physicochem_average <- physicochem_average %>% top_n(round(number_of_rows * select_n_percent))
  print(paste0("Number of samples in the top ",select_n_percent * 100,"%: ", nrow(top_physicochem_average)))

  # 2. Select the bottom 50 percent tumor samples
  bottom_physicochem_average <- physicochem_average %>% top_n(-(number_of_rows - round(number_of_rows * select_n_percent)))
  print(paste0("Number of samples in the bottom ",(1 - select_n_percent) * 100,"%: ", nrow(bottom_physicochem_average)))

  # write out the top and bottom 50 percent tumor samples.
  write.csv(top_physicochem_average, file = paste0(output_dir, cancer,"_",analysis_strategy,"_top",select_n_percent * 100,".csv"), row.names = FALSE)
  write.csv(bottom_physicochem_average, file = paste0(output_dir, cancer,"_",analysis_strategy,"_bottom",(1 - select_n_percent) * 100,".csv"), row.names = FALSE)

}

