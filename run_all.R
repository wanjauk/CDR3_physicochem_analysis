
# cancers under analysis
# cancers <- c("LUSC", "KIRC", "CESC",
#             "LUAD", "PAAD", "COAD",
#             "BLCA", "BRCA", "PRAD",
#             "SKCM")

cancers <- "LUSC"

genes <- c("CD274", "CD33")

# input_data_dirs <- c("Exome_data", "RNA-seq_data")
input_data_dirs <- c("RNA-seq_data")

expression_profiles_dir <- "Expression_profiles"


for (input_data_dir in input_data_dirs){
  
  # get the name of the type of data
  analysis_data <- sub("\\_.*", "", input_data_dir)
  
  # parameters
  mutation_classification <- "Missense_Mutation"
  receptors <- c("TRA", "TRB")
  aromatic_aa <- c("F", "W", "Y")
  
  # 1. Tumor
  sample <- "tumor"
  sample_type <- c("Metastatic","Primary Tumor")
  
      ###################
      # Aromaticity
      ###################
      physicochemical_property <- "aromaticity"
      analysis_strategy <- "aromaticity_average"
      source("run_CDR3_physicochem_analysis.R")
      
      analysis_category <- "average"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      
      
      ##########################
      # Uversky Hydropathy
      ##########################
      physicochemical_property <- "uversky_hydropathy"
      analysis_strategy <- "uversky_hydropathy_average"
      source("run_CDR3_physicochem_analysis.R")
      
      analysis_category <- "average"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      # script for summing physicochem properties 
      source("physicochemical_property_sum.R")
      
      physicochemical_property <- "physicochem"
      analysis_category <- "average_sum"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor_sum"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
  
  # 2. Blood
  
  sample <- "blood"
  sample_type <- "Blood Derived Normal"
  
  if ((analysis_data == "RNA-seq") & (sample == "blood")) {
    next
  }
  
      ###################
      # Aromaticity
      ###################
      physicochemical_property <- "aromaticity"
      analysis_strategy <- "aromaticity_average"
      source("run_CDR3_physicochem_analysis.R")
      
      
      analysis_category <- "average"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      ##########################
      # Uversky Hydropathy
      ##########################
      physicochemical_property <- "uversky_hydropathy"
      analysis_strategy <- "uversky_hydropathy_average"
      source("run_CDR3_physicochem_analysis.R")
      
      analysis_category <- "average"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      # script for summing physicochem properties 
      source("physicochemical_property_sum.R")
      
      physicochemical_property <- "physicochem"
      analysis_category <- "average_sum"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
      analysis_category <- "jacki_factor_sum"
      source("run_survival_analysis.R")
      
      source("expression_profiles.R")
      
}

