library(tidyverse)

#  A function to calculate the jacki factor
calculate_jacki_factor <- function(mutect_missense_mutations = mutect_missense_mutations,
                                   cdr3_physicochem = cdr3_physicochem,
                                   sample_type = sample_type,
                                   receptor = receptor,
                                   aromatic_aa = aromatic_aa,
                                   physicochemical_property = physicochemical_property,
                                   analysis_strategy = analysis_strategy){
  
  # subset the mutect missense mutation to obtain aromatic amino acid mutants
  mutect_aromatic_aa <- mutect_missense_mutations[mutect_missense_mutations$Mutant_AA %in% aromatic_aa, ]
  
  # filter the physicochemical properties to obtain receptors of interest and limit the data to primary and metastatic tumors
  physicochem <- subset(cdr3_physicochem, Sample.Type %in% sample_type & Receptor %in% receptor)

  # subset the receptors physicochemical properties to obtain only the physicochemical property of interest
  physicochem <- physicochem[,c("Filename", physicochemical_property)]

  # group by file name before finding average value of the physicochemical property
  physicochem <- physicochem %>% group_by(Filename)

  # get the average value for physicochemical property of interest
  # physicochem_max <- slice_max(physicochem, order_by = physicochemical_property, n = 1)
  physicochem_average <- summarise(physicochem, physicochem_property_avg = mean(.data[[physicochemical_property]]))

  # filter out duplicated filenames (samples) if present
  physicochem_average <- subset(physicochem_average, !duplicated(physicochem_average[,"Filename"]))

  # add the counts of aromatic amino acids tumor samples IDs from mutect file alongside matching file names in physicochem_average (Recovered CDR3s)
  cdr3_mutect_overlap <- mutect_aromatic_aa %>% 
    count(Tumor_Sample_ID) %>% 
    merge(.,physicochem_average, by.x= "Tumor_Sample_ID", by.y="Filename") 
  
  # rename the column with counts of Tumor sample files with mutant amino acids
  cdr3_mutect_overlap <- cdr3_mutect_overlap %>% rename(aromatic_mutants_count = n)
  
  # multiply the counts of mutant tumor samples with maximum fraction of aromatic amino acids
  cdr3_mutect_overlap$jacki_factor <- cdr3_mutect_overlap[,"aromatic_mutants_count"] * cdr3_mutect_overlap[,"physicochem_property_avg"]
  
  # give the physicochemical_property_avg column an informative name
  cdr3_mutect_overlap <- cdr3_mutect_overlap %>% rename_at("physicochem_property_avg",~analysis_strategy)
  
  return(cdr3_mutect_overlap)
}





