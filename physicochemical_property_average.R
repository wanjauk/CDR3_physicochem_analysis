#required package
library(tidyverse)

physicochemical_property_average <- function(cdr3_physicochem = cdr3_physicochem,
                                              receptor = receptor,
                                              sample_type = sample_type,
                                              physicochemical_property = physicochemical_property,
                                             analysis_strategy = analysis_strategy) {
  
  # filter the physicochemical properties to obtain receptors of interest and limit the data to primary and metastatic tumors
  physicochem <- subset(cdr3_physicochem, Sample.Type %in% sample_type & Receptor %in% receptor)
  
  # subset the receptors physicochemical properties to obtain only the physicochemical property of interest
  physicochem <- physicochem[,c("Filename", physicochemical_property)]
  
  # group by file name before finding average value of the physicochemical property
  physicochem <- physicochem %>% group_by(Filename)
  
  # get the average value for physicochemical property of interest
  physicochem_average <- summarise(physicochem, physicochem_property_avg = mean(.data[[physicochemical_property]]))
  
  # give the physicochemical_property_avg column an informative name
  physicochem_average <- physicochem_average %>% rename_at("physicochem_property_avg",~analysis_strategy)
  
  
  return(physicochem_average)
  
}



