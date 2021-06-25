# function to download clinical data from cbioportal

get_clinical_data <- function(cancer = cancer){
  
  cancer <- tolower(cancer)
  
  url <- paste0("http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id=",cancer,"_tcga_all")
  req <- httr::GET(url)
  clinical_data <- httr::content(req,
                                 type = 'text/tab-separated-values',
                                 col_names = T,
                                 col_types = NULL)
  return(clinical_data)
}

## convert empty strings -> NA values
convert_blank_to_na <- function(x) {
  if (!purrr::is_character(x)) {
    warning('input vector is not character - returning original input')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
