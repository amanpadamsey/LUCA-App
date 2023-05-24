#this function uses the functions in functions.R and does the downstream processing
#this is also a ground for testing the whole workflow
source("functions.R")
source("libraries.R")


# name == main functionality in R -----------------------------------------


main <- function() {
  
  edds_datapath <<- "./tmdd testing/edds_test.csv" #comment out when finished testing
  flowjo_datapath <<- "./tmdd testing/"

  
}

if (getOption('run.main', default=TRUE)) {
  main()
}


# end main script ---------------------------------------------------------


#process edds file and give warnings if anything is amiss
edds_process <- function(edds_datapath) {
  
  edds <- edds_read(edds_datapath)
 
  dis_vals <-
    edds %>% 
    distinct(`Experiment date`, `Well number`, `Plate number`) %>%
    nrow() # find distinct rows for each distinct well 
  
  if(!(dis_vals == nrow(edds))) {
    dups <- edds %>% 
      group_by(`Experiment date`, `Well number`, `Plate number`) %>% 
      count() %>% filter(n>1)
    as.character(dups)
    
    validate(paste("Following entries were non-unique",capture.output(print(dups)), collapse = "\n"))
    }
  # shinyFeedback::feedbackWarning('edds', !(dis_vals == nrow(edds)), 'Duplicate rows found')
  
  dups <- edds %>% filter(`Tapir ID_unlabeled molecule (parent)` != 'MOCK') %>% 
    group_by(`Tapir ID_unlabeled molecule (parent)`) %>%  summarise(n = n()) %>% ungroup() %>% 
    distinct(n) %>% nrow() > 1 #test for unequal replicates for molecules
  
  if(dups){ #print and display those unequal duplicates
    
    dups <- edds %>% filter(`Tapir ID_unlabeled molecule (parent)` != 'MOCK') %>% 
      group_by(`Tapir ID_unlabeled molecule (parent)`) %>%  summarise(n = n()) %>% ungroup()
    
    message(paste("Following TAPIR IDS were unequal in replicates",
                  capture.output(print(dups)), collapse = "\n"))
    
  }
  
  # shinyFeedback::feedbackWarning('edds',dups,'Date cannot be parsed')
  
  return(edds)

}


# flowjo processing -------------------------------------------------------

flowjo_processing <- function(flowjo_datapaths) {
  
  fj <- flow_jo_clean(flowjo_datapaths)
  
  na_vals <- fj %>% complete.cases() %>% any() #check for missing values
  if (na_vals) {
    na_row <- fj[!(fj %>% complete.cases()),]
    #map(fj,is.na) %>% map(\(.x) unlist(.x) %>% which(isTRUE(.x))) %>% unlist() %>% as.character()
  }
  shinyFeedback::feedbackWarning('flowjo',na_vals,'Flowjo file contains missing values')
  req(!na_vals) 
  
  numeric_cols <- fj %>% ungroup() %>% 
    select(c("Cell count_total", "Cell count_morphology_live", "Cell count_morphology", "Geometric Mean : pHAb-A")) %>% 
    map_lgl(is.numeric) %>% 
    unlist() %>% 
    all() # check whether following columns are numeric
  
  shinyFeedback::feedbackWarning('flowjo',!numeric_cols,'Flowjo file contains non-numeric columns')
  req(numeric_cols) 
  
  
  
  fj %>% mutate(across(where(is.character), toupper))

}



