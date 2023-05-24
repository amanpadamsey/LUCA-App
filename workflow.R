#this function uses the functions in functions.R and does the downstream processing
#this is also a ground for testing the whole workflow
source("functions.R")
source("libraries.R")


# name == main functionality in R -----------------------------------------


main <- function() {
  
  edds_datapath <<- "./tmdd testing/edds_test.csv" #comment out when finished testing
  flowjo_datapath_list <<- "./test flowjo export/" %>% list.files(full.names = TRUE) 
  dosing_datapath <<-  "./test tecan readout/" %>% list.files(full.names = TRUE) 
  control_mabs <<- c("P1AF1537","P1AA4006")
  mfi_choices <<- "Geometric Mean : pHAb-A"
}

if (getOption('run.main', default=TRUE)) {
  main()
}


# end main script ---------------------------------------------------------

# read and process edds file ----------------------------------------------


edds_read <- function(edds_datapath){
  read_delim(edds_datapath,
             trim_ws = TRUE, na =  c("", "NA")) |>
    mutate(across(where(is.character), str_trim),
          `Plate number` = as.character(`Plate number`)) |> 
    mutate(`Experiment date` = ymd(`Experiment date`), #parse date
           `Tapir ID_unlabeled molecule (parent)` = 
             toupper(`Tapir ID_unlabeled molecule (parent)`),
           `Incubation time` = as.numeric(`Incubation time`)) %>% 
    tryCatch(warning = function(err){cat("Could not parse date in EDDS, please
                                          check whether format is y-m-d")
      message(err)}, 
      error = function(err){cat("Could not parse date in EDDS, please
                                          check whether format is y-m-d")
        message(err)}) # stop if date parsing error or warning
  
}



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

flowjo_processing <- function(flowjo_datapath_list) {
  
  fj <- flow_jo_clean(flowjo_datapath_list)
  
  na_vals <- fj %>% complete.cases() %>% any() #check for missing values
  if (na_vals) {
    na_row <- fj[!(fj %>% complete.cases()),]
    #map(fj,is.na) %>% map(\(.x) unlist(.x) %>% which(isTRUE(.x))) %>% unlist() %>% as.character()
    message(paste("Following rows had missing values",
                  capture.output(print(na_row)), collapse = "\n"))
    message("Ensure the naming convention of the wells are as follows: Plate xyz_B4_B04")
  }
  # shinyFeedback::feedbackWarning('flowjo',na_vals,'Flowjo file contains missing values')
  # req(!na_vals) 
  
  # numeric_cols <- fj %>% ungroup() %>% 
  #   select(c("Cell count_total", "Cell count_morphology_live", "Cell count_morphology", "Geometric Mean : pHAb-A")) %>% 
  #   map_lgl(is.numeric) %>% 
  #   unlist() %>% 
  #   all() # check whether following columns are numeric
  fj <- fj %>% select(-ultimate_gate) #remove unnecesary column from data
  
  if (fj[,-c(1,2)] %>% map(is.numeric) %>% unlist() %>% all() %>% !.) {
    message("Flowjo file cannot be parsed as numbers, please check format")
  }
  
  
  # shinyFeedback::feedbackWarning('flowjo',!numeric_cols,'Flowjo file contains non-numeric columns')
  # req(numeric_cols) 
  
  
  
  flowjo <- fj %>% mutate(across(where(is.character), toupper)) #make everything upper case
  
  return(flowjo)
}


# dosing solution processing ----------------------------------------------

dosing_sol_clean <- function(dosing_datapath) {
  
  #lapply(dosing_sol_file_path, function)
  
  dosing_sol_file_read <- readxl::read_excel(dosing_datapath,skip = 30,trim_ws = TRUE,
                                             col_names = c('well','Fluorescence_dosing solution','Tapir ID_unlabeled molecule (parent)')) |> 
    mutate(across(where(is.character), str_trim)) |> 
    mutate(`Tapir ID_unlabeled molecule (parent)` = toupper(`Tapir ID_unlabeled molecule (parent)`))
  
  #|> mutate(`Experiment date` = case_when(well == 'End Time:' ~ parse_date(`Tapir ID_unlabeled molecule (parent)`,str_extract(,'\\d+/\\d+/\\d+'),'%m/%d/%Y') ))
  dosing_sol_file_read['Experiment date'] <- mdy(str_extract(dosing_sol_file_read[nrow(dosing_sol_file_read),2],'\\d+/\\d+/\\d+')) 
  
  dosing_sol_file_read$`Fluorescence_dosing solution` <- as.numeric(dosing_sol_file_read$`Fluorescence_dosing solution`) 
  dosing_sol_file_read[!is.na(dosing_sol_file_read[,3]),] %>% select(-well)
}


#process the dosing solution file and message warnings

dosing_processing <- function(dosing_datapath) {
  ds <- lapply(dosing_datapath, dosing_sol_clean) |> bind_rows(.id = "column_label")
  
  dups <- ds %>% group_by(`Tapir ID_unlabeled molecule (parent)`) %>%  summarise(n = n()) %>% ungroup() %>% 
    distinct(n,.keep_all = TRUE) #test for uneven duplicates in dataframe
  if(nrow(dups)>1) validate("TAPIR names in names duplicated")
  if(ds$`Fluorescence_dosing solution`%>% is.numeric() %>% !.) validate("could not parse numeric in dosing file")
  
  # shinyFeedback::feedbackWarning('dosing',nrow(dups)==1,'Uneven duplicates in TAPIR ID')
  # shinyFeedback::feedbackWarning('dosing',!(ds$`Fluorescence_dosing solution` %>% is.numeric()),'Dosing solution not numeric')
  # print(!(ds$`Experiment date` %>% is.Date()))
  # shinyFeedback::feedbackWarning('dosing',(ds$`Experiment date` %>% is.na() %>% unlist() %>% any()),'Date cannot be parsed')
  
  # 
  # dis_vals <- ds %>% distinct(`Experiment date`,`Tapir ID_unlabeled molecule (parent)`) %>% nrow == nrow(ds)
  # shinyFeedback::feedbackWarning('dosing',!dis_vals,'Duplicate values found')
  # 
  dosing <- ds %>% mutate(across(where(is.character), toupper)) #capitalize all words
}


# appending edds file with flowjo and dosing ------------------------------

EDDS_combined_processing <- function(edds,mfi_choices) {
  if (c(mfi_choices, 'Fluorescence_dosing solution') %in% colnames(edds) %>% all()) {
    edds_combined <- edds
    
  } else if (mfi_choices %in% colnames(edds)) {
    req(dosing)
    edds_combined <- left_join(
      edds,
      dosing,
      by = c("Experiment date", "Tapir ID_unlabeled molecule (parent)")
    )
    
  } else {
    req(flowjo, dosing)
    edds_combined <- left_join(edds,
                               dosing,
                               by = c("Experiment date", "Tapir ID_unlabeled molecule (parent)")) |>
      left_join(flowjo,
                by = c('Plate number', 'Well number')) #|> View())
  }
  
  return(edds_combined)
}

  



