#!/usr/bin/Rscript
#this function uses the functions in functions.R and does the downstream processing
#this is also a ground for testing the whole workflow
source("functions.R")
source("libraries.R")


# name == main functionality in R -----------------------------------------

if (sys.nframe() == 0){
  # edds_datapath <<- "./tmdd testing/edds_test.csv" #comment out when finished testing
  # flowjo_datapath_list <<- "./test flowjo export/" %>% list.files(full.names = TRUE)
  # dosing_datapath <<-  "./test tecan readout/" %>% list.files(full.names = TRUE) 
  # control_mabs <<- c("P1AF1537","P1AA4006")
  # mfi_choices <<- "Geometric Mean : pHAb-A"
  # 
  edds_datapath <<- "./testing Jo data/20230524_edds_nextCD3.csv" #comment out when finished testing
  flowjo_datapath_list <<- "./testing Jo data/" %>% list.files(pattern = "Flow",full.names = TRUE)
  dosing_datapath <<-  "./testing Jo data/" %>% list.files(full.names = TRUE,pattern = ".xlsx") 
  control_mabs <<- c("P1AF1537","P1AA4006")
  mfi_choices <<- "Geometric Mean : pHAb-A"
  
  
  
  
}


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
      showNotification(err)}, 
      error = function(err){cat("Could not parse date in EDDS, please
                                          check whether format is y-m-d")
        showNotification(err)}) # stop if date parsing error or warning
  
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
    
    showNotification(paste("Following TAPIR IDS were unequal in replicates",
                  capture.output(print(dups)), collapse = "\n"))
    
  }
  
  # shinyFeedback::feedbackWarning('edds',dups,'Date cannot be parsed')
  
  return(edds)

}




# flowjo cleaning ---------------------------------------------------------


flow_jo_clean <- function(flowjo_datapath_list) { #enter list of flowjo files (path)
  
  
  semi_clean <- lapply(flowjo_datapath_list,readxl::read_excel,col_types = c('text','text','numeric','numeric'),trim_ws = TRUE) |>  #read all flowjo files
    
    # map2(str_extract(list_flowjo_names,'\\d+-\\w+-\\d+') %>% lapply(.,parse_date,format='%d-%b-%Y'),\(x,y) x |>  mutate(`Experiment date` = y)) |>  #parse date and add column
    # map2(str_extract(list_flowjo_names,'(\\d)hr',group=1),\(x,y) x |>  mutate(`Incubation time (hrs)` = as.numeric(y))) |> #add Incubation time (hrs) from filename
    # 
    bind_rows(.id = "column_label") |> #combine read files into one big dataframe
    
    mutate(across(where(is.character), str_trim)) |> 
    
    mutate(`Well number` = str_extract(Name,'(...)\\.fcs$', group=1),
           `Plate number` = str_extract(Name,regex("plate.(\\w+)_.*_",ignore_case = TRUE), group=1)) |> #add well number and plate number
    
    fill(`Well number`,`Plate number`) #filldown 
  
  
  uni_depth <- unique(semi_clean$Depth) #find depth levels in df
  
  stats_flow <- grepl('=',semi_clean$Name) #find all rows with '='
  
  semi_clean$Var <- NA
  
  
  semi_clean$Var[stats_flow] <- semi_clean$Name[stats_flow] %>%  #find all rows with '='
    map(~unlist(str_split(.,' = '))[1]) |> unlist() #split strings and select variable name for new column Var
  
  
  semi_clean |> mutate(`Cell count_total` = case_when(is.na(Depth)~`#Cells`), 
                       `Cell count_morphology_live` = case_when(Depth == uni_depth[length(uni_depth)-1] ~ `#Cells`),
                       `Cell count_morphology` = case_when(Depth == uni_depth[length(uni_depth)-2] ~ `#Cells`),
                       `Plate number` = as.character(`Plate number`),
                       ultimate_gate = case_when(Depth == uni_depth[length(uni_depth)-1] ~ str_extract(Name, '[^/](\\w*)\\S$'))) |> 
    
    pivot_wider(names_from = 'Var',values_from = 'Statistic')  |> 
    
    group_by(`Well number`,`Plate number`) |> 
    
    fill(everything(), .direction='downup') |> 
    
    select(-c(Depth:`#Cells`,`NA`,column_label)) |> 
    
    unique() 
  
  
}






# flowjo processing -------------------------------------------------------




flowjo_processing <- function(flowjo_datapath_list) {
  
  fj <- flow_jo_clean(flowjo_datapath_list)
  
  na_vals <- fj %>% complete.cases() %>% all() #check for missing values, equals TRUE when no missing values present
  shinyFeedback::feedbackWarning("flowjo",!na_vals,"flowjo contains missing values")
  # if(na_vals) showNotification("hello there")
  if (!na_vals) {
    na_row <- fj[!(fj %>% complete.cases()),]
    #map(fj,is.na) %>% map(\(.x) unlist(.x) %>% which(isTRUE(.x))) %>% unlist() %>% as.character()
    showNotification(fj[fj %>% complete.cases() %>% !. ,c("Well number","Plate number")] %>% print() %>% capture.output() %>% paste(collapse = "\n"),
                     duration = NULL)
    
    # showNotification("Ensure the naming convention of the wells are as follows: Plate xyz_B4_B04")
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
    showNotification("Flowjo file cannot be parsed as numbers, please check format")
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


#process the dosing solution file and showNotification warnings

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

EDDS_combined_processing <- function(edds,mfi_choices,flowjo,dosing) {
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
    #check if there are no common elements 
    # if (match( 
    #   edds$`Tapir ID_unlabeled molecule (parent)`,
    #   dosing$`Tapir ID_unlabeled molecule (parent)`
    # ) %>%
    # is.na() %>% any()) {
    #   showNotification(paste(
    #     "Following rows had missing values",
    #     capture.output(print(na_row)),
    #     collapse = "\n"
    #   ))
    #   
    # }
    edds_combined <- left_join(edds,
                               dosing,
                               by = c("Experiment date", "Tapir ID_unlabeled molecule (parent)")) |>
      left_join(flowjo,
                by = c('Plate number', 'Well number')) #|> View())
  }
  
  na_vals_mfi <- edds_combined %>% filter(`Tapir ID_unlabeled molecule (parent)` != "MOCK") %>% 
    select(mfi_choices) %>% is.na(.) %>% any() #checks whether there are any NA values in the geomean mfi column, returns TRUE if none
  
  
  
  if(!na_vals_mfi) showNotification("combined file has NA values, please check your input files",duration = NULL)
  
  return(edds_combined)
}



# Mathematical manipulations on EDDS --------------------------------------




edds_analysis <- function(edds_combined,mfi_choices,control_mabs) {
  
  edds_dn <- edds_combined
  
  edds_dn <-  edds_dn |> 
    mutate(`Viability` = `Cell count_morphology_live` / `Cell count_morphology`,
           `Cell fraction_gated` = `Cell count_morphology_live` /`Cell count_total`) |> 
    group_by(`Experiment date`, `Biosample ID`, `Incubation time`) %>%
    do(mock_sub(.,mfi_choices)) |>
    mutate(
      '{mfi_choices}_BG subtracted_dose normalized' := .data[[paste0(mfi_choices, '_BG subtracted')]] /
        `Fluorescence_dosing solution`
    ) 
  
  
  edds_dn %>% filter(!str_detect(`Tapir ID_unlabeled molecule (parent)`,
                                 regex('mock', ignore_case = TRUE))) %>%
    group_by(`Experiment date`,
             `Tapir ID_unlabeled molecule (parent)`,
             `Biosample ID`) %>%
    do(lm_wo_normpoint = lm(.[[paste0(mfi_choices, '_BG subtracted_dose normalized')]] ~ 0 + .[['Incubation time']],data = .)
    ) -> model_wo.pointnorm
  
  edds_dn <- open_model(model_wo.pointnorm,'lm_wo_normpoint',edds_dn)
  
  
  edds_dn <- edds_dn %>% ungroup() %>% group_by(`Experiment date`,`Biosample ID`) %>%
    do(min_max_day(.,control_mabs,'estimate__lm_wo_normpoint')) %>% ungroup()
  
  
  
  edds_dn <-
    edds_dn %>% 
    group_by(`Experiment date`, `Biosample ID`) %>% 
    do(point_norm(., mfi_choices, control_mabs)) %>% ungroup()
  
  model_w.pointnorm <- edds_dn %>% filter(!str_detect(`Tapir ID_unlabeled molecule (parent)`,
                                                      regex('mock', ignore_case = TRUE))) %>% 
    group_by(`Tapir ID_unlabeled molecule (parent)`,`Biosample ID`) %>%
    do(lm_w_normpoint = lm(.[[paste0(mfi_choices, '_BG subtracted_dose_control normalized')]] ~ 0 + .[['Incubation time']],data = .))
  
  edds_dn <- open_model(model_w.pointnorm,'lm_w_normpoint',edds_dn)

}


# Write to file for copy in graphpad --------------------------------------

# create function to alternate columns 
blend_df <- function(df) {
  half <- df[,1:(ncol(df)/2)]
  other_half <- df[,((ncol(df)/2)+1):ncol(df)]
  neworder <- order(c(2*(seq_along(half) - 1) + 1,
                      2*seq_along(other_half)))
  cbind(half, other_half)[,neworder]
}


if(sys.nframe() == 0)type_lm = "wo"

slope_stderr_format <- function(edds_dn, type_lm) { 
  
  slope <- paste0("estimate__lm_",type_lm,"_normpoint")
  std_err <- paste0("std.error__lm_",type_lm,"_normpoint")
  
  slope_std_err <- edds_dn %>%  filter(`Tapir ID_unlabeled molecule (parent)` != "MOCK") %>%
    select(`Tapir ID_unlabeled molecule (parent)`,
           !!as.name(slope), !!as.name(std_err)) %>%
    distinct(`Tapir ID_unlabeled molecule (parent)`, .keep_all = TRUE) %>%
    pivot_wider(names_from = `Tapir ID_unlabeled molecule (parent)`, 
                values_from = c(!!as.name(slope), !!as.name(std_err))) %>% 
    blend_df()
  
  names(slope_std_err) <- map_chr(names(slope_std_err), ~ str_remove(.x,regex("(wo)|(w)")))

  return(slope_std_err)
}



slope_stderr_combine <- function(edds_dn) {
  
  indi_vals <- edds_dn %>% group_by(`Experiment date`) %>% 
    do(slope_stderr_format(.,"wo")) %>% 
    ungroup()
  indi_vals$`Experiment date` <- as.character(indi_vals$`Experiment date`)
  
  
  comb_vals <- edds_dn %>% do(slope_stderr_format(.,"w"))
  comb_vals <- comb_vals %>% mutate(`Experiment date` = "day1+day2")
  combined_vals <- bind_rows(indi_vals,comb_vals)
  
  names(combined_vals)[-1] <- map_chr(names(combined_vals[-1]), ~ str_extract(.x,"point_(.*)",group = 1))
  
  return(combined_vals)
}




normalized_days <- function(edds_dn) {

  
  edds_dn %>%  filter(`Tapir ID_unlabeled molecule (parent)` != "MOCK") %>% select(
    `Tapir ID_unlabeled molecule (parent)`,
    min_max_estimate__lm_wo_normpoint,
    `Experiment date`
  ) %>% distinct(`Tapir ID_unlabeled molecule (parent)`,
                 `Experiment date`,
                 .keep_all = TRUE) %>% pivot_wider(names_from = `Tapir ID_unlabeled molecule (parent)`,
                                                   values_from = min_max_estimate__lm_wo_normpoint,
                                                   names_sep = "_")
  
  
  
}

#raw data 

raw_data_edds <- function(edds_dn, mfi_choices, col_raw) {
  
  raw_data_edds_ <- edds_dn %>% ungroup() %>% arrange(`Experiment date`) %>% 
    filter(`Tapir ID_unlabeled molecule (parent)` != "MOCK") %>%
    select(
    `Tapir ID_unlabeled molecule (parent)`,
    `Incubation time`,
    # !!as.name(paste0(mfi_choices,"_BG subtracted")),
    !!as.name(col_raw)
  ) %>% pivot_wider(
    names_from = `Tapir ID_unlabeled molecule (parent)`,
    values_from = c(
      # !!as.name(paste0(mfi_choices,"_BG subtracted")),
      !!as.name(col_raw)
    ),values_fn = list) %>% unnest_wider(col = everything(), names_sep = "_")
  
  names(raw_data_edds_)[-1] <- map_chr(names(raw_data_edds_[-1]), ~ str_extract(.x,"(.*)_",group = 1))
  
  return(raw_data_edds_)
}


raw_data_edds_combine <- function(edds_dn, mfi_choices) {
  
  list(
    raw_data_edds(
      edds_dn,
      mfi_choices,
      paste0(mfi_choices, "_BG subtracted_dose normalized")
    ),
    raw_data_edds(
      edds_dn,
      mfi_choices,
      paste0(mfi_choices, "_BG subtracted_dose_control normalized")
    )
  )
  
}


#make big list of all data to be written to graphpad

files_to_write <- function(edds_dn, mfi_choices) { #returns list of files that can be graphpad friendly
  
  total_graphpad <- raw_data_edds_combine(edds_dn,mfi_choices) %>% 
    append(list(normalized_days(edds_dn))) %>% 
    append(list(slope_stderr_combine(edds_dn)))
  
  #total_graphpad %>% length()
  names(total_graphpad) <- c("rawdata_BG subtracted_dose normalized.csv",
                             "rawdata_BG subtracted_dose_control normalized.csv",
                             "")

}




