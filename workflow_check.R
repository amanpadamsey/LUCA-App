library(datapasta)
# Testing workflow using Johannes data ------------------------------------


# flowjo cleaning ---------------------------------------------------------

  
  #add error if NA in any of the columns
  fj <- flow_jo_clean(list("drive-download-20230419T055152Z-001/15-Mar-2023.xls","drive-download-20230419T055152Z-001/16-Mar-2023.xls"))
  # fj %>% ungroup()
  fj_cc <- fj %>% complete.cases() 
  fj_cc %>% all() #determine if there are NA vals anywhere, FALSE if NA somewhere
  which(!fj_cc) %>% unlist() %>%  paste(collapse  = ', ') #paste all NA value rows together
  
  
  
  
  fj %>% distinct(`Well number`,`Plate number`) %>% nrow() == nrow(fj) # check for uniqueness of each row
  fj %>% group_by(`Well number`,`Plate number`) %>% summarise(count = n()) %>% (.$count > 1) #
  
  which((test$count > 1))
  
  
  fj %>% ungroup() %>% 
    select(c("Cell count_total", "Cell count_morphology_live", "Cell count_morphology", "Geometric Mean : pHAb-A")) %>% 
    map_lgl(is.numeric) %>% 
    unlist() %>% 
    all() #test whether specific columns are numeric or not
  
  fj <- fj %>% mutate(across(where(is.character), toupper)) #capitalize all words
  

# dosing solution  --------------------------------------------------------

  ds <- dosing_sol_clean("drive-download-20230419T055152Z-001/NGCD3TCB Dosing solution d1.xlsx") %>% 
  bind_rows(dosing_sol_clean("drive-download-20230419T055152Z-001/NGCD3TCB Dosing solution d2.xlsx"))  
    
  ds$`Fluorescence_dosing solution` %>% is.numeric()
  ds$`Experiment date` %>% is.Date()
  
  ds %>% group_by(`Tapir ID_unlabeled molecule (parent)`) %>%  summarise(n = n()) %>% ungroup() %>% 
     distinct(n) %>% nrow() > 1 # TRUE if you are missing some tapir ids in either day of dosing sol

    
  ds %>% 
    group_by(`Tapir ID_unlabeled molecule (parent)`,`Experiment date`) %>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    distinct(n) %>% nrow() > 1 # TRUE if data repeated for group mentioned
  

# EDDS --------------------------------------------------------------------

edds <- edds_read("./tmdd testing/edds_abc.csv")

all(is.Date(edds$`Experiment date`), is.numeric(edds$`Incubation time`))    #check these two imp parameters

edds <- edds %>% mutate(across(where(is.character), toupper))
  
edds %>% distinct(`Experiment date`,`Well number`,`Plate number`) %>% nrow() == nrow(edds) #test whether all rows are distinct

edds %>% filter(`Tapir ID_unlabeled molecule (parent)` != 'MOCK') %>% 
  group_by(`Tapir ID_unlabeled molecule (parent)`) %>%  summarise(n = n()) %>% ungroup() %>% View()
  # distinct(n,.keep_all = TRUE) %>% nrow() > 1 


edds %>% mutate(`Plate number` = as.character(`Plate number`)) -> edds
# stitch together ---------------------------------------------------------


edds %>% left_join(fj) %>% left_join(ds) -> days_combined#write_csv('20230419_edds_aso.csv')


# math --------------------------------------------------------------------
#check whether there are mocks for both time points

mfi_variable <- 'Geometric Mean : pHAb-A'
edds$'Fluorescence_dosing solution' <- 30000

edds_dn <-  edds |> 
  mutate(`Viability` = `Cell count_morphology_live` / `Cell count_morphology`,
         `Cell fraction_gated` = `Cell count_morphology_live` /`Cell count_total`) |> 
  group_by(`Experiment date`, `Biosample ID`, `Incubation time`) %>%
  do(mock_sub(.,mfi_variable)) |>
  mutate(
    '{mfi_variable}_BG subtracted_dose normalized' := .data[[paste0(mfi_variable, '_BG subtracted')]] /
      `Fluorescence_dosing solution`
  ) 


edds_dn %>% filter(!str_detect(`Alias_labeled molecule`,
                               regex('mock', ignore_case = TRUE))) %>%
  group_by(`Experiment date`,
           `Tapir ID_unlabeled molecule (parent)`,
           `Biosample ID`) %>%
  do(lm_wo_normpoint = lm(.[[paste0(mfi_variable, '_BG subtracted_dose normalized')]] ~ 0 + .[['Incubation time']],data = .)
                          ) -> model_wo.pointnorm

edds_dn <- open_model(model_wo.pointnorm,'lm_wo_normpoint',edds_dn)


edds_dn <- edds_dn %>% ungroup() %>% group_by(`Experiment date`,`Biosample ID`) %>%
  do(min_max_day(.,control_mabs,'estimate__lm_wo_normpoint')) %>% ungroup()



edds_dn <-
  edds_dn %>% 
  group_by(`Experiment date`, `Biosample ID`) %>% 
  do(point_norm(., mfi_variable, control_mabs)) %>% ungroup()

model_w.pointnorm <- edds_dn %>% filter(!str_detect(`Alias_labeled molecule`,
                              regex('mock', ignore_case = TRUE))) %>% 
  group_by(`Tapir ID_unlabeled molecule (parent)`,`Biosample ID`) %>%
  do(lm_w_normpoint = lm(.[[paste0(mfi_variable, '_BG subtracted_dose_control normalized')]] ~ 0 + .[['Incubation time']],data = .))

edds_dn <- open_model(model_w.pointnorm,'lm_w_normpoint',edds_dn)

# write_csv(edds_dn,'edds_test.csv')



# final check -------------------------------------------------------------
edds_dn <- read_csv('./tmdd testing/edds_abc.csv')
edds_dn[['Fluorescence_dosing solution']] <- 30000
control_mabs <- c('hfhg')
mfi_variable <- 'Geometric Mean : pHAb-A'


edds_dn <-  edds_dn |> 
  mutate(`Viability` = `Cell count_morphology_live` / `Cell count_morphology`,
         `Cell fraction_gated` = `Cell count_morphology_live` /`Cell count_total`) |> 
  group_by(`Experiment date`, `Biosample ID`, `Incubation time`) %>%
  do(mock_sub(.,mfi_variable)) |>
  mutate(
    '{mfi_variable}_BG subtracted_dose normalized' := .data[[paste0(mfi_variable, '_BG subtracted')]] /
      `Fluorescence_dosing solution`
  ) 


edds_dn %>% filter(!str_detect(`Alias_labeled molecule`,
                               regex('mock', ignore_case = TRUE))) %>%
  group_by(`Experiment date`,
           `Tapir ID_unlabeled molecule (parent)`,
           `Biosample ID`) %>%
  do(lm_wo_normpoint = lm(.[[paste0(mfi_variable, '_BG subtracted_dose normalized')]] ~ 0 + .[['Incubation time']],data = .)
  ) -> model_wo.pointnorm

edds_dn <- open_model(model_wo.pointnorm,'lm_wo_normpoint',edds_dn)


edds_dn <- edds_dn %>% ungroup() %>% group_by(`Experiment date`,`Biosample ID`) %>%
  do(min_max_day(.,control_mabs,'estimate__lm_wo_normpoint')) %>% ungroup()



edds_dn <-
  edds_dn %>% 
  group_by(`Experiment date`, `Biosample ID`) %>% 
  do(point_norm(., mfi_variable, control_mabs)) %>% ungroup()

model_w.pointnorm <- edds_dn %>% filter(!str_detect(`Alias_labeled molecule`,
                                                    regex('mock', ignore_case = TRUE))) %>% 
  group_by(`Tapir ID_unlabeled molecule (parent)`,`Biosample ID`) %>%
  do(lm_w_normpoint = lm(.[[paste0(mfi_variable, '_BG subtracted_dose_control normalized')]] ~ 0 + .[['Incubation time']],data = .))

edds_dn <- open_model(model_w.pointnorm,'lm_w_normpoint',edds_dn)

