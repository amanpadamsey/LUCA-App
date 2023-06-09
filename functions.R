
# open model (reframe) ----------------------------------------------------

open_model <- function(df, column_model, ori_df) {
  ori_df <- ori_df %>% ungroup()
  
  df %>% reframe(broom::tidy(!!as.name(column_model)), broom::glance(!!as.name(column_model))) %>% rename_with(., function(x) {
    paste(x, column_model, sep = "__")
  }, .cols = everything()) %>% 
    bind_cols(df[, 1:ncol(df)-1], .) %>% 
    left_join(ori_df,.)
  
}


# minmax normalization ----------------------------------------------------

min_max_day <- function(df, control_mabs, col_est) {
  
  if (length(control_mabs) == 2) {
    
  mota <-
    df %>% 
    filter(`Tapir ID_unlabeled molecule (parent)` == control_mabs[1]) %>% 
    select(starts_with(col_est)) %>% .[1, 1] %>% as.numeric()
  
  cd_20 <-
    df %>% filter(`Tapir ID_unlabeled molecule (parent)` == control_mabs[2]) %>% 
    select(starts_with(col_est)) %>% .[1, 1] %>% as.numeric()
  
  df <- df %>% 
    mutate('min_max_{col_est}' := (!!as.name(col_est) - mota) / (cd_20 - mota))
  
  } else if (length(control_mabs) == 1) {
    
    mota <-
      df %>% 
      filter(`Tapir ID_unlabeled molecule (parent)` == control_mabs[1]) %>% 
      select(starts_with(col_est)) %>% .[1, 1] %>% as.numeric()
    
    df <- df %>% 
      mutate('min_max_{col_est}' := (!!as.name(col_est)) / (mota))
    
    
  } else {
    stop(' Please enter at least one control mab ')
  }
  
  df
  
}




# mock subtract -----------------------------------------------------------

mock_sub <- function(df,mfi_variable) {
  # enter grouped data to be subtracted
  mock <- mean((df |> 
                  filter(str_detect(`Tapir ID_unlabeled molecule (parent)`, 
                                    regex('^mock',ignore_case = TRUE))))[[mfi_variable]], na.rm=TRUE)
  
  df[[paste0(mfi_variable,'_BG subtracted')]] <- df[[mfi_variable]] - mock
  
  df
}


# Normalize points --------------------------------------------------------
# do linear regression of controls and find 'fitted' points for each control
point_norm  <- function(df,mfi_variable,control_mabs) {
  #group by experiment date and biosample id before this
  #column_of_interest <- colnames(df) %>% str_detect('subtracted_dose normalized') %>% colnames(df)[.]
  
  control_slopes <- map(control_mabs, \(x) df %>% 
                          filter(`Tapir ID_unlabeled molecule (parent)` == x) %>% 
                          select(estimate__lm_wo_normpoint) %>% 
                          .[1, 1] %>% as.numeric() )
  
  names(control_slopes) <- control_mabs
  
  if (length(control_mabs) == 2){
    
    df[paste0(mfi_variable, '_BG subtracted_dose_control normalized')] <- 
      (df[[paste0(mfi_variable, '_BG subtracted_dose normalized')]] - df[['Incubation time']] *
         control_slopes[[1]]) / (control_slopes[[2]] - control_slopes[[1]])
    
  } else if (length(control_mabs) == 1){
    
    df[paste0(mfi_variable, '_BG subtracted_dose_control normalized')] <- 
      (df[[paste0(mfi_variable, '_BG subtracted_dose normalized')]]) / control_slopes[[1]]
    
  }
  
  df
}


point_norm_2 <- function(df,control_mabs,mfi_variable) {
  mota <-
    df %>% filter(`Tapir ID_unlabeled molecule (parent)` == control_mabs[1]) %>% select(estimate__lm_wo_normpoint) %>% .[1, 1] %>% as.numeric()
  
  cd_20 <-
    df %>% filter(`Tapir ID_unlabeled molecule (parent)` == control_mabs[2]) %>% select(estimate__lm_wo_normpoint) %>% .[1, 1] %>% as.numeric()
  
  df %>% mutate(
    `Norm Geometric Mean : pHAb-A_BG subtracted_dose normalized` := (
      `Geometric Mean : pHAb-A_BG subtracted_dose normalized` - `Incubation time` *
        mota
    ) / (cd_20 - mota)
  )
  
}


# linear regression for all mabs-------------------------------------------------------

linear_reg_grp  <- function(df,mfi_variable) {
  #group by mabs before this
  model <- lm(df[[paste0(mfi_variable,'_BG subtracted_dose_control normalized')]] ~ 0 + df[['Incubation time']])
  
  model_open <- cbind(broom::glance(model)[,-c(4:9)],broom::tidy(model))
  
  #left_join(df ,model_open,by=c(mfi_variable,'Incubation time (hrs)'))
  cbind(df,model_open)
  
}



# convert to graphpad file ------------------------------------------------
graphpad <- function(df, mfi_variable) {
  df <- df |> ungroup()
  slope_tapir <- map(df |> distinct(`Biosample ID`,`Tapir ID_unlabeled molecule (parent)`, slope) |> 
    filter(!is.na(`Tapir ID_unlabeled molecule (parent)`)) %>%
    split(f = as.factor(.$`Biosample ID`)),
    
    ~pivot_wider(.x,names_from = `Tapir ID_unlabeled molecule (parent)`, values_from = 'slope')
  )
  
  names(slope_tapir) <- paste0('combined_norm_slope_',names(slope_tapir))
  
  # capture all raw data into graphpad file
  
  df_date_bio_id <- df %>% split(list(.$`Experiment date`,.$`Biosample ID`))
  
  vars_to_pzfx <- names(df)[names(df) |> grep(pattern=mfi_variable)]
  vals_to_pzfx <- list()
  
  
  for (var in vars_to_pzfx) {
    
    vals_to_pzfx_temp <- map(
      df_date_bio_id,
      
      ~ .x |>
        arrange(`Tapir ID_unlabeled molecule (parent)`) |>
        filter(!is.na(`Tapir ID_unlabeled molecule (parent)`)) |>
        select(`Incubation time`, `Tapir ID_unlabeled molecule (parent)`, var) |>
        pivot_wider(names_from = `Tapir ID_unlabeled molecule (parent)`, values_from = var) |>
        unnest_wider(col = -1, names_sep = '_')
    ) 
    
    names(vals_to_pzfx_temp) <- paste0(var,' ',names(vals_to_pzfx_temp)) 
    #vals_to_pzfx_temp <- map(vals_to_pzfx_temp,.x)
    
    vals_to_pzfx <- append(vals_to_pzfx,vals_to_pzfx_temp)
  }
  
  
 to_pzfx <- map(c(slope_tapir,vals_to_pzfx), ~ .x |> mutate_all(as.numeric))
 print(to_pzfx)
 to_pzfx
}


# test final workflow -----------------------------------------------------
# mfi_variable <- 'Geometric Mean : pHAb-A'
# 
# edds <- edds_read('~/rshiny test files/EDDS_dummy - EDDS_dummy 2.csv')
# flowjo <- flow_jo_clean(list("C:\\Users\\PADAMSEA\\Documents\\rshiny test files\\flowjo_files\\21-Dec-2022_2hr.xls",
#                              "C:\\Users\\PADAMSEA\\Documents\\rshiny test files\\flowjo_files\\21-Dec-2022_4hr.xls",
#                              "C:\\Users\\PADAMSEA\\Documents\\rshiny test files\\flowjo_files\\22-Dec-2022_2hr.xls"))
# 
# dosing <- dosing_sol_clean("C:\\Users\\PADAMSEA\\Documents\\rshiny test files\\NGCD3TCB Dosing solution day 2.xlsx")
# 
# 
# left_join(edds,
#           dosing,
#           by = c("Experiment date", "Tapir ID_unlabeled molecule (parent)")) |>
#   left_join(flowjo,
#             by = c('Plate number', 'Well number')) |>
#   mutate(`Viability` = `Cell count_morphology_live` / `Cell count_morphology`,
#               `Cell fraction_gated` = `Cell count_morphology_live` /`Cell count_total`) |>
#   group_by(`Experiment date`, `Biosample ID`) %>%
#   do(mock_sub(.,'Geometric Mean : pHAb-A')) |>
#   mutate(
#     '{mfi_variable}_BG subtracted_dose normalized' := .data[[paste0(mfi_variable, '_BG subtracted')]] /
#       `Fluorescence_dosing solution`
#   )  %>%
#   do(point_norm(.,control_mabs = 'CD44 var 1 PGLALA',mfi_variable = 'Geometric Mean : pHAb-A')) |> 
#   group_by(`Experiment date`, `Biosample ID`,`Tapir ID_unlabeled molecule (parent)`) -> df
# 
# 
# na_vals <- is.na(df[[paste0(mfi_variable,'_BG subtracted_dose_control normalized')]])
# df[!na_vals,]  %>%
#   do(linear_reg_grp(.,mfi_variable = mfi_variable)) |> bind_rows(df[na_vals,]) -> df
# 
# #   
#   
#   # do(point_norm(.)) |> View()

# df <- read_csv("C:\\Users\\PADAMSEA\\Downloads\\EDDS_dummy - EDDS_dummy 2 (6).csv")
  