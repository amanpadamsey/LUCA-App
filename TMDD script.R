#HUMAN TREM2 TMDD repeat
source("libraries.R")
source("workflow.R")
library(dplyr)
library(purrr)
library(tidyr)

#enter location of the flowjo files 
folder <-
  "Z:\\eADME\\Experiments\\LUCA\\results\\2023-05-08 Repeat huMacorphages TMDD Trem2 for Manuscript\\Titration + Donors Experiment\\Aman analysis using R"

flowjo_datapath_list <-
  folder %>% list.files(pattern = ".xls", full.names = TRUE)
flowjo <- flowjo_processing(flowjo_datapath_list)
flowjo <- flowjo %>% arrange(`Plate number`, `Well number`)

flowjo_original <- flowjo
#build edds from flojo file already prepared

mabs_tested <- "P1AF1094
P1AG4483
P1AE3306
P1AF8008
"
mabs_tested <-
  mabs_tested %>% str_extract_all("\\w+", simplify = TRUE) #%>% list()

pattern_mabs_1_2_5_6 <-
  map(mabs_tested[-4], \(x) append(rep(x, each = 9), "mock") %>% append(., .)) %>% unlist(recursive = FALSE)

flowjo[1:(2 * length(pattern_mabs_1_2)), "Tapir ID_unlabeled molecule (parent)"] <-
  pattern_mabs_1_2 %>% rep(2)
# flowjo %>% View()

pattern_mabs_3_4_7_8 <-
  map(mabs_tested[-3], \(x) append(rep(x, each = 9), "mock") %>% append(., .)) %>% unlist(recursive = FALSE)

part1 <- flowjo %>% filter(`Plate number` %in% c(1, 2, 5, 6))
part2 <- flowjo %>% filter(`Plate number` %in% c(3, 4, 7, 8))
part1[, "Tapir ID_unlabeled molecule (parent)"] <-
  pattern_mabs_1_2_5_6 %>% rep(4)
part2[, "Tapir ID_unlabeled molecule (parent)"] <-
  pattern_mabs_3_4_7_8 %>% rep(4)

flowjo <- part1 %>% bind_rows(part2)


#enter donor names

flowjo <-
  flowjo %>% mutate("Biosample ID" = if_else(`Plate number` %in% c(1, 2), "A", NA)) %>%
  mutate(`Biosample ID` = if_else(`Plate number` %in% c(5, 6), "C", `Biosample ID`))


flowjo <-
  flowjo %>% mutate(`Biosample ID` = if_else(
    str_detect(`Well number`, "F|G") &
      `Plate number` %in% c(7, 8),
    "C",
    `Biosample ID`
  )) %>%
  mutate(`Biosample ID` = if_else(
    str_detect(`Well number`,"(F|G)") &
      `Plate number` %in% c(3, 4),
    "A",
    `Biosample ID`
  ))

flowjo <-
  flowjo %>% mutate(`Biosample ID` = if_else(
    !str_detect(`Well number`, "F|G") &
      `Plate number` %in% c(7, 8),
    "D",
    `Biosample ID`
  )) %>%
  mutate(`Biosample ID` = if_else(
    !str_detect(`Well number`, "F|G") &
      `Plate number` %in% c(3, 4),
    "B",
    `Biosample ID`
  ))


#add incubation times
flowjo <-
  flowjo %>% mutate("Incubation time" = if_else(`Plate number` %in% c(1, 5, 3, 7), 10, 20))

#add quenched or not

pattern_quenched <- rep(c("Q", "NQ"), each = 10)
flowjo$Quenched <-
  rep(pattern_quenched, length(unique(flowjo$`Plate number`)) * 3)

flowjo <- flowjo %>% mutate(`Incubation compound 2` = if_else(`Tapir ID_unlabeled molecule (parent)` %in% mabs_tested[3:4],"10mg/ml IVIVG",NA))


#add concentration of labelled mab added to compound
flowjo[,"Fluorescent label"] <- "FITC"

pattern_conc <- c(200/(2.5**(0:8)),0)

flowjo <- flowjo %>% arrange(`Plate number`,`Well number`)
flowjo[,"Concentration_compound 1 (nM)"] <- pattern_conc %>% rep(48)

flowjo <- flowjo %>% mutate("Experiment date" = if_else(`Plate number` %in% c(1:4),ymd(20230531),ymd(20230601))) 


edds_hu_mac <- flowjo %>% mutate("Assay type" = "TMDD", "Project name" = "TREM2")
# write_csv(edds_hu_mac,paste0(folder,"\\20230601 EDDS_huMacrophage_repeat_2.csv"))



#retrieve ABC files
bead_vals <- folder %>% list.files(pattern = "bead",full.names = TRUE)
bead_vals <- bead_vals %>% 
  map( ~ read_csv(.x) %>%
         mutate(Date = mdy(Date))) %>%
  bind_rows()

bead_vals[-seq(1,nrow(bead_vals),4),] %>% ggplot(aes(log(ABC),log(MFI))) + geom_smooth(method = "lm") +
  geom_point() + 
  facet_grid(cols = vars(`TAPIR`), rows = vars(`Date`))



# fits  -------------------------------------------------------------------
# bead_valsr <- bead_vals %>% rename(`Geometric Mean : FITC-A` = MFI,
#                      `Experiment date` = Date,
#                      `Tapir ID_unlabeled molecule (parent)` = TAPIR)

# bead_vals <- bead_vals %>% mutate(Date_f = as.factor(Date))
# 
# fits_3 <- lmList(log(ABC) ~ log(MFI)*Date|TAPIR,data = bead_vals[-seq(1,nrow(bead_vals),4),])
# fits_4 <- lmList(log(ABC) ~ log(MFI)*Date|TAPIR,data = bead_vals)
# 
# bead_vals$fits3 <- predict(fits_3,bead_vals %>% select(MFI,TAPIR,Date_f)) %>% exp()
# bead_vals$fits4 <- predict(fits_4,bead_vals %>% select(MFI,TAPIR,Date_f)) %>% exp()
# 
# 
# 
# split_beads <- bead_vals[-seq(1,nrow(bead_vals),4),] %>% group_by(Date) %>% group_split()
# fits_3_1 <- lmList(log(ABC) ~ log(MFI)|TAPIR,data = split_beads[[1]])
# fits_3_2 <- lmList(log(ABC) ~ log(MFI)|TAPIR,data = split_beads[[2]])
# fits_3_2 %>% coef()
# fits_3_1 %>% coef()
# 
# predict(fits_3_1,bead_vals %>% filter(Date == split_beads[[1]]$Date[1]) %>% select(MFI,TAPIR)) %>% exp()
# predict(fits_3_2,bead_vals %>% filter(Date == split_beads[[2]]$Date[1]) %>% select(MFI,TAPIR)) %>% exp()
# 
ggplot(bead_vals,aes(ABC,fits3)) + geom_point() + geom_smooth(method = "lm")
ggplot(bead_vals,aes(ABC,fits4)) + geom_point() + geom_smooth(method = "lm")
# 
# 

# make function to do grouped lm ------------------------------------------
edds_hu_mac <- edds_hu_mac %>% rename(MFI = `Geometric Mean : FITC-A`,
                                      TAPIR = `Tapir ID_unlabeled molecule (parent)`,
                         
                                                   Date = `Experiment date`) %>% ungroup()
edds_hu_mac <- edds_hu_mac %>% mutate(TAPIR = if_else(TAPIR == "mock",
                                                      NA,TAPIR)) %>%
  fill(TAPIR,.direction = "down")


edds_hu_mac_mdoel <- bead_vals[-seq(1,nrow(bead_vals),4),] %>% group_by(Date,TAPIR) %>% 
  do(model = lm(log(ABC) ~ log(MFI),data = .)) %>% left_join(edds_hu_mac)

edds_hu_mac_mdoel <- edds_hu_mac_mdoel %>% mutate(pred_ind = predict(model,list("MFI"=MFI)) %>% exp()) 


# predictions -------------------------------------------------------------
#rename edds for predictions
# 
# 
# edds_hu_mac <- edds_hu_mac %>% ungroup() 
# # 
# # fits_3 <- lmList(log(ABC) ~ log(MFI)*Date|TAPIR,data = bead_vals[-seq(1,nrow(bead_vals),4),])
# # fits_4 <- lmList(log(ABC) ~ log(MFI)*Date|TAPIR,data = bead_vals)
# 


# edds_hu_mac[,"ABC_pred"] <- predict(fits_3,edds_hu_mac %>% select(MFI,Date,TAPIR)) %>% exp()
# 
# 
# 
# edds_hu_mac %>% ggplot(aes(ABC_pred,`Well number`)) + geom_col() +
#   facet_wrap("Date")


# subtract mock  ----------------------------------------------------------



  # enter grouped data to be subtracted
edds_hu_mac_mdoel_ms <- edds_hu_mac_mdoel %>% filter(`Concentration_compound 1 (nM)` == 0) %>% 
    group_by(`Biosample ID`,Date,TAPIR) %>% summarise(mean_mock = mean(pred_ind)) %>% 
  ungroup() %>% 
  left_join(edds_hu_mac_mdoel) %>% mutate(pred_ind_womock = pred_ind - mean_mock)
  
edds_hu_mac_mdoel_ms %>% ggplot(aes(pred_ind_womock,`Well number`)) + geom_col() +
  facet_wrap(vars(Date))


# estimate quenched vs non quenched ---------------------------------------

#formula : internalized = (quenched - (1-QE)*non_quenched)/QE

internalized_estimate <- function(df,qe){
  q <- df %>% filter(Quenched == "Q") %>% select(pred_ind_womock) %>% as.numeric()
  nq <- df %>% filter(Quenched == "NQ") %>% select(pred_ind_womock) %>% as.numeric()
  internalized <- (q-(1-qe)*nq)/qe
  df$internalized <- internalized
  df$surface <- nq - internalized
  df
  
}


edds_hu_mac_mdoel_ms <- edds_hu_mac_mdoel_ms %>% group_by(`Biosample ID`,Date,TAPIR,`Incubation time`,`Concentration_compound 1 (nM)`) %>%
  do(internalized_estimate(.,0.85)) %>% ungroup()

edds_hu_mac_mdoel_ms <- edds_hu_mac_mdoel_ms %>% mutate(internalized_per_time = internalized/`Incubation time`) 

edds_hu_mac_mdoel_ms %>% ggplot(aes(internalized_per_time,`Well number`)) + 
  geom_col() + facet_wrap(vars(Date))




# export ------------------------------------------------------------------
order <- "P1AF1094	P1AG4483	P1AE3306	P1AF8008" %>% str_split("\t",simplify = TRUE)
rshiny_folder <- "G:\\Shared drives\\PS_iSafe__in vitro ADME LM gDrive\\Projects\\TREM2\\Publication Figures & Data & Draft\\GraphPad Files\\Rshiny\\"
edds_hu_mac_mdoel_ms %>%  group_by(`Biosample ID`) %>% 
  distinct(`Biosample ID`,TAPIR,`Concentration_compound 1 (nM)`,`Incubation time`,.keep_all = TRUE) %>% 
  select(`Biosample ID`,TAPIR,`Concentration_compound 1 (nM)`,internalized_per_time) %>% 
  group_split() %>%
  map(~ .x %>% 
  pivot_wider(names_from = TAPIR,
    values_from = internalized_per_time) %>% 
    unnest_wider(col = -c(1:2),names_sep = "_") %>%
    select(c(1,2),starts_with(order)) %>% 
    write_csv(paste0(rshiny_folder,
                     .[1,1],".csv"))
)# %>% .[[1]] %>% View()




write_csv(edds_hu_mac_mdoel_ms,paste0(rshiny_folder,"\\edds_complete.csv"))



# fit modified mm equation ------------------------------------------------
edds_hu_mac_mdoel_ms %>% group_by(`Biosample ID`,TAPIR) %>%
  distinct(TAPIR,`Incubation time`,`Concentration_compound 1 (nM)`,.keep_all = TRUE) %>% 
  group_split() -> edds_hu_mac_mdoel_ms_split


fitting_nls <- function(df) {
  tryCatch(
    
  expr = {
  model_lm <-nls(internalized_per_time ~ (v * `Concentration_compound 1 (nM)`/ (b+`Concentration_compound 1 (nM)`)) + c*`Concentration_compound 1 (nM)`,
  data = df, start=list(v=1,b=1,c=1),algorithm = "port",lower = list(v=0.001,b=0.001,c=0.001)) 
  
  return(model_lm)
  },
  
  error = function(e) {
    
    model_lm <- (lm(internalized_per_time ~ `Concentration_compound 1 (nM)`, data = df))  
    return(model_lm)
  } 
  
  )
  
  # model_nls %>% broom::augment()
}
# edds_hu_mac_mdoel_ms_split[[4]] -> df

# map(edds_hu_mac_mdoel_ms_split, fitting_nls)
# 
# fitting_nls(edds_hu_mac_mdoel_ms_split[[4]])
# 
# 


 test_fitting_nls <- edds_hu_mac_mdoel_ms  %>% filter(Quenched == "Q") %>% group_by(TAPIR,`Biosample ID`) %>% 
  do(model_nls = fitting_nls(.)) %>% left_join(edds_hu_mac_mdoel_ms)

 
 df_row <- test_fitting_nls[1,]
 open_nlm <- function(df_row) {
   stopifnot(df_row %>% nrow() == 1) #stops function if more than one row passed
   # df_row <- as.data.frame(df_row)
   stopifnot(is.data.frame(df_row))
   tidy_m <- df_row$model_nls[[1]]  %>% broom::tidy(.)
   glance_m <- df_row$model_nls[[1]]  %>% broom::glance(.)
   # summnls <- summary(df_row$model_nls[[1]])
   # summnls$coefficients %>% as.data.frame()
   bind_cols(df_row,tidy_m,glance_m)
   # as.data.frame(df_row)
   # df_row %>% reframe(broom::tidy(model_nls), broom::glance(model_nls))# %>% 
     # bind_cols(df_row)
   # 
 }
 
 test_fitting_nls_opened <- test_fitting_nls %>%
   distinct(TAPIR, Date, `Biosample ID`, .keep_all = TRUE) %>% 
   group_by(row_number()) %>%  
   do(open_nlm(.)) #extra columns produced for lm model automatically assigned
   
 
 test_fitting_nls_combined <- test_fitting_nls %>% left_join(test_fitting_nls_opened, by = c("TAPIR","Date","Biosample ID"),
                                relationship = "many-to-many")
 
 test_fitting_nls_combined %>% write_csv(folder %>% paste0("\\edds_w_fits.csv"))
 # pmap(test_fitting_nls,open_nlm)
#  
#  test_fitting_nls %>% mutate()
#  
#  
#  test_fitting_nls$model_nls[[1]] %>% broom::tidy(.)
#  
#  test_fitting_nls %>% reframe(broom::tidy(model_nls), broom::glance(model_nls)) %>% rename_with(., function(x) {
#    paste(x, model_nls, sep = "__")
#  }, .cols = everything()) %>% 
#    bind_cols(test_fitting_nls[, 1:ncol(test_fitting_nls)-1], .) %>% 
#    left_join(ori_test_fitting_nls,.)
#  
# open_model(test_fitting_nls,"model_nls",edds_hu_mac_mdoel_ms)
