#read quickcal
quickcal <- 
  readxl::read_excel('./16459_BD.xls')[7:10, 3:4] %>% mutate(across(everything(), as.numeric))
colnames(quickcal) <- c('abc', 'bead_fl')
#make model
model_abc <- lm(log(abc) ~ log(bead_fl), data = quickcal)

#predict input file
to_conv_list <-
  list.files('/Users/padamsea/Downloads/20230421 bangs lab quick cal/',
             pattern = 'to_conv',
             full.names = TRUE)

to_conv_list_read <-
  map(to_conv_list,
      \(S) read_csv(S, col_names = FALSE) %>%
        mutate(across(everything(), as.numeric)))


# do math -----------------------------------------------------------------






to_conv_list_read_converted <-
  to_conv_list_read %>% map(\(S) S %>% map(\(v) predict(model_abc, newdata = list('bead_fl' = v)) %>% exp())) %>% map(as.data.frame)


# Display





# # check -------------------------------------------------------------------
# temp_directory <- tempdir()
# dir.create(tempdir())
# 
# to_conv_list_read_converted %>%  imap(function(x, y) {
#   if (!is.null(x)) {
#     file_name <- glue("{y}_data.csv")
#     readr::write_csv(x, file.path(temp_directory, file_name), col_names = FALSE)
#   }
# })
# 
# 
# 

# calculations with edds file ---------------------------------------------
bead_fl_test <- read_csv('bead_fl_test.csv')

models_mabs <- map(bead_fl_test, \(S) lm(log(quickcal$abc) ~ log(S)))

edds_abc <-
  edds_read(edds_path = '20230420_edds_maskedaso_abc testing.csv')



to_abc <- function(edds_abc_split, models_mabs) {
  
  trans_df <- data.frame()
  
  for (df in edds_abc_split) {
    if (df$`Tapir ID_unlabeled molecule (parent)`[1] %in% names(models_mabs)) {
      abc <-
        predict(models_mabs[[df$`Tapir ID_unlabeled molecule (parent)`[1]]],
                list(x = df$`Geometric Mean : pHAb-A`)) %>% exp()
      
      df$abc <- abc
      
      trans_df <- trans_df %>% bind_rows(df)
    } else {
      trans_df <- trans_df %>% bind_rows(df)
    }
    
  }

  return(trans_df)
}


edds_abc <- edds_abc %>% 
  group_by(`Tapir ID_unlabeled molecule (parent)`) %>%
  group_split() %>%
  to_abc(.,models_mabs)

write_csv(edds_abc,'edds_abc.csv')


# fit Thomas' modified equation -------------------------------------------
read_delim(clipboard(),locale = locale(decimal_mark = ',')) %>% 
  datapasta::dpasta()




non_blocked_t30 <- tibble::tribble(
  ~ Concentration..nM.,
  ~ Mota.wt.,
  ~ P1AF1101,
  ~ P1AH0116,
  200,
  33802.67414,
  45624.54393,
  48785.20615,
  80,
  12776.46665,
  21731.89628,
  25839.40518,
  32,
  5004.047234,
  12494.32852,
  14643.28543,
  12.8,
  3529.115086,
  7425.60194,
  9746.472514,
  5.12,
  1555.45862,
  5771.62831,
  5861.954641,
  2.048,
  631.595361,
  3084.711493,
  4360.33194,
  0.8192,
  371.2504239,
  2979.979887,
  1728.350378,
  0.32768,
  577.3842272,
  1331.407381,
  270.5399167,
  0.131072,
  -27.40174409,
  472.6906299,
  1245.620064
)

library(drc)

m1 <- drm(P1AF1101 ~ Concentration..nM. + Concentration..nM. * e, data = non_blocked_t30, fct = MM.2())


plot(m1, log = '', pch = 17, main = "Fitted MM curve")



mm <- structure(list(S = c(3.6, 1.8, 0.9, 0.45, 0.225, 0.1125, 3.6, 
                           1.8, 0.9, 0.45, 0.225, 0.1125, 3.6, 1.8, 0.9, 0.45, 0.225, 0.1125, 
                           0), v = c(0.004407692, 0.004192308, 0.003553846, 0.002576923, 
                                     0.001661538, 0.001064286, 0.004835714, 0.004671429, 0.0039, 0.002857143, 
                                     0.00175, 0.001057143, 0.004907143, 0.004521429, 0.00375, 0.002764286, 
                                     0.001857143, 0.001121429, 0)), .Names = c("S", "v"), class = "data.frame", row.names = c(NA, 
                                                                                                                              -19L))

model.drm <- drm (v ~ S, data = non_blocked_t30, fct = MM.2())

model.nls <- nls(v ~ Vm * S/(K+S) + e*S, data = mm, 
                 start = list(K = max(mm$v)/2, Vm = max(mm$v), e=0))


library(nls2)
model.nls2 <- nls2(v ~ Vm * S/(K+S) , data = non_blocked_t30, start = list(K = max(mm$v)/2, Vm = max(mm$v)) )

selfStart(model = )
  # 
# model.nls <- nls(P1AF1101 ~ Vm * Concentration..nM./(K+Concentration..nM.) + e*Concentration..nM., data = non_blocked_t30, 
#                  start = list(K = max(non_blocked_t30$P1AF1101)/2, Vm = max(non_blocked_t30$P1AF1101), e=0))

colnames(non_blocked_t30)[c(1,3)] <- c('S','v')

summary(model.drm)
summary(model.nls)
summary(model.nls2)

# plot regression ---------------------------------------------------------
model_test <- models_mabs[[1]]


for (model_test in models_mabs) {
  model_test %>% broom::augment() %>% ggplot(aes(`log(quickcal$abc)`, `log(x)`)) +
    geom_point() +
    geom_smooth(method = 'lm')
  
}



# test lm -----------------------------------------------------------------


fit  <- lm(v~S, data = non_blocked_t30)
#second degree
fit2 <- lm(v~poly(S,2,raw=TRUE),data = non_blocked_t30)
#third degre,data = non_blocked_t30e
fit3 <- lm(v~poly(S,3,raw=TRUE),data = non_blocked_t30)
#fourth degre,data = non_blocked_t30e
fit4 <- lm(v~poly(S,4,raw=TRUE),data = non_blocked_t30)
#generate range of 50 numbers starting from 30 and ending at 160
xx <- seq(1,250, length=30)
plot(non_blocked_t30$S,non_blocked_t30$v,pch=19,ylim=c(0,60000))
lines(xx, predict(fit, data.frame(S=xx)), col="red")
lines(xx, predict(fit2, data.frame(S=xx)), col="green")
lines(xx, predict(fit3, data.frame(S=xx)), col="blue")
lines(xx, predict(fit4, data.frame(S=xx)), col="purple")
