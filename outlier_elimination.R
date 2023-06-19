# OUTLIER ELIMINATION TESTING
library(outliers)

if(sys.nframe() == 0){
  data_test <- c(42.542, 100, 109, 90, 95, 98, 99, 89)
  
  
  test <- grubbs.test(c(100, 109, 90, 1000006))
  test$method
  
  
  chisq.out.test(data_test)
  
  out_method  <- function(method, ...) {
    method(data_test, ...)
  }
  
  out_method(chisq.out.test)
  
  out_method(dixon.test)
  out_method(grubbs.test)
  
  
  # col <- "Geometric Mean : pHAb-A_BG subtracted_dose normalized"
  # df <- edds_dn %>% group_by(`Tapir ID_unlabeled molecule (parent)`,`Incubation time`,`Biosample ID`) %>% 
  #   group_split() %>% .[[1]]
  p <- 0.6
  df <- mtcars %>% mutate(Outlier_detected = 'Y')
  col <- 'hp'
  df <- df %>% mutate(SrNo = 1:nrow(df))
  df_orm <- outlier_rm(df,col,0.6)
  df$Outlier_detected[df_orm$SrNo] <- 'N'
  # where(is.na(df$Outlier_detected))
}

# recursive fucntion to remove outliers until grubbs test is below p value
outlier_rm <- function(df, col, p) {
  # stopifnot("column not found in df" == col  %in% colnames(mtcars))
  if (sum(!is.na(df[[col]])) < 3) #return df if number of elements in col is less than 3 
    return(df)
  
  # p <- 0.2
  # df <- df
  
  out_test <-
    df[[col]] %>%
    grubbs.test()
  
  if (out_test$p.value < p) { 
    if ("lowest" %in% str_split(out_test$alternative, " ", simplify = TRUE)) {
      df <-
        df[-which.min(df[[col]]), ]
      # print("jello_min")
      # outlier_rm(df, col, p)
    } else {
      
      df <-
        df[-which.max(df[[col]]), ]
      # print("jello_max")
      
    }
    
    # print("jello1")
    
    
    return(outlier_rm(df, col, p))
    # return(df)
  } else return(df)
  
  
  # print("jello2")
  
}


# edds_dn %>% group_by(`Tapir ID_unlabeled molecule (parent)`,`Incubation time`,`Biosample ID`) %>% 
#   group_split() -> test_edds
# outlier_rm(edds_dn,col,0.2)

if (sys.nframe() == 0) {
  library(MASS)
  library(foreign)
  
  cdata <- read.dta('C:\\Users\\padamsea\\Downloads\\crime.dta')
  rr.bisquare <- rlm(crime ~ poverty + single, data=cdata, psi = psi.bisquare)
  lm_fit <- lm(crime ~ poverty + single, data=cdata)
  
  summary(rr.bisquare)
  
  
  biweights <- data.frame(state = cdata$state, resid = rr.bisquare$resid, weight = rr.bisquare$w)
  biweights2 <- biweights[order(rr.bisquare$w), ]
  biweights2[1:15, ]
  plot(rr.bisquare)
  plot(lm_fit)
  
  
  edds_dn_1 <- edds_dn %>% filter(`Tapir ID_unlabeled molecule (parent)` == 'P1AJ6195-010')
  edds_dn_1 %>% ggplot(aes(`Incubation time`,`Geometric Mean : pHAb-A_BG subtracted_dose normalized`)) +
    geom_point(aes(col=as.factor(`Biosample ID`)))
  
  edds_dn_13 <- edds_dn %>% filter(`Tapir ID_unlabeled molecule (parent)` == 'P1AJ6195-010',
                                   `Biosample ID` == '3')
  
  
  rr.bisquare <- rlm(edds_dn_1$`Geometric Mean : pHAb-A_BG subtracted_dose normalized` ~  0 + edds_dn_1$`Incubation time` + edds_dn_1$`Biosample ID`, psi = psi.bisquare)
  lm_fit <- lm(edds_dn_1$`Geometric Mean : pHAb-A_BG subtracted_dose normalized` ~  0 + edds_dn_1$`Incubation time` + edds_dn_1$`Biosample ID`)
  plot(rr.bisquare)
  plot(lm_fit)
  biweights <- cbind.data.frame(edds_dn_1,resid = rr.bisquare$resid, weight = rr.bisquare$w)
  biweights2 <- biweights[order(rr.bisquare$w), ]
  
  rr.bisquare <- rlm(edds_dn_13$`Geometric Mean : pHAb-A_BG subtracted_dose normalized` ~  0 + edds_dn_13$`Incubation time` + edds_dn_13$`Biosample ID`, psi = psi.bisquare)
  lm_fit <- lm(edds_dn_13$`Geometric Mean : pHAb-A_BG subtracted_dose normalized` ~  0 + edds_dn_13$`Incubation time` + edds_dn_13$`Biosample ID`)
  plot(rr.bisquare)
  plot(lm_fit)
  biweights <- cbind.data.frame(edds_dn_13,resid = rr.bisquare$resid, weight = rr.bisquare$w)
  biweights2 <- biweights[order(rr.bisquare$w), ]
  
  }
