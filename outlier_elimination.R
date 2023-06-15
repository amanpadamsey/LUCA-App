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
  
  
  col <- "Geometric Mean : pHAb-A_BG subtracted_dose normalized"
  df <- edds_dn %>% group_by(`Tapir ID_unlabeled molecule (parent)`,`Incubation time`,`Biosample ID`) %>% 
    group_split() %>% .[[1]]
  p <- 0.2
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
