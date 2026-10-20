##----------------------------------------------------------------------------------
## Fucntion to perform paiewise Fisher's test between 2 columns with binary values
##----------------------------------------------------------------------------------

## Arguments:
#       1. df_data: data frame with the columns
#       2. col_1: colummn name of the first variable for fishers test
#       3. col_2: column name of the second variable for fishers test
#       4. col_name: the value that needs to be added to a column in the output file

## libraries
  library(tidyverse)
  library(tibble)

pairwise_fishers <- function(df_data, 
                               col_1= 1,
                               col_2= 2,
                               col_name = "name"){
    
    df <- df_data[, c(col_1, col_2)]
    colnames(df) <- c("Var1", "Var2")
    
    
      a = sum(df$Var1 == 1 & df$Var2 == 1)
      b = sum(df$Var1 == 1 & df$Var2 == 0)
      c = sum(df$Var1 == 0 & df$Var2 == 1)
      d = sum(df$Var1 == 0 & df$Var2 == 0)
      
      if (c == 0) c = 1
      if (d == 0) d = 1
      
      fish_mat <- matrix(data = c(a,b,c,d), byrow = TRUE, nrow = 2)
      f_test <- fisher.test(fish_mat)
      
      result_df <- data.frame(
        cell_type = col_name,
        p_val = f_test$p.value,
        odds_ratio = as.numeric(f_test$estimate),
        contingency_a = a,
        contingency_b = b,
        contingency_c = c,
        contingency_d = d)
      
      return(result_df)
      
    
    
}    
