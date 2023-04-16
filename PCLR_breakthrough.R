# Requires: breakthrough_sample_soma.rds, breakthrough_proteins.rds, Breakthrough_PCLR_min/,
# and `Breakthrough_PCLR_1SE/`.

# Run PCLR on 1 random bootstrap resample of breakthrough_sample_soma.rds dataset.
# The results after multiple repetitions are to be compiled in two directories: 
# `Breakthrough_PCLR_min` and `Breakthrough_PCLR_1SE`. 

ii <- commandArgs()[6]

library(tidyverse)
library(clogitL1)

select <- dplyr::select

#' Boostrap a dataset while maintaining strata structure
#' @param df Dataframe of data to bootstrap
#' @param strat Name of stratum variable
#' @param seed Seed for reproducibility
#' @return Dataframe of bootstrap resampled data
strata_bootstrap <- function(df, strat, seed=NA){
  set.seed(seed)
  strat_levels <- df %>% pull(all_of(strat)) %>% unique()
  
  # Split train and test sets
  shuffled_i <- sample(strat_levels, replace=T)
  
  df_list <- list()
  for (i in 1:length(shuffled_i)){
    df_list[[i]] <- df %>%
      filter(get(strat) == shuffled_i[i]) %>%
      mutate(!!rlang::sym(strat) := i)
  }
  bs_set <- do.call("rbind", df_list)
  bs_set <- bs_set %>% mutate(!!rlang::sym(strat) := !!rlang::sym(strat))
  
  return(bs_set)
}


#' Extract non-zero coefficients from cross-validated clogitL1 object
#' @param cv_out cv.clogitL1 object
#' @param prot_df Dataframe of analytes used to label the non-zero coefficients
#' @param s Which model to extract coefficients from. Usually "min" or "1se_alt"
#' @return Dataframe of labelled, non-zero coefficients
clogit_nzc <- function(cv_out, prot_df, s="min"){
  cv_smry <- summary.cv.clogitL1(cv_out)
  
  if (s=="min"){
    prot_df$Coef <- cv_smry$beta_minCV
  } else if (s=="1se"){
    prot_df$Coef <- cv_smry$beta_minCV1se
  } else if (s=="1se_alt"){
    # Best model 1SE away, but with MORE predictors
    adj_dev <- which(min(cv_out$mean_cv) - (cv_out$mean_cv - cv_out$se_cv) > 0)# -
    alt.1se.i <- adj_dev[length(adj_dev)]
    cv_out$nz_beta[alt.1se.i]
    
    prot_df$Coef <- cv_out$beta[alt.1se.i,]
  }
  
  prot_df <- prot_df %>%
    filter(Coef != 0) %>%
    select(c(AptName,
             TargetFullName,
             UniProt,
             EntrezGeneSymbol,
             Coef))
  
  return(prot_df)
}


#' Main function used to retrieve non-zero coefficients from a (bootstrapped) set of data.
#' @param sample_soma Proteomic dataset. Requires columns "Group" and "GroupId".
#' @return List of two dataframes from s="min" and s="1se_alt" from clogit_nzc.
#'         Additional columns are provided for alternative coefficient interpretation. 
main_func <- function(sample_soma, seed){
  set.seed(seed)
  # Run initial CLR
  pclr_out <- clogitL1(
    y=sample_soma %>% pull(Group),
    x=sample_soma %>% select(-all_of(c("Group", "GroupId"))),
    strata=sample_soma %>% pull(GroupId)
  )
  
  # Run cross-validation CLR
  cv_out <- cv.clogitL1(pclr_out)
  
  # Get cross-validated coefficients
  beta_df <- data.frame(cv_out$beta)
  names(beta_df) <- names(sample_soma)[3:dim(sample_soma)]
  
  # 1SE Coefficients
  nzcs_1se <- clogit_nzc(cv_out, proteins, s="1se_alt")
  nzcs_1se <- nzcs_1se %>% mutate(
    Odds=exp(Coef),
    UnitRiskPercentage=(exp(Coef) - 1)*100
  )
  
  # Min coefficients
  nzcs_min <- clogit_nzc(cv_out, proteins, s="min")
  nzcs_min <- nzcs_min %>% mutate(
    Odds=exp(Coef),
    UnitRiskPercentage=(exp(Coef) - 1)*100
  )
  
  out <- list(nzcs_min, nzcs_1se)
  names(out) <- c("nzcs.min", "nzcs.1se")
  
  return(out)
}


# Read in pre-prepared data (1)
sample_soma <- readRDS("./breakthrough_sample_soma.rds")
sample_soma <- sample_soma %>% mutate(Group=case_when(
  Group == "Control" ~ 0,
  Group == "Case" ~ 1))
# Read in pre-prepared data (2)
proteins <- readRDS("./breakthrough_proteins.rds")

# Produce a bootstrap dataset
soma_bs <- strata_bootstrap(sample_soma, "GroupId", seed=as.integer(ii))

out <- main_func(soma_bs, seed=as.integer(ii)) # Same seed as CLI index

write.csv(out$nzcs.min, paste0("./Breakthrough_PCLR_min/", toString(ii), ".csv"), row.names=F)
write.csv(out$nzcs.1se, paste0("./Breakthrough_PCLR_1SE/", toString(ii), ".csv"), row.names=F)
