library(readxl)
library(tidyverse)
library(viridis)
library(bio3d)
library(writexl)

analyse_uncertain_new_subs <- function(structure,date){
  
  tryCatch({setwd(paste0("/Users/user/Documents/Case\ study\ paper/3mo_loos_ensemble/case_study_3moloos_unaligned_50/spike_filtered_",date,"/Master_tables/",structure))
  # read in predictions made by ensemble model
  alltrain_count <- read.csv("rf_50_ensemble_counts.csv")
  alltrain_count_probs <- read.csv("rf_50_ensemble_count_probs.csv")
  
  # check if this date includes 12 month test set
  test_sets <- list.files(pattern="*TEST_month12_*")
  if(length(test_sets)==0){
    full_tests<-F
  } else {full_tests<-T}
  
  
  if(full_tests==T){
# if we have full test sets (including 12 month timepoint) we analyse all --------
  # read in test sets -------------------------------------------------------
  # read in 1 month test set
  # read in master table of scores
  master_conservation_1 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month1_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_1$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_1_50 <- master_conservation_1 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  rm(master_conservation_1)
  
  # read in 3 month test set
  # read in master table of scores
  master_conservation_3 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month3_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_3$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_3_50 <- master_conservation_3 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  rm(master_conservation_3)
  
  # read in 6 month test set
  # read in master table of scores
  master_conservation_6 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month6_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_6$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_6_50 <- master_conservation_6 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  rm(master_conservation_6)
  
  # read in 12 month test set
  # read in master table of scores
  master_conservation_12 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month12_",date,".xlsx"))
  
  # calculate quantiles of observed probability to define classes later
  quantiles <- quantile(master_conservation_12$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_12_50 <- master_conservation_12 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted"))
  
  rm(master_conservation_12)
  
  # filter for individual predict sets
  alltrain_count_probs_1 <- alltrain_count_probs %>% 
    filter(Group=="1month")
  
  alltrain_count_probs_3 <- alltrain_count_probs %>% 
    filter(Group=="3month")
  
  alltrain_count_probs_6 <- alltrain_count_probs %>% 
    filter(Group=="6month")
  
  alltrain_count_probs_12 <- alltrain_count_probs %>% 
    filter(Group=="12month")
  
  # 1 month
  # make a new tibble with just replacement column to join by and the observed class of test set
  master_conservation_1_50$Replacement <- paste0(aa321(master_conservation_1_50$WTAA),master_conservation_1_50$Chain,master_conservation_1_50$Resno,master_conservation_1_50$ReplacementAA)

  # join to the probability results from the ensemble model
  alltrain_count_probs_1 <- left_join(alltrain_count_probs_1,master_conservation_1_50)
  
  # 3 month
  # make a new tibble with just replacement column to join by and the observed class of test set
  master_conservation_3_50$Replacement <- paste0(aa321(master_conservation_3_50$WTAA),master_conservation_3_50$Chain,master_conservation_3_50$Resno,master_conservation_3_50$ReplacementAA)

  # join to the probability results from the ensemble model
  alltrain_count_probs_3 <- left_join(alltrain_count_probs_3,master_conservation_3_50)
  
  # 6 month
  # make a new tibble with just replacement column to join by and the observed class of test set
  master_conservation_6_50$Replacement <- paste0(aa321(master_conservation_6_50$WTAA),master_conservation_6_50$Chain,master_conservation_6_50$Resno,master_conservation_6_50$ReplacementAA)

  
  # join to the probability results from the ensemble model
  alltrain_count_probs_6 <- left_join(alltrain_count_probs_6,master_conservation_6_50)
  
  # 12 month
  # make a new tibble with just replacement column to join by and the observed class of test set
  master_conservation_12_50$Replacement <- paste0(aa321(master_conservation_12_50$WTAA),master_conservation_12_50$Chain,master_conservation_12_50$Resno,master_conservation_12_50$ReplacementAA)

  # join to the probability results from the ensemble model
  alltrain_count_probs_12 <- left_join(alltrain_count_probs_12,master_conservation_12_50)
  
  # split the replacement column into composite parts to inspect character of the substitutions
  alltrain_count_probs_1 <- alltrain_count_probs_1 %>% 
    separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
    select(-X)
  
  alltrain_count_probs_3 <- alltrain_count_probs_3 %>% 
    separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
    select(-X)
  
  alltrain_count_probs_6 <- alltrain_count_probs_6 %>% 
    separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
    select(-X)
  
  alltrain_count_probs_12 <- alltrain_count_probs_12 %>% 
    separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
    select(-X)
  
  
  # class switching subs --------------------------------------------------------------
  
  # filter for class switching substitutions
  # read in all train sets
  trainset_reader <- function(count){
    train_set <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TRAIN_loosample_",count,"_",date,".xlsx"))
    
    # calculate quantiles of observed probability to define classes later
    quantiles <- quantile(train_set$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    # remove all potentially troublesome punctuation from column names
    train_set <- train_set %>% 
      dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)
    
    train_set <- train_set %>% 
      dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)
    
    train_set <- train_set %>% 
      dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)
    
    train_set <- train_set %>% 
      dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)
    
    train_set <- train_set %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted"))
    
    train_set$Replacement <- paste0(aa321(train_set$WTAA),train_set$Chain,train_set$Resno,train_set$ReplacementAA)
    
    train_set <- train_set %>% 
      arrange(Replacement)
    
    return(train_set)
  }
  
  no_samples <- seq(1:10)
  all_train <- lapply(no_samples,trainset_reader)
  
  # make new tables for only new subs
  alltrain_count_probs_new_1 <- alltrain_count_probs_1
  
  alltrain_count_probs_new_3 <- alltrain_count_probs_3
  
  alltrain_count_probs_new_6 <- alltrain_count_probs_6
  
  alltrain_count_probs_new_12 <- alltrain_count_probs_12
  
  alltrain_count_probs_new_1 <- alltrain_count_probs_new_1 %>% 
    arrange(Replacement) %>% 
    filter(Observed_class != all_train[[1]]$Observed_class | 
           Observed_class != all_train[[2]]$Observed_class |
           Observed_class != all_train[[3]]$Observed_class |
           Observed_class != all_train[[4]]$Observed_class |
           Observed_class != all_train[[5]]$Observed_class |
           Observed_class != all_train[[6]]$Observed_class |
           Observed_class != all_train[[7]]$Observed_class |
           Observed_class != all_train[[8]]$Observed_class |
           Observed_class != all_train[[9]]$Observed_class |
           Observed_class != all_train[[10]]$Observed_class)
  
  alltrain_count_probs_new_3 <- alltrain_count_probs_new_3 %>% 
    arrange(Replacement) %>% 
    filter(Observed_class != all_train[[1]]$Observed_class | 
             Observed_class != all_train[[2]]$Observed_class |
             Observed_class != all_train[[3]]$Observed_class |
             Observed_class != all_train[[4]]$Observed_class |
             Observed_class != all_train[[5]]$Observed_class |
             Observed_class != all_train[[6]]$Observed_class |
             Observed_class != all_train[[7]]$Observed_class |
             Observed_class != all_train[[8]]$Observed_class |
             Observed_class != all_train[[9]]$Observed_class |
             Observed_class != all_train[[10]]$Observed_class)
  
  alltrain_count_probs_new_6 <- alltrain_count_probs_new_6 %>% 
    arrange(Replacement) %>% 
    filter(Observed_class != all_train[[1]]$Observed_class | 
             Observed_class != all_train[[2]]$Observed_class |
             Observed_class != all_train[[3]]$Observed_class |
             Observed_class != all_train[[4]]$Observed_class |
             Observed_class != all_train[[5]]$Observed_class |
             Observed_class != all_train[[6]]$Observed_class |
             Observed_class != all_train[[7]]$Observed_class |
             Observed_class != all_train[[8]]$Observed_class |
             Observed_class != all_train[[9]]$Observed_class |
             Observed_class != all_train[[10]]$Observed_class)
  
  alltrain_count_probs_new_12 <- alltrain_count_probs_new_12 %>% 
    arrange(Replacement) %>% 
    filter(Observed_class != all_train[[1]]$Observed_class | 
             Observed_class != all_train[[2]]$Observed_class |
             Observed_class != all_train[[3]]$Observed_class |
             Observed_class != all_train[[4]]$Observed_class |
             Observed_class != all_train[[5]]$Observed_class |
             Observed_class != all_train[[6]]$Observed_class |
             Observed_class != all_train[[7]]$Observed_class |
             Observed_class != all_train[[8]]$Observed_class |
             Observed_class != all_train[[9]]$Observed_class |
             Observed_class != all_train[[10]]$Observed_class)
  
  full_stats
  new_stats
  
  
  full_means_1 <- c(mean(alltrain_count_probs_1$DeltaDeltaG),mean(alltrain_count_probs_1$total_score),mean(alltrain_count_probs_1$LogLikelihood),mean(alltrain_count_probs_1$ESST_probability))
  new_means_1 <- c(mean(alltrain_count_probs_new_1$DeltaDeltaG),mean(alltrain_count_probs_new_1$total_score),mean(alltrain_count_probs_new_1$LogLikelihood),mean(alltrain_count_probs_new_1$ESST_probability))
  
  full_means_3 <- c(mean(alltrain_count_probs_3$DeltaDeltaG),mean(alltrain_count_probs_3$total_score),mean(alltrain_count_probs_3$LogLikelihood),mean(alltrain_count_probs_3$ESST_probability))
  new_means_3 <- c(mean(alltrain_count_probs_new_3$DeltaDeltaG),mean(alltrain_count_probs_new_3$total_score),mean(alltrain_count_probs_new_3$LogLikelihood),mean(alltrain_count_probs_new_3$ESST_probability))
  
  full_means_6 <- c(mean(alltrain_count_probs_6$DeltaDeltaG),mean(alltrain_count_probs_6$total_score),mean(alltrain_count_probs_6$LogLikelihood),mean(alltrain_count_probs_6$ESST_probability))
  new_means_6 <- c(mean(alltrain_count_probs_new_6$DeltaDeltaG),mean(alltrain_count_probs_new_6$total_score),mean(alltrain_count_probs_new_6$LogLikelihood),mean(alltrain_count_probs_new_6$ESST_probability))
  
  full_means_12 <- c(mean(alltrain_count_probs_12$DeltaDeltaG),mean(alltrain_count_probs_12$total_score),mean(alltrain_count_probs_12$LogLikelihood),mean(alltrain_count_probs_12$ESST_probability))
  new_means_12 <- c(mean(alltrain_count_probs_new_12$DeltaDeltaG),mean(alltrain_count_probs_new_12$total_score),mean(alltrain_count_probs_new_12$LogLikelihood),mean(alltrain_count_probs_new_12$ESST_probability))
  
  
  full_out <- rbind(alltrain_count_probs_new_1,alltrain_count_probs_new_3,alltrain_count_probs_new_6,alltrain_count_probs_new_12)
  
  full_out$Date <- date
  
  rm(all_train)
  
  # setwd to print results
  setwd("/Users/user/Documents/Paper\ planning/results+figures/class_switch_data")
  
  write_xlsx(full_out,paste0(structure,"_all_class_switching_subs_",date,".xlsx"))
  
  ### uncomment the code block below if you want to output PDBs to visualise substitutions on the 3D structure ###
  
  # # writing pdbs ------------------------------------------------------------
  # # incorrectly classified substitutions
  # # incorrect 1
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_incorrect_1 <- left_join(pdb_frame, incorrect_resno_counts_1, by = c("chain","resno"))
  # pdb_frame_incorrect_1$b.x <- pdb_frame_incorrect_1$b.y
  # # drop extra column
  # pdb_frame_incorrect_1 <- select(pdb_frame_incorrect_1, -c("b.y"))
  # 
  # pdb_frame_incorrect_1 <- dplyr::rename(pdb_frame_incorrect_1, "b" = "b.x")
  # 
  # # replace NAs(no new/incorrect classification counts) with 0
  # pdb_frame_incorrect_1$b[is.na(pdb_frame_incorrect_1$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_incorrect_1
  # write.pdb(pdb_raw, paste0(structure,"_incorrect1_preds_count.pdb"))
  # 
  # # if you want to look at specific residues with many mistakes...
  # # top_incorrect_1 <- filter(pdb_frame_incorrect_1,b==10|b==9|b==8)
  # 
  # # incorrect 3
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_incorrect_3 <- left_join(pdb_frame, incorrect_resno_counts_3, by = c("chain","resno"))
  # pdb_frame_incorrect_3$b.x <- pdb_frame_incorrect_3$b.y
  # # drop extra column
  # pdb_frame_incorrect_3 <- select(pdb_frame_incorrect_3, -c("b.y"))
  # 
  # pdb_frame_incorrect_3 <- dplyr::rename(pdb_frame_incorrect_3, "b" = "b.x")
  # 
  # # replace NAs(no new/incorrect classification counts) with 0
  # pdb_frame_incorrect_3$b[is.na(pdb_frame_incorrect_3$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_incorrect_3
  # write.pdb(pdb_raw, paste0(structure,"_incorrect3_preds_count.pdb"))
  # 
  # 
  # # incorrect 6
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_incorrect_6 <- left_join(pdb_frame, incorrect_resno_counts_6, by = c("chain","resno"))
  # pdb_frame_incorrect_6$b.x <- pdb_frame_incorrect_6$b.y
  # # drop extra column
  # pdb_frame_incorrect_6 <- select(pdb_frame_incorrect_6, -c("b.y"))
  # 
  # pdb_frame_incorrect_6 <- dplyr::rename(pdb_frame_incorrect_6, "b" = "b.x")
  # 
  # # replace NAs(no new/incorrect classification counts) with 0
  # pdb_frame_incorrect_6$b[is.na(pdb_frame_incorrect_6$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_incorrect_6
  # write.pdb(pdb_raw, paste0(structure,"_incorrect6_preds_count.pdb"))
  # 
  # 
  # # incorrect 12
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_incorrect_12 <- left_join(pdb_frame, incorrect_resno_counts_12, by = c("chain","resno"))
  # pdb_frame_incorrect_12$b.x <- pdb_frame_incorrect_12$b.y
  # # drop extra column
  # pdb_frame_incorrect_12 <- select(pdb_frame_incorrect_12, -c("b.y"))
  # 
  # pdb_frame_incorrect_12 <- dplyr::rename(pdb_frame_incorrect_12, "b" = "b.x")
  # 
  # # replace NAs(no new/incorrect classification counts) with 0
  # pdb_frame_incorrect_12$b[is.na(pdb_frame_incorrect_12$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_incorrect_12
  # write.pdb(pdb_raw, paste0(structure,"_incorrect12_preds_count.pdb"))
  # 
  # 
  # # newly classified substitutions
  # # new1
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_new_1 <- left_join(pdb_frame, new_resno_counts_1, by = c("chain","resno"))
  # pdb_frame_new_1$b.x <- pdb_frame_new_1$b.y
  # # drop extra column
  # pdb_frame_new_1 <- select(pdb_frame_new_1, -c("b.y"))
  # 
  # pdb_frame_new_1 <- dplyr::rename(pdb_frame_new_1, "b" = "b.x")
  # 
  # # replace NAs(no new/newclassification counts) with 0
  # pdb_frame_new_1$b[is.na(pdb_frame_new_1$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_new_1
  # write.pdb(pdb_raw, paste0(structure,"_new1_preds_count.pdb"))
  # 
  # 
  # # new3
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_new_3 <- left_join(pdb_frame, new_resno_counts_3, by = c("chain","resno"))
  # pdb_frame_new_3$b.x <- pdb_frame_new_3$b.y
  # # drop extra column
  # pdb_frame_new_3 <- select(pdb_frame_new_3, -c("b.y"))
  # 
  # pdb_frame_new_3 <- dplyr::rename(pdb_frame_new_3, "b" = "b.x")
  # 
  # # replace NAs(no new/newclassification counts) with 0
  # pdb_frame_new_3$b[is.na(pdb_frame_new_3$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_new_3
  # write.pdb(pdb_raw, paste0(structure,"_new3_preds_count.pdb"))
  # 
  # 
  # # new6
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_new_6 <- left_join(pdb_frame, new_resno_counts_6, by = c("chain","resno"))
  # pdb_frame_new_6$b.x <- pdb_frame_new_6$b.y
  # # drop extra column
  # pdb_frame_new_6 <- select(pdb_frame_new_6, -c("b.y"))
  # 
  # pdb_frame_new_6 <- dplyr::rename(pdb_frame_new_6, "b" = "b.x")
  # 
  # # replace NAs(no new/newclassification counts) with 0
  # pdb_frame_new_6$b[is.na(pdb_frame_new_6$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_new_6
  # write.pdb(pdb_raw, paste0(structure,"_new6_preds_count.pdb"))
  # 
  # 
  # # new12
  # # read in pdb and replace b factor column with counts
  # pdb_raw <- read.pdb(paste0(structure))
  # pdb_frame <- pdb_raw$atom
  # # bind
  # pdb_frame_new_12 <- left_join(pdb_frame, new_resno_counts_12, by = c("chain","resno"))
  # pdb_frame_new_12$b.x <- pdb_frame_new_12$b.y
  # # drop extra column
  # pdb_frame_new_12 <- select(pdb_frame_new_12, -c("b.y"))
  # 
  # pdb_frame_new_12 <- dplyr::rename(pdb_frame_new_12, "b" = "b.x")
  # 
  # # replace NAs(no new/newclassification counts) with 0
  # pdb_frame_new_12$b[is.na(pdb_frame_new_12$b)] <- 0
  # 
  # # write output
  # pdb_raw$atom <- pdb_frame_new_12
  # write.pdb(pdb_raw, paste0(structure,"_new12_preds_count.pdb"))

  } else {
    # if we don't have the 12 month timepoint we just analyse 1, 3 and 6 month timepoints --------
    # read in test sets -------------------------------------------------------
    # read in 1 month test set
    # read in master table of scores
    master_conservation_1 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month1_",date,".xlsx"))
    
    # calculate quantiles of observed probability to defin classes later
    quantiles <- quantile(master_conservation_1$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_1_50 <- master_conservation_1 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    rm(master_conservation_1)
    
    # read in 3 month test set
    # read in master table of scores
    master_conservation_3 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month3_",date,".xlsx"))
    
    # calculate quantiles of observed probability to defin classes later
    quantiles <- quantile(master_conservation_3$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_3_50 <- master_conservation_3 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    rm(master_conservation_3)
    
    # read in 6 month test set
    # read in master table of scores
    master_conservation_6 <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TEST_month6_",date,".xlsx"))
    
    # calculate quantiles of observed probability to defin classes later
    quantiles <- quantile(master_conservation_6$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_6_50 <- master_conservation_6 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    rm(master_conservation_6)
    
    # filter for individual predict sets
    alltrain_count_probs_1 <- alltrain_count_probs %>% 
      filter(Group=="1month")
    
    alltrain_count_probs_3 <- alltrain_count_probs %>% 
      filter(Group=="3month")
    
    alltrain_count_probs_6 <- alltrain_count_probs %>% 
      filter(Group=="6month")

    
    # 1 month
    # make a new tibble with just replacement column to join by and the observed class of test set
    master_conservation_1_50$Replacement <- paste0(aa321(master_conservation_1_50$WTAA),master_conservation_1_50$Chain,master_conservation_1_50$Resno,master_conservation_1_50$ReplacementAA)
    
    master_conservation_1_50 <- master_conservation_1_50 %>% 
      select(Observed_class,Replacement)
    
    # join to the probability results from the ensemble model
    alltrain_count_probs_1 <- left_join(alltrain_count_probs_1,master_conservation_1_50)
    
    # 3 month
    # make a new tibble with just replacement column to join by and the observed class of test set
    master_conservation_3_50$Replacement <- paste0(aa321(master_conservation_3_50$WTAA),master_conservation_3_50$Chain,master_conservation_3_50$Resno,master_conservation_3_50$ReplacementAA)
    
    master_conservation_3_50 <- master_conservation_3_50 %>% 
      select(Observed_class,Replacement)
    
    # join to the probability results from the ensemble model
    alltrain_count_probs_3 <- left_join(alltrain_count_probs_3,master_conservation_3_50)
    
    # 6 month
    # make a new tibble with just replacement column to join by and the observed class of test set
    master_conservation_6_50$Replacement <- paste0(aa321(master_conservation_6_50$WTAA),master_conservation_6_50$Chain,master_conservation_6_50$Resno,master_conservation_6_50$ReplacementAA)
    
    master_conservation_6_50 <- master_conservation_6_50 %>% 
      select(Observed_class,Replacement)
    
    # join to the probability results from the ensemble model
    alltrain_count_probs_6 <- left_join(alltrain_count_probs_6,master_conservation_6_50)
    
    # split the replacement column into composite parts to inspect character of the substitutions
    alltrain_count_probs_1 <- alltrain_count_probs_1 %>% 
      separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
      select(-X)
    
    alltrain_count_probs_3 <- alltrain_count_probs_3 %>% 
      separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
      select(-X)
    
    alltrain_count_probs_6 <- alltrain_count_probs_6 %>% 
      separate(Replacement,into=c("WTAA","Chain","Resno","ReplacementAA"),sep=c(1,2,-1),remove=F) %>% 
      select(-X)

    
    # 'new' subs --------------------------------------------------------------
    
    # filter for 'new' substitutions
    # read in all train sets
    trainset_reader <- function(count){
      train_set <- read_xlsx(paste0(structure,"_spike_master_scores_ML_TRAIN_loosample_",count,"_",date,".xlsx"))
      
      # calculate quantiles of observed probability to define classes later
      quantiles <- quantile(train_set$Observed_probability)
      
      mean_classifier <- quantiles[3]
      
      # remove all potentially troublesome punctuation from column names
      train_set <- train_set %>% 
        dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)
      
      train_set <- train_set %>% 
        dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)
      
      train_set <- train_set %>% 
        dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)
      
      train_set <- train_set %>% 
        dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)
      
      train_set <- train_set %>% 
        mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted"))
      
      train_set$Replacement <- paste0(aa321(train_set$WTAA),train_set$Chain,train_set$Resno,train_set$ReplacementAA)
      
      train_set <- train_set %>% 
        select(Replacement,Observed_class) %>% 
        arrange(Replacement)
      
      return(train_set)
    }
    
    no_samples <- seq(1:10)
    all_train <- lapply(no_samples,trainset_reader)
    
    # make new tables for only new subs
    alltrain_count_probs_new_1 <- alltrain_count_probs_1
    
    alltrain_count_probs_new_3 <- alltrain_count_probs_3
    
    alltrain_count_probs_new_6 <- alltrain_count_probs_6
    
    alltrain_count_probs_new_1 <- alltrain_count_probs_new_1 %>% 
      arrange(Replacement) %>% 
      filter(Observed_class != all_train[[1]]$Observed_class | 
               Observed_class != all_train[[2]]$Observed_class |
               Observed_class != all_train[[3]]$Observed_class |
               Observed_class != all_train[[4]]$Observed_class |
               Observed_class != all_train[[5]]$Observed_class |
               Observed_class != all_train[[6]]$Observed_class |
               Observed_class != all_train[[7]]$Observed_class |
               Observed_class != all_train[[8]]$Observed_class |
               Observed_class != all_train[[9]]$Observed_class |
               Observed_class != all_train[[10]]$Observed_class)
    
    alltrain_count_probs_new_3 <- alltrain_count_probs_new_3 %>% 
      arrange(Replacement) %>% 
      filter(Observed_class != all_train[[1]]$Observed_class | 
               Observed_class != all_train[[2]]$Observed_class |
               Observed_class != all_train[[3]]$Observed_class |
               Observed_class != all_train[[4]]$Observed_class |
               Observed_class != all_train[[5]]$Observed_class |
               Observed_class != all_train[[6]]$Observed_class |
               Observed_class != all_train[[7]]$Observed_class |
               Observed_class != all_train[[8]]$Observed_class |
               Observed_class != all_train[[9]]$Observed_class |
               Observed_class != all_train[[10]]$Observed_class)
    
    alltrain_count_probs_new_6 <- alltrain_count_probs_new_6 %>% 
      arrange(Replacement) %>% 
      filter(Observed_class != all_train[[1]]$Observed_class | 
               Observed_class != all_train[[2]]$Observed_class |
               Observed_class != all_train[[3]]$Observed_class |
               Observed_class != all_train[[4]]$Observed_class |
               Observed_class != all_train[[5]]$Observed_class |
               Observed_class != all_train[[6]]$Observed_class |
               Observed_class != all_train[[7]]$Observed_class |
               Observed_class != all_train[[8]]$Observed_class |
               Observed_class != all_train[[9]]$Observed_class |
               Observed_class != all_train[[10]]$Observed_class)
    
    
    full_out <- rbind(alltrain_count_probs_new_1,alltrain_count_probs_new_3,alltrain_count_probs_new_6)
    
    full_out$Date <- date

    rm(all_train)
    
    # setwd to print results
    setwd("/Users/user/Documents/Paper\ planning/results+figures/class_switch_data")
    
    write_xlsx(full_out,paste0(structure,"_all_class_switching_subs_",date,".xlsx"))
    
    # # writing pdbs ------------------------------------------------------------
    # # incorrectly classified substitutions
    # # incorrect 1
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_incorrect_1 <- left_join(pdb_frame, incorrect_resno_counts_1, by = c("chain","resno"))
    # pdb_frame_incorrect_1$b.x <- pdb_frame_incorrect_1$b.y
    # # drop extra column
    # pdb_frame_incorrect_1 <- select(pdb_frame_incorrect_1, -c("b.y"))
    # 
    # pdb_frame_incorrect_1 <- dplyr::rename(pdb_frame_incorrect_1, "b" = "b.x")
    # 
    # # replace NAs(no new/incorrect classification counts) with 0
    # pdb_frame_incorrect_1$b[is.na(pdb_frame_incorrect_1$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_incorrect_1
    # write.pdb(pdb_raw, paste0(structure,"_incorrect1_preds_count.pdb"))
    # 
    # # if you want to look at specific residues with many mistakes...
    # # top_incorrect_1 <- filter(pdb_frame_incorrect_1,b==10|b==9|b==8)
    # 
    # # incorrect 3
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_incorrect_3 <- left_join(pdb_frame, incorrect_resno_counts_3, by = c("chain","resno"))
    # pdb_frame_incorrect_3$b.x <- pdb_frame_incorrect_3$b.y
    # # drop extra column
    # pdb_frame_incorrect_3 <- select(pdb_frame_incorrect_3, -c("b.y"))
    # 
    # pdb_frame_incorrect_3 <- dplyr::rename(pdb_frame_incorrect_3, "b" = "b.x")
    # 
    # # replace NAs(no new/incorrect classification counts) with 0
    # pdb_frame_incorrect_3$b[is.na(pdb_frame_incorrect_3$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_incorrect_3
    # write.pdb(pdb_raw, paste0(structure,"_incorrect3_preds_count.pdb"))
    # 
    # 
    # # incorrect 6
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_incorrect_6 <- left_join(pdb_frame, incorrect_resno_counts_6, by = c("chain","resno"))
    # pdb_frame_incorrect_6$b.x <- pdb_frame_incorrect_6$b.y
    # # drop extra column
    # pdb_frame_incorrect_6 <- select(pdb_frame_incorrect_6, -c("b.y"))
    # 
    # pdb_frame_incorrect_6 <- dplyr::rename(pdb_frame_incorrect_6, "b" = "b.x")
    # 
    # # replace NAs(no new/incorrect classification counts) with 0
    # pdb_frame_incorrect_6$b[is.na(pdb_frame_incorrect_6$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_incorrect_6
    # write.pdb(pdb_raw, paste0(structure,"_incorrect6_preds_count.pdb"))
    # 
    # 
    # 
    # # newly classified substitutions
    # # new1
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_new_1 <- left_join(pdb_frame, new_resno_counts_1, by = c("chain","resno"))
    # pdb_frame_new_1$b.x <- pdb_frame_new_1$b.y
    # # drop extra column
    # pdb_frame_new_1 <- select(pdb_frame_new_1, -c("b.y"))
    # 
    # pdb_frame_new_1 <- dplyr::rename(pdb_frame_new_1, "b" = "b.x")
    # 
    # # replace NAs(no new/newclassification counts) with 0
    # pdb_frame_new_1$b[is.na(pdb_frame_new_1$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_new_1
    # write.pdb(pdb_raw, paste0(structure,"_new1_preds_count.pdb"))
    # 
    # 
    # # new3
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_new_3 <- left_join(pdb_frame, new_resno_counts_3, by = c("chain","resno"))
    # pdb_frame_new_3$b.x <- pdb_frame_new_3$b.y
    # # drop extra column
    # pdb_frame_new_3 <- select(pdb_frame_new_3, -c("b.y"))
    # 
    # pdb_frame_new_3 <- dplyr::rename(pdb_frame_new_3, "b" = "b.x")
    # 
    # # replace NAs(no new/newclassification counts) with 0
    # pdb_frame_new_3$b[is.na(pdb_frame_new_3$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_new_3
    # write.pdb(pdb_raw, paste0(structure,"_new3_preds_count.pdb"))
    # 
    # 
    # # new6
    # # read in pdb and replace b factor column with counts
    # pdb_raw <- read.pdb(paste0(structure))
    # pdb_frame <- pdb_raw$atom
    # # bind
    # pdb_frame_new_6 <- left_join(pdb_frame, new_resno_counts_6, by = c("chain","resno"))
    # pdb_frame_new_6$b.x <- pdb_frame_new_6$b.y
    # # drop extra column
    # pdb_frame_new_6 <- select(pdb_frame_new_6, -c("b.y"))
    # 
    # pdb_frame_new_6 <- dplyr::rename(pdb_frame_new_6, "b" = "b.x")
    # 
    # # replace NAs(no new/newclassification counts) with 0
    # pdb_frame_new_6$b[is.na(pdb_frame_new_6$b)] <- 0
    # 
    # # write output
    # pdb_raw$atom <- pdb_frame_new_6
    # write.pdb(pdb_raw, paste0(structure,"_new6_preds_count.pdb"))
    # 
  }
  }, error = function(e) {
    NULL
  })
}

# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){
  # expand a grid of all argument combinations
  l <- expand.grid(..., stringsAsFactors=FALSE)
  
  # apply the function
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)
}

structure_list <- c("6vxx","7lws","7v7n","7wp9")

date_list <- c("DEC_20","JAN_21","FEB_21","MAR_21","APR_21","MAY_21","JUN_21","JUL_21","AUG_21","SEP_21","OCT_21","NOV_21","DEC_21",
               "JAN_22","FEB_22","MAR_22","APR_22","MAY_22","JUN_22","JUL_22","AUG_22","SEP_22","OCT_22","NOV_22","DEC_22")

# note we expect this to throw errors for dates that do not use all 4 structures
# as long as structures are in chronological order this is ok
cmapply(analyse_uncertain_new_subs,structure_list,date_list)
