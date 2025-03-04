require(readxl)
require(tidyverse)
require(Biostrings)
require(bio3d)
require(pracma)
require(caret)
require(caTools)
require(writexl)
require(MLeval)
require(doParallel)
require(Matrix)

# read in the correct variant structure and date from directory name
cwd <- getwd()
variant <- str_sub(cwd,start=-4,end=nchar(cwd))
date <- str_split(cwd,"\\/")
date <- date[[1]][8]
date <- str_sub(date,start=-6,end=nchar(date))

# get number of cores from environment variable set in jobscript
numCoresAllowed <- as.numeric(Sys.getenv("MC_CORES", unset=1))
#
ML_train <- function(train_data){
# read in training set
# read in master table of scores
master_conservation <- read_xlsx(train_data)

# calculate quantiles of observed probability to define classes later
quantiles <- quantile(master_conservation$Observed_probability)

mean_classifier <- quantiles[3]

# remove all potentially troublesome punctuation from column names
master_conservation <- master_conservation %>% 
  dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)

master_conservation <- master_conservation %>% 
  dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)

master_conservation <- master_conservation %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)

master_conservation <- master_conservation %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)


master_conservation_50 <- master_conservation %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

# convert all character columns to factors
master_conservation_50[sapply(master_conservation_50, is.character)] <- lapply(master_conservation_50[sapply(master_conservation_50, is.character)], 
                                                                                         as.factor)
# create a vector of outcomes
outcomes_50 <- master_conservation_50$Observed_class

# read in 1 month test set
# read in master table of scores
master_conservation_1 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month1_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_1$Observed_probability)

mean_classifier <- quantiles[3]

# remove all potentially troublesome punctuation from column names
master_conservation_1 <- master_conservation_1 %>% 
  dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)

master_conservation_1 <- master_conservation_1 %>% 
  dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)

master_conservation_1 <- master_conservation_1 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)

master_conservation_1 <- master_conservation_1 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)

master_conservation_1_50 <- master_conservation_1 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_1)

# add train set obs column to test master tables
master_conservation_1_50$train_obs <- master_conservation_50$Observed_class

outcomes_1_50 <- master_conservation_1_50$Observed_class

# read in 3 month test set
# read in master table of scores
master_conservation_3 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month3_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_3$Observed_probability)

mean_classifier <- quantiles[3]

# remove all potentially troublesome punctuation from column names
master_conservation_3 <- master_conservation_3 %>% 
  dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)

master_conservation_3 <- master_conservation_3 %>% 
  dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)

master_conservation_3 <- master_conservation_3 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)

master_conservation_3 <- master_conservation_3 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)

master_conservation_3_50 <- master_conservation_3 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_3)

# add train set obs column to test master tables
master_conservation_3_50$train_obs <- master_conservation_50$Observed_class

outcomes_3_50 <- master_conservation_3_50$Observed_class

# read in 6 month test set
# read in master table of scores
master_conservation_6 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month6_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_6$Observed_probability)

mean_classifier <- quantiles[3]

# remove all potentially troublesome punctuation from column names
master_conservation_6 <- master_conservation_6 %>% 
  dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)

master_conservation_6 <- master_conservation_6 %>% 
  dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)

master_conservation_6 <- master_conservation_6 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)

master_conservation_6 <- master_conservation_6 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)

master_conservation_6_50 <- master_conservation_6 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_6)

# add train set obs column to test master tables
master_conservation_6_50$train_obs <- master_conservation_50$Observed_class

outcomes_6_50 <- master_conservation_6_50$Observed_class

# read in 12 month test set
# read in master table of scores
master_conservation_12 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month12_",date,".xlsx"))

# calculate quantiles of observed probability to define classes later
quantiles <- quantile(master_conservation_12$Observed_probability)

mean_classifier <- quantiles[3]

# remove all potentially troublesome punctuation from column names
master_conservation_12 <- master_conservation_12 %>% 
  dplyr::rename("hydrogen_bond_to_other_sidechain_heterogen" = `hydrogen_bond_to_other_sidechain/heterogen`)

master_conservation_12 <- master_conservation_12 %>% 
  dplyr::rename("cis_peptide_bond" = `cis-peptide_bond`)

master_conservation_12 <- master_conservation_12 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_amide" = `mainchain_to_mainchain_hydrogen_bonds_(amide)`)

master_conservation_12 <- master_conservation_12 %>% 
  dplyr::rename("mainchain_to_mainchain_hydrogen_bonds_carbonyl" = `Mainchain_to_mainchain_hydrogen_bonds_(carbonyl)`)

master_conservation_12_50 <- master_conservation_12 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_12)

# add train set obs column to test master tables
master_conservation_12_50$train_obs <- master_conservation_50$Observed_class

outcomes_12_50 <- master_conservation_12_50$Observed_class

# flag all 'new' substitutions for each test set and store identifying names - these will be used later during analysis
# filter for only changed classes
# add train set obs column to test master tables
master_conservation_1_new_50 <- master_conservation_1_50
master_conservation_3_new_50 <- master_conservation_3_50
master_conservation_6_new_50 <- master_conservation_6_50
master_conservation_12_new_50 <- master_conservation_12_50

master_conservation_1_new_50$train_obs <- master_conservation_50$Observed_class
master_conservation_3_new_50$train_obs <- master_conservation_50$Observed_class
master_conservation_6_new_50$train_obs <- master_conservation_50$Observed_class
master_conservation_12_new_50$train_obs <- master_conservation_50$Observed_class

master_conservation_1_new_50 <- master_conservation_1_50 %>% 
  filter(Observed_class != train_obs)

master_conservation_3_new_50 <- master_conservation_3_50 %>% 
  filter(Observed_class != train_obs)

master_conservation_6_new_50 <- master_conservation_6_50 %>% 
  filter(Observed_class != train_obs)

master_conservation_12_new_50 <- master_conservation_12_50 %>% 
  filter(Observed_class != train_obs)

# make a unique naming vector for every mutation (these are the same for all sets so any can be used)
new_1_index <- paste0(aa321(master_conservation_1_new_50$WTAA),master_conservation_1_new_50$Chain,master_conservation_1_new_50$Resno,master_conservation_1_new_50$ReplacementAA)
new_3_index <- paste0(aa321(master_conservation_3_new_50$WTAA),master_conservation_3_new_50$Chain,master_conservation_3_new_50$Resno,master_conservation_3_new_50$ReplacementAA)
new_6_index <- paste0(aa321(master_conservation_6_new_50$WTAA),master_conservation_6_new_50$Chain,master_conservation_6_new_50$Resno,master_conservation_6_new_50$ReplacementAA)
new_12_index <- paste0(aa321(master_conservation_12_new_50$WTAA),master_conservation_12_new_50$Chain,master_conservation_12_new_50$Resno,master_conservation_12_new_50$ReplacementAA)

# remove the 'new' sub tibbles
rm(master_conservation_1_new_50,master_conservation_3_new_50,master_conservation_6_new_50,master_conservation_12_new_50)

sample_id <- str_sub(train_data,start=-23,end=-13)

# write the indices of new subs
write_csv(as_tibble(new_1_index),paste0("new_subs_1_index_",sample_id,".csv"),col_names = F)
write_csv(as_tibble(new_3_index),paste0("new_subs_3_index_",sample_id,".csv"),col_names = F)
write_csv(as_tibble(new_6_index),paste0("new_subs_6_index_",sample_id,".csv"),col_names = F)
write_csv(as_tibble(new_12_index),paste0("new_subs_12_index_",sample_id,".csv"),col_names = F)

# make a unique naming vector for every mutation (these are the same for all sets so any can be used)
test_names <- paste0(aa321(master_conservation_6_50$WTAA),master_conservation_6_50$Chain,master_conservation_6_50$Resno,master_conservation_6_50$ReplacementAA)

# we have raw train and test tables, now we need some additional processing
# filter for only predictors selected by feature selection
master_conservation_50 <- master_conservation_50 %>% 
  select(ReplacementAA,ESST_probability,ref,fa_dun,entropy_sidechain,fa_sol,LogLikelihood,
           pred_epitope_rank,fa_atr,entropy_mainchain,Solvation_Polar,backbone_clash,hbond_sc,Solvation_Hydrophobic,
           Van_der_Waals,fa_rep,fa_pair,p_aa_pp,Backbone_Hbond,total_score,torsional_clash,WTAA,hbond_bb_sc,Electrostatics,
           Ooi_number,Sidechain_Hbond,Van_der_Waals_clashes,energy_Ionisation,dslf_ss_dih,electrostatic_kon)

master_conservation_1_50 <- master_conservation_1_50 %>% 
  select(ReplacementAA,ESST_probability,ref,fa_dun,entropy_sidechain,fa_sol,LogLikelihood,
         pred_epitope_rank,fa_atr,entropy_mainchain,Solvation_Polar,backbone_clash,hbond_sc,Solvation_Hydrophobic,
         Van_der_Waals,fa_rep,fa_pair,p_aa_pp,Backbone_Hbond,total_score,torsional_clash,WTAA,hbond_bb_sc,Electrostatics,
         Ooi_number,Sidechain_Hbond,Van_der_Waals_clashes,energy_Ionisation,dslf_ss_dih,electrostatic_kon)


master_conservation_3_50 <- master_conservation_3_50 %>% 
  select(ReplacementAA,ESST_probability,ref,fa_dun,entropy_sidechain,fa_sol,LogLikelihood,
         pred_epitope_rank,fa_atr,entropy_mainchain,Solvation_Polar,backbone_clash,hbond_sc,Solvation_Hydrophobic,
         Van_der_Waals,fa_rep,fa_pair,p_aa_pp,Backbone_Hbond,total_score,torsional_clash,WTAA,hbond_bb_sc,Electrostatics,
         Ooi_number,Sidechain_Hbond,Van_der_Waals_clashes,energy_Ionisation,dslf_ss_dih,electrostatic_kon)


master_conservation_6_50 <- master_conservation_6_50 %>% 
  select(ReplacementAA,ESST_probability,ref,fa_dun,entropy_sidechain,fa_sol,LogLikelihood,
         pred_epitope_rank,fa_atr,entropy_mainchain,Solvation_Polar,backbone_clash,hbond_sc,Solvation_Hydrophobic,
         Van_der_Waals,fa_rep,fa_pair,p_aa_pp,Backbone_Hbond,total_score,torsional_clash,WTAA,hbond_bb_sc,Electrostatics,
         Ooi_number,Sidechain_Hbond,Van_der_Waals_clashes,energy_Ionisation,dslf_ss_dih,electrostatic_kon)


master_conservation_12_50 <- master_conservation_12_50 %>% 
  select(ReplacementAA,ESST_probability,ref,fa_dun,entropy_sidechain,fa_sol,LogLikelihood,
         pred_epitope_rank,fa_atr,entropy_mainchain,Solvation_Polar,backbone_clash,hbond_sc,Solvation_Hydrophobic,
         Van_der_Waals,fa_rep,fa_pair,p_aa_pp,Backbone_Hbond,total_score,torsional_clash,WTAA,hbond_bb_sc,Electrostatics,
         Ooi_number,Sidechain_Hbond,Van_der_Waals_clashes,energy_Ionisation,dslf_ss_dih,electrostatic_kon)


# set seed for consistency
set.seed(4681)

# create cross validation folds
cvFolds_5 <- createMultiFolds(outcomes_50, k=10,times=25)

# define train control with precision/recall or ROC summary function
fitControl_cv <- trainControl(
  method="repeatedcv",
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = TRUE,
  savePredictions = TRUE,
  returnResamp = "all",
  index = cvFolds_5
)
# RF -----------------------------------------------------------

# rf model with predictors selected by RFE (final set selected through rank aggregation)
model_rf <- train(x=master_conservation_50,y=outcomes_50,
                  metric="ROC",
                  method="ranger",
                  tuneGrid = expand.grid(
                    .mtry = c(6),
                    .splitrule = "extratrees",
                    .min.node.size = c(1)
                  ),
                  trControl=fitControl_cv,
                  num.threads=1
)

name <- str_sub(train_data,29,nchar(train_data)-5)

write.csv(model_rf$results,paste0(name,"_rf_results.csv"))

# predict -----------------------------------------------------------------
# predict on test sets
probs_preds_1 <- predict(model_rf,master_conservation_1_50,type="prob")  
probs_preds_1$Prediction <- predict(model_rf,master_conservation_1_50)                                                                               
probs_preds_1$Group <- "1month"
probs_preds_1$Replacement <- test_names

probs_preds_3 <- predict(model_rf,master_conservation_3_50,type="prob")  
probs_preds_3$Prediction <- predict(model_rf,master_conservation_3_50)                                                                               
probs_preds_3$Group <- "3month"
probs_preds_3$Replacement <- test_names

probs_preds_6 <- predict(model_rf,master_conservation_6_50,type="prob")  
probs_preds_6$Prediction <- predict(model_rf,master_conservation_6_50)                                                                               
probs_preds_6$Group <- "6month"
probs_preds_6$Replacement <- test_names

probs_preds_12 <- predict(model_rf,master_conservation_12_50,type="prob") 
probs_preds_12$Prediction <- predict(model_rf,master_conservation_12_50)                                                                               
probs_preds_12$Group <- "12month"
probs_preds_12$Replacement <- test_names

out_preds <- rbind(probs_preds_1,probs_preds_3,probs_preds_6,probs_preds_12)

return(out_preds)

}

train_sets <- list.files(pattern="*_master_scores_ML_TRAIN_*")

# open multithread cluster
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
# train models
alltrain_out <- lapply(train_sets,ML_train)
# close cluster
stopCluster(cl)

# # uncomment this code if you want results for individual sets printed
# for(i in 1:10){
# write.csv(alltrain_out[i],paste0("rf_50_ensemble_rawpreds_set",i,".csv"),row.names = F)
# }

alltrain_full <- bind_rows(alltrain_out, .id = "train_set_id")

# count the number of each prediction for each replacement
alltrain_count <- alltrain_full %>% 
  group_by(Group,Replacement) %>% 
  count(Prediction,.drop=F)

# calculate the percentage of these probabilities
alltrain_count_probs <- alltrain_count %>% 
  mutate(Ensemble_vote = case_when(Prediction=="Permitted" ~ n/10,Prediction=="Forbidden" ~ 0+n/10)) 

alltrain_count_probs <- alltrain_count_probs %>% 
  select(-n) %>% 
  pivot_wider(names_from = Prediction,values_from=Ensemble_vote)


write.csv(alltrain_count,"rf_50_ensemble_counts.csv")
write.csv(alltrain_count_probs,"rf_50_ensemble_count_probs.csv")


# Analysis ----------------------------------------------------------------

# analyse the results of the ensemble model
# read in the test sets to assess aggregated results
# read in 1 month test set
# read in master table of scores
master_conservation_1 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month1_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_1$Observed_probability)

mean_classifier <- quantiles[3]

master_conservation_1_50 <- master_conservation_1 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_1)

outcomes_1_50 <- master_conservation_1_50$Observed_class

# read in 3 month test set
# read in master table of scores
master_conservation_3 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month3_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_3$Observed_probability)

mean_classifier <- quantiles[3]

master_conservation_3_50 <- master_conservation_3 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_3)

outcomes_3_50 <- master_conservation_3_50$Observed_class


# read in 6 month test set
# read in master table of scores
master_conservation_6 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month6_",date,".xlsx"))

# calculate quantiles of observed probability to defin classes later
quantiles <- quantile(master_conservation_6$Observed_probability)

mean_classifier <- quantiles[3]

master_conservation_6_50 <- master_conservation_6 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_6)

outcomes_6_50 <- master_conservation_6_50$Observed_class

# read in 12 month test set
# read in master table of scores
master_conservation_12 <- read_xlsx(paste0(variant,"_spike_master_scores_ML_TEST_month12_",date,".xlsx"))

# calculate quantiles of observed probability to define classes later
quantiles <- quantile(master_conservation_12$Observed_probability)

mean_classifier <- quantiles[3]

master_conservation_12_50 <- master_conservation_12 %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

rm(master_conservation_12)

outcomes_12_50 <- master_conservation_12_50$Observed_class


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

master_conservation_1_50 <- master_conservation_1_50 %>% 
  select(Observed_class,Replacement)

# join to the probability results from the ensemble model
alltrain_count_probs_full_1 <- left_join(alltrain_count_probs_1,master_conservation_1_50)

# 3 month
# make a new tibble with just replacement column to join by and the observed class of test set
master_conservation_3_50$Replacement <- paste0(aa321(master_conservation_3_50$WTAA),master_conservation_3_50$Chain,master_conservation_3_50$Resno,master_conservation_3_50$ReplacementAA)

master_conservation_3_50 <- master_conservation_3_50 %>% 
  select(Observed_class,Replacement)

# join to the probability results from the ensemble model
alltrain_count_probs_full_3 <- left_join(alltrain_count_probs_3,master_conservation_3_50)

# 6 month
# make a new tibble with just replacement column to join by and the observed class of test set
master_conservation_6_50$Replacement <- paste0(aa321(master_conservation_6_50$WTAA),master_conservation_6_50$Chain,master_conservation_6_50$Resno,master_conservation_6_50$ReplacementAA)

master_conservation_6_50 <- master_conservation_6_50 %>% 
  select(Observed_class,Replacement)

# join to the probability results from the ensemble model
alltrain_count_probs_full_6 <- left_join(alltrain_count_probs_6,master_conservation_6_50)

# 12 month
# make a new tibble with just replacement column to join by and the observed class of test set
master_conservation_12_50$Replacement <- paste0(aa321(master_conservation_12_50$WTAA),master_conservation_12_50$Chain,master_conservation_12_50$Resno,master_conservation_12_50$ReplacementAA)

master_conservation_12_50 <- master_conservation_12_50 %>% 
  select(Observed_class,Replacement)

# join to the probability results from the ensemble model
alltrain_count_probs_full_12 <- left_join(alltrain_count_probs_12,master_conservation_12_50)

# set up tibble for PR-gain curve calculation with Mleval
prg_tbl_1 <- tibble("Permitted"=alltrain_count_probs_full_1$Permitted,"Forbidden"=alltrain_count_probs_full_1$Forbidden,"obs"=alltrain_count_probs_full_1$Observed_class,"Group"=alltrain_count_probs_full_1$Group)
prg_tbl_3 <- tibble("Permitted"=alltrain_count_probs_full_3$Permitted,"Forbidden"=alltrain_count_probs_full_3$Forbidden,"obs"=alltrain_count_probs_full_3$Observed_class,"Group"=alltrain_count_probs_full_3$Group)
prg_tbl_6 <- tibble("Permitted"=alltrain_count_probs_full_6$Permitted,"Forbidden"=alltrain_count_probs_full_6$Forbidden,"obs"=alltrain_count_probs_full_6$Observed_class,"Group"=alltrain_count_probs_full_6$Group)
prg_tbl_12 <- tibble("Permitted"=alltrain_count_probs_full_12$Permitted,"Forbidden"=alltrain_count_probs_full_12$Forbidden,"obs"=alltrain_count_probs_full_12$Observed_class,"Group"=alltrain_count_probs_full_12$Group)

prg_tbl <- rbind(prg_tbl_1,prg_tbl_3,prg_tbl_6,prg_tbl_12)

prg_tbl$Group <- factor(prg_tbl$Group,levels=c("1month","3month","6month","12month"))

full_eval <- evalm(data.frame(prg_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12)


# assess performance on class switching substitutions
# read in indices of new substitution
# any substitution that changes class between at least 1 training samples and the test sample we consider to have switched class

# 1 month
newsub_filelist <- list.files(pattern="new_subs_1_index*")

newsub_all <- lapply(newsub_filelist,read_csv,col_names=F)

# using one of the other replacement columns, filter for subs that meet our criteria

newsub_filtered <- list()

for(i in 1:length(newsub_all)){
newsub_filtered[[i]] <- master_conservation_1_50 %>% 
  filter(Replacement %in% newsub_all[[i]]$X1)

newsub_filtered_tbl <- bind_rows(newsub_filtered)

newsub_filtered_tbl_1 <- unique(newsub_filtered_tbl)
}


# 3 month
newsub_filelist <- list.files(pattern="new_subs_3_index*")

newsub_all <- lapply(newsub_filelist,read_csv,col_names=F)

# using one of the other replacement columns, filter for subs that meet our criteria

newsub_filtered <- list()

for(i in 1:length(newsub_all)){
  newsub_filtered[[i]] <- master_conservation_3_50 %>% 
    filter(Replacement %in% newsub_all[[i]]$X1)
  
  newsub_filtered_tbl <- bind_rows(newsub_filtered)
  
  newsub_filtered_tbl_3 <- unique(newsub_filtered_tbl)
}


# 6 month
newsub_filelist <- list.files(pattern="new_subs_6_index*")

newsub_all <- lapply(newsub_filelist,read_csv,col_names=F)

# using one of the other replacement columns, filter for subs that meet our criteria

newsub_filtered <- list()

for(i in 1:length(newsub_all)){
  newsub_filtered[[i]] <- master_conservation_6_50 %>% 
    filter(Replacement %in% newsub_all[[i]]$X1)
  
  newsub_filtered_tbl <- bind_rows(newsub_filtered)
  
  newsub_filtered_tbl_6 <- unique(newsub_filtered_tbl)
}


# 12 month
newsub_filelist <- list.files(pattern="new_subs_12_index*")

newsub_all <- lapply(newsub_filelist,read_csv,col_names=F)

# using one of the other replacement columns, filter for subs that meet our criteria

newsub_filtered <- list()

for(i in 1:length(newsub_all)){
  newsub_filtered[[i]] <- master_conservation_12_50 %>% 
    filter(Replacement %in% newsub_all[[i]]$X1)
  
  newsub_filtered_tbl <- bind_rows(newsub_filtered)
  
  newsub_filtered_tbl_12 <- unique(newsub_filtered_tbl)
}


# join to the probability results from the ensemble model
alltrain_count_probs_new_1 <- drop_na(left_join(alltrain_count_probs_1,newsub_filtered_tbl_1))

alltrain_count_probs_new_3 <- drop_na(left_join(alltrain_count_probs_3,newsub_filtered_tbl_3))

alltrain_count_probs_new_6 <- drop_na(left_join(alltrain_count_probs_6,newsub_filtered_tbl_6))

alltrain_count_probs_new_12 <- drop_na(left_join(alltrain_count_probs_12,newsub_filtered_tbl_12))


# set up tibble for PR-gain curve calculation with Mleval
prg_tbl_1 <- tibble("Permitted"=alltrain_count_probs_new_1$Permitted,"Forbidden"=alltrain_count_probs_new_1$Forbidden,"obs"=alltrain_count_probs_new_1$Observed_class,"Group"=alltrain_count_probs_new_1$Group)
prg_tbl_3 <- tibble("Permitted"=alltrain_count_probs_new_3$Permitted,"Forbidden"=alltrain_count_probs_new_3$Forbidden,"obs"=alltrain_count_probs_new_3$Observed_class,"Group"=alltrain_count_probs_new_3$Group)
prg_tbl_6 <- tibble("Permitted"=alltrain_count_probs_new_6$Permitted,"Forbidden"=alltrain_count_probs_new_6$Forbidden,"obs"=alltrain_count_probs_new_6$Observed_class,"Group"=alltrain_count_probs_new_6$Group)
prg_tbl_12 <- tibble("Permitted"=alltrain_count_probs_new_12$Permitted,"Forbidden"=alltrain_count_probs_new_12$Forbidden,"obs"=alltrain_count_probs_new_12$Observed_class,"Group"=alltrain_count_probs_new_12$Group)

prg_tbl <- rbind(prg_tbl_1,prg_tbl_3,prg_tbl_6,prg_tbl_12)

prg_tbl$Group <- factor(prg_tbl$Group,levels=c("1month","3month","6month","12month"))

new_eval <- evalm(data.frame(prg_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12)


# assess the substitutions with modeled uncertainty (that is those that do not have unanimous agreement across the ensemble)
# calculate the percentage for the non-unanimously agreeing substitutions
alltrain_uncertain <- alltrain_count %>% 
  filter(n!=10) %>% 
  filter(n!=0)

# calculate ratios of uncertain vs certain predictions and new vs total permitted substitutions
uncertainty_ratio <- nrow(alltrain_uncertain)/nrow(alltrain_count)

new_ratio_1 <- nrow(alltrain_count_probs_new_1)/nrow(master_conservation_1_50)
new_ratio_3 <- nrow(alltrain_count_probs_new_3)/nrow(master_conservation_3_50)
new_ratio_6 <- nrow(alltrain_count_probs_new_6)/nrow(master_conservation_6_50)
new_ratio_12 <- nrow(alltrain_count_probs_new_12)/nrow(master_conservation_12_50)


alltrain_uncertain_count_probs <- alltrain_uncertain %>% 
  mutate(Ensemble_vote = case_when(Prediction=="Permitted" ~ n/10,Prediction=="Forbidden" ~ 0+n/10)) 

alltrain_uncertain_count_probs <- alltrain_uncertain_count_probs %>% 
  select(-n) %>% 
  pivot_wider(names_from = Prediction,values_from=Ensemble_vote)

alltrain_uncertain_count_probs_1 <- alltrain_uncertain_count_probs %>% 
  filter(Group=="1month")

alltrain_uncertain_count_probs_3 <- alltrain_uncertain_count_probs %>% 
  filter(Group=="3month")

alltrain_uncertain_count_probs_6 <- alltrain_uncertain_count_probs %>% 
  filter(Group=="6month")

alltrain_uncertain_count_probs_12 <- alltrain_uncertain_count_probs %>% 
  filter(Group=="12month")

# join to the probability results from the ensemble model
alltrain_uncertain_count_probs_full_1 <- left_join(alltrain_uncertain_count_probs_1,master_conservation_1_50)
alltrain_uncertain_count_probs_full_3 <- left_join(alltrain_uncertain_count_probs_3,master_conservation_3_50)
alltrain_uncertain_count_probs_full_6 <- left_join(alltrain_uncertain_count_probs_6,master_conservation_6_50)
alltrain_uncertain_count_probs_full_12 <- left_join(alltrain_uncertain_count_probs_12,master_conservation_6_50)

# set up tibble for PR-gain curve calculation with Mleval
prg_tbl_1 <- tibble("Permitted"=alltrain_uncertain_count_probs_full_1$Permitted,"Forbidden"=alltrain_uncertain_count_probs_full_1$Forbidden,"obs"=alltrain_uncertain_count_probs_full_1$Observed_class,"Group"=alltrain_uncertain_count_probs_full_1$Group)
prg_tbl_3 <- tibble("Permitted"=alltrain_uncertain_count_probs_full_3$Permitted,"Forbidden"=alltrain_uncertain_count_probs_full_3$Forbidden,"obs"=alltrain_uncertain_count_probs_full_3$Observed_class,"Group"=alltrain_uncertain_count_probs_full_3$Group)
prg_tbl_6 <- tibble("Permitted"=alltrain_uncertain_count_probs_full_6$Permitted,"Forbidden"=alltrain_uncertain_count_probs_full_6$Forbidden,"obs"=alltrain_uncertain_count_probs_full_6$Observed_class,"Group"=alltrain_uncertain_count_probs_full_6$Group)
prg_tbl_12 <- tibble("Permitted"=alltrain_uncertain_count_probs_full_12$Permitted,"Forbidden"=alltrain_uncertain_count_probs_full_12$Forbidden,"obs"=alltrain_uncertain_count_probs_full_12$Observed_class,"Group"=alltrain_uncertain_count_probs_full_12$Group)

prg_tbl <- rbind(prg_tbl_1,prg_tbl_3,prg_tbl_6,prg_tbl_12)

prg_tbl$Group <- as.factor(prg_tbl$Group)

full_uncertain_eval <- evalm(data.frame(prg_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12)

# assess predictions for uncertain and class switching substitutions
# join to the probability results from the ensemble model
alltrain_uncertain_count_probs_new_1 <- drop_na(left_join(alltrain_uncertain_count_probs_1,alltrain_count_probs_new_1))
alltrain_uncertain_count_probs_new_3 <- drop_na(left_join(alltrain_uncertain_count_probs_3,alltrain_count_probs_new_3))
alltrain_uncertain_count_probs_new_6 <- drop_na(left_join(alltrain_uncertain_count_probs_6,alltrain_count_probs_new_6))
alltrain_uncertain_count_probs_new_12 <- drop_na(left_join(alltrain_uncertain_count_probs_12,alltrain_count_probs_new_12))

# set up tibble for PR-gain curve calculation with Mleval
prg_tbl_1 <- tibble("Permitted"=alltrain_uncertain_count_probs_new_1$Permitted,"Forbidden"=alltrain_uncertain_count_probs_new_1$Forbidden,"obs"=alltrain_uncertain_count_probs_new_1$Observed_class,"Group"=alltrain_uncertain_count_probs_new_1$Group)
prg_tbl_3 <- tibble("Permitted"=alltrain_uncertain_count_probs_new_3$Permitted,"Forbidden"=alltrain_uncertain_count_probs_new_3$Forbidden,"obs"=alltrain_uncertain_count_probs_new_3$Observed_class,"Group"=alltrain_uncertain_count_probs_new_3$Group)
prg_tbl_6 <- tibble("Permitted"=alltrain_uncertain_count_probs_new_6$Permitted,"Forbidden"=alltrain_uncertain_count_probs_new_6$Forbidden,"obs"=alltrain_uncertain_count_probs_new_6$Observed_class,"Group"=alltrain_uncertain_count_probs_new_6$Group)
prg_tbl_12 <- tibble("Permitted"=alltrain_uncertain_count_probs_new_12$Permitted,"Forbidden"=alltrain_uncertain_count_probs_new_12$Forbidden,"obs"=alltrain_uncertain_count_probs_new_12$Observed_class,"Group"=alltrain_uncertain_count_probs_new_12$Group)

prg_tbl <- rbind(prg_tbl_1,prg_tbl_3,prg_tbl_6,prg_tbl_12)

prg_tbl$Group <- factor(prg_tbl$Group,levels=c("1month","3month","6month","12month"))

new_uncertain_eval <- evalm(data.frame(prg_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12)

# calculate ratio of new mutations in the uncertain dataset
new_uncertain_ratio_1 <- nrow(alltrain_uncertain_count_probs_new_1)/nrow(alltrain_uncertain_count_probs_full_1)
new_uncertain_ratio_3 <- nrow(alltrain_uncertain_count_probs_new_3)/nrow(alltrain_uncertain_count_probs_full_3)
new_uncertain_ratio_6 <- nrow(alltrain_uncertain_count_probs_new_6)/nrow(alltrain_uncertain_count_probs_full_6)
new_uncertain_ratio_12 <- nrow(alltrain_uncertain_count_probs_new_12)/nrow(alltrain_uncertain_count_probs_full_12)

# finally, assess performance on 'old' mutations (those that have not changed between aggregated train and test data)

# filter for only changed classes
master_conservation_1_old_50 <- master_conservation_1_50 %>% 
  filter(!(Replacement %in% alltrain_count_probs_new_1$Replacement))

master_conservation_3_old_50 <- master_conservation_3_50 %>% 
  filter(!(Replacement %in% alltrain_count_probs_new_3$Replacement))

master_conservation_6_old_50 <- master_conservation_6_50 %>% 
  filter(!(Replacement %in% alltrain_count_probs_new_6$Replacement))

master_conservation_12_old_50 <- master_conservation_12_50 %>% 
  filter(!(Replacement %in% alltrain_count_probs_new_12$Replacement))

# join to the probability results from the ensemble model
alltrain_count_probs_old_1 <- drop_na(left_join(alltrain_count_probs_1,master_conservation_1_old_50))

alltrain_count_probs_old_3 <- drop_na(left_join(alltrain_count_probs_3,master_conservation_3_old_50))

alltrain_count_probs_old_6 <- drop_na(left_join(alltrain_count_probs_6,master_conservation_6_old_50))

alltrain_count_probs_old_12 <- drop_na(left_join(alltrain_count_probs_12,master_conservation_12_old_50))


# set up tibble for PR-gain curve calculation with Mleval
prg_tbl_1 <- tibble("Permitted"=alltrain_count_probs_old_1$Permitted,"Forbidden"=alltrain_count_probs_old_1$Forbidden,"obs"=alltrain_count_probs_old_1$Observed_class,"Group"=alltrain_count_probs_old_1$Group)
prg_tbl_3 <- tibble("Permitted"=alltrain_count_probs_old_3$Permitted,"Forbidden"=alltrain_count_probs_old_3$Forbidden,"obs"=alltrain_count_probs_old_3$Observed_class,"Group"=alltrain_count_probs_old_3$Group)
prg_tbl_6 <- tibble("Permitted"=alltrain_count_probs_old_6$Permitted,"Forbidden"=alltrain_count_probs_old_6$Forbidden,"obs"=alltrain_count_probs_old_6$Observed_class,"Group"=alltrain_count_probs_old_6$Group)
prg_tbl_12 <- tibble("Permitted"=alltrain_count_probs_old_12$Permitted,"Forbidden"=alltrain_count_probs_old_12$Forbidden,"obs"=alltrain_count_probs_old_12$Observed_class,"Group"=alltrain_count_probs_old_12$Group)

prg_tbl <- rbind(prg_tbl_1,prg_tbl_3,prg_tbl_6,prg_tbl_12)

prg_tbl$Group <- factor(prg_tbl$Group,levels=c("1month","3month","6month","12month"))

old_eval <- evalm(data.frame(prg_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12)

setwd("./std_res")

# write all the outputs of interest
write.csv(full_eval$stdres,"full_eval_ROC_50.csv")
write.csv(new_eval$stdres,"new_eval_ROC_50.csv")
write.csv(full_uncertain_eval$stdres,"full_uncertain_eval_ROC_50.csv")
write.csv(new_uncertain_eval$stdres,"new_uncertain_eval_ROC_50.csv")
write.csv(old_eval$stdres,"old_eval_ROC_50.csv")


write.csv(c(new_ratio_1,new_ratio_3,new_ratio_6,new_ratio_12),"new_sub_ratios.csv")
write.csv(c(new_uncertain_ratio_1,new_uncertain_ratio_3,new_uncertain_ratio_6,new_uncertain_ratio_12),"new_uncertain_sub_ratios.csv")
write.csv(uncertainty_ratio,"uncertainty_ratio.csv")

png(filename="full_eval_ROC_cc.png",type="cairo")
full_eval$cc
dev.off()

png(filename="new_eval_ROC_cc.png",type="cairo")
new_eval$cc
dev.off()

png(filename="full_uncertain_eval_ROC_cc.png",type="cairo")
full_uncertain_eval$cc
dev.off()

png(filename="new_uncertain_eval_ROC_cc.png",type="cairo")
new_uncertain_eval$cc
dev.off()

png(filename="old_eval_ROC_cc.png",type="cairo")
old_eval$cc
dev.off()

