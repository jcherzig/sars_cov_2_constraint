library(tidyverse)
library(bio3d)
library(readxl)
library(writexl)

master_table_writer <- function(structure){
  # MAC paths
  setwd("/Users/user/Documents/SARS-CoV-2")
  
  # Windows paths
  #setwd("C:/Users/jamie/OneDrive - The University of Manchester/SARS-CoV-2/Structure")
  
  # Read PDB file to create list of mutations -------------------------------
  
  # read pdb file
  ## REMEMBER to set the pdb_code object
  raw_pdb <- read.pdb(paste(structure))
  
  pdb_frame <- raw_pdb$atom
  rm(raw_pdb)
  
  pdb_frame <- as_tibble(pdb_frame)
  
  #remove NAG
  pdb_frame <- pdb_frame %>% 
    filter(resid != "NAG")
  
  # filter for Calpha atoms only
  pdb_frame <- pdb_frame %>% 
    filter(elety == "CA")
  
  #list aa replacement to join to later
  aa_list <- c("A", "C", "D", "E", "F", "G", "H", "I",
               "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  ReplacementAA <- rep(aa_list, times=nrow(pdb_frame))
  
  #create empty master tibble to contain all data with nrows = to number of mutations
  master_conservation <- tibble("Resno" = rep(pdb_frame$resno,each=20),
                                "WTAA" = rep(pdb_frame$resid,each=20),
                                "Chain" = rep(pdb_frame$chain,each=20),
                                "ReplacementAA" = ReplacementAA,
                                "BFactor" = rep(pdb_frame$b,each=20),
                                .rows=(20*nrow(pdb_frame)))
  
  # # Read FoldX DeltaDeltaGs from buildmodel (simultaneous chain mutation) --------
  #
  setwd(paste("/Users/user/Documents/SARS-CoV-2/Structure/Foldx/foldx5_buildmodel_all_",structure,sep=""))
  
  foldx_6vxx_buildmodel_master <- read_csv(paste(structure,"_buildmodel_output_master_ML.csv",sep=""))
  
  # replicate each row 3 times
  foldx_6vxx_buildmodel_master <- foldx_6vxx_buildmodel_master %>% dplyr::slice(rep(1:n(), each = 3))
  # add a chain column
  foldx_6vxx_buildmodel_master$Chain <- rep(c("A","B","C"),times=nrow(foldx_6vxx_buildmodel_master)/3)
  
  master_conservation <- right_join(foldx_6vxx_buildmodel_master, master_conservation)
  
  
  # Read ESSTs from Crescendo dumpouts --------------------------------------
  
  #set wd to parent
  setwd(paste("/Users/user/Documents/SARS-CoV-2/Structure/Crescendo/",structure,sep=""))
  
  # read all values from single file
  esst_all <- read_csv(paste("PROCESSED_crescendo_dumpouts_",structure,".csv",sep=""))
  
  colnames(esst_all) <- c("Resno","WTAA","ReplacementAA","ESST_probability","Chain")
  
  # this is necessary because the Blundell lab assigns 'J' to cysteines which do not form disulfide bonds
  esst_all$WTAA[esst_all$WTAA == "J"] <- "C"
  
  esst_all$WTAA <- aa123(esst_all$WTAA)
  
  master_conservation <- right_join(esst_all, master_conservation)
  
  # Read environmental AA information from processed TEM file --------------------------------------
  
  # read all values from single file
  TEM_full <- read_csv(paste0(structure,"_TEM_output.csv"))
  
  TEM_full <- TEM_full %>% dplyr::rename("WTAA" = sequence)
  
  # replicate each row 20 times
  TEM_full <- TEM_full %>% dplyr::slice(rep(1:n(), each = 20))
  TEM_full$ReplacementAA <- master_conservation$ReplacementAA
  
  master_conservation <- right_join(TEM_full, master_conservation)
  
  
  # Read PSSMs from processed CSVs ------------------------------------------
  # PSSMs manually extracted from CDD ---------------------------------------
  # #read processed PSSMs from csv
  # setwd("./S1_NTD")
  # s1_ntd_pssm <- read_csv("PROCESSED_allbeta_cd21527_pssm.csv")
  # setwd("..")
  # 
  # setwd("./S1_RBD")
  # s1_rbd_pssm<- read_csv("PROCESSED_all_beta_cd21470_pssm.csv")
  # setwd("..")
  # 
  # setwd("./S2_and_S1_S2_cleavage")
  # s2_cleavage_pssm<- read_csv("PROCESSED_allbeta_cd22370_pssm.csv")
  # setwd("..")
  
  
  # PSSMs from DELTABLAST ---------------------------------------------------
  
  #set wd to parent
  setwd(paste("/Users/user/Documents/SARS-CoV-2/BLAST/DELTABLAST/",structure,sep=""))
  
  pssm_all <- read_csv(paste("PROCESSED_DELTABLAST_",structure,"_pssm.csv",sep=""))
  
  ### NOTE ####
  # the code commented out below is only there to convert the log likelihood values to a linear likelihood
  # no need for it unless you want these values for some reason
  
  # # scale log-odds values back to 1
  # ### NOTE: this step is only necessary if using PSSMs directly from the CDD database
  # ### if DELTABLAST or PSIBLAST is used values are already at scaling 1
  # # pssm_all <- pssm_all %>% 
  # #   mutate(LogLikelihood = LogLikelihood/100)
  # 
  # # extract lambda score from PSSM file
  # # read PSSM file from DELTA BLAST
  # spike_pssm_raw <- read_file("dblast_6vxx_pssm.out")
  # 
  # # extract ungapped lambda score
  # ### SHOULD WE BE USING UNGAPPED VALUE OR NOT???
  # ### Need to modify regex if not
  # spike_pssm <- str_match(spike_pssm_raw, "(?s).*lambdaUngapped\\s*\\{\\s*(.*)\\},\\s*kappaUngapped.*")
  # 
  # spike_pssm <- spike_pssm[2]
  # spike_pssm <- str_split(spike_pssm, ",\\s*", simplify=T)
  # # lambda <- as.numeric(spike_pssm[1])
  # # exponent <- as.numeric(spike_pssm[3])
  # # 
  # # #convert log likelihoods to likelihoods
  # # pssm_all <- pssm_all %>% 
  # #   mutate(Likelihood = exp(LogLikelihood*lambda*10^exponent))
  
  # join to master
  master_conservation <- inner_join(master_conservation, pssm_all)
  
  # calculate conservation score
  # equations in https://doi.org/10.1016/j.jmb.2004.08.022
  
  # add single replacement column to test against later
  master_conservation$Replacement <- paste(aa321(master_conservation$WTAA),master_conservation$Resno,master_conservation$ReplacementAA,sep="")
  
  # # Read Rosetta scores from processed CSVs -----------------------------------------------------
  
  setwd(paste("/Users/user/Documents/SARS-CoV-2/Structure/Rosetta/",structure,"_full_run",sep=""))
  
  # hard_rep_scores_min <- read_csv("scores_fixbb_6vxx_hardrep_minimise.csv")
  #
  # hard_rep_scores_min <- hard_rep_scores_min %>%
  #   select(!c(Score,WTScore)) %>%
  #   rename(deltaRosetta_hard = deltaRosetta)
  #
  # master_conservation <- inner_join(master_conservation,hard_rep_scores_min)
  
  soft_rep_scores_min <- read_csv(paste("scores_fixbb_",structure,"_softrep_min_dif_ML.csv",sep=""))
  
  soft_rep_scores_min <- soft_rep_scores_min %>%
    dplyr::rename(Resno = res_nums, ReplacementAA = replacement_aa)
  
  # replicate each row 3 times
  soft_rep_scores_min <- soft_rep_scores_min %>% dplyr::slice(rep(1:n(), each = 3))
  # add a chain column
  soft_rep_scores_min$Chain <- rep(c("A","B","C"),times=nrow(soft_rep_scores_min)/3)
  
  
  master_conservation <- inner_join(master_conservation,soft_rep_scores_min)
  
  
  # Read in epitope predictions ---------------------------------------------
  setwd("/Users/user/Documents/SARS-CoV-2/Epitope_mapping")
  
  epitope_tbl <- read_xlsx(paste(structure,"_epitopes_Saini.xlsx",sep=""))
  
  colnames(epitope_tbl) <- c("WTAA","Resno","Chain","pred_epitope_rank","obs_epitope_prev")
  
  # replicate each row 20 times
  epitope_tbl <- epitope_tbl[rep(1:nrow(epitope_tbl), each = 20), ]
  
  # add a replacement AA column
  epitope_tbl$ReplacementAA <- rep(unique(master_conservation$ReplacementAA),times=nrow(epitope_tbl)/20)
  
  master_conservation <- inner_join(master_conservation,epitope_tbl)
  # Read matrix of observed substitutions -----------------------------------
  
  setwd("/Users/user/Documents/Case study paper/1mo_ensemble/case_study_filtering_allmonths_12mo_50/spike_filtered_NOV_22")
  # if calculating full range of test sets switch code below
  for(i in c(1,3,6,12)){
  # for(i in c(1,3)){
  ### TRAINING SET
  obs_matrix <- read_csv(paste0("TEST_NOV_22_month",i,"_consensus_matrix.csv"))
  
  obs_matrix <- obs_matrix %>% 
    dplyr::rename("Resno" = "...1") %>% 
    select(-'*')
  
  # calculate observed frequencies as proportions
  obs_matrix_prop <- obs_matrix
  
  # sum total observed at each position
  obs_matrix_prop <- obs_matrix_prop %>%
    select(-Resno) %>% 
    mutate(Total_sum = rowSums(across(where(is.numeric))))
  
  # calculate proportion and drop sum column
  obs_matrix_prop <- obs_matrix_prop %>%
    mutate_at(vars(1:21),.funs=~((./Total_sum)*100)) %>% 
    select(-Total_sum)
  
  # add back resno and pivot long to join to master table
  obs_matrix_prop$Resno <- 1:nrow(obs_matrix_prop)
  
  obs_matrix_prop <- obs_matrix_prop %>% 
    pivot_longer(cols=-Resno,names_to = "ReplacementAA", values_to = "Observed_probability") %>% 
    filter(ReplacementAA != "X")
  
  check_sums <- obs_matrix_prop %>% group_by(Resno) %>% summarise(total = sum(Observed_probability))
  
  #join to master
  master_conservation_train <- left_join(master_conservation, obs_matrix_prop)
  
  setwd(paste("./Master_tables/",structure,sep=""))
  
  # write master conservation file
  write_xlsx(master_conservation_train, paste0(structure,"_spike_master_scores_ML_TEST_month",i,"_NOV_22.xlsx"))
  setwd("..")
  setwd("..")
  }
}

structure_list <- c("6vxx","7lws","7v7n","7wp9")

# note we expect this to throw errors for dates that do not include all 4 structures
# as long as structures are in chronological order this is ok
lapply(structure_list,master_table_writer)
