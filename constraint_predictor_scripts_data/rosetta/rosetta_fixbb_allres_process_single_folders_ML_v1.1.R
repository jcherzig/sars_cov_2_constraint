library(tidyverse)
library(viridis)
library(bio3d)
library(writexl)

# read rosetta output with each position in separate folder
# set wd to file containing all output folders and relevant PDB file

# Read PDB file to create list of mutations -------------------------------


# read pdb file
## REMEMBER to set the pdb_code object
# read pdb file
raw_pdb <- read.pdb("6vxx_relax_soft.pdb")

pdb_frame <- raw_pdb$atom

pdb_frame <- as_tibble(pdb_frame)

#remove NAG
pdb_frame <- pdb_frame %>% 
  filter(resid != "NAG")

# filter for Calpha atoms only
pdb_frame <- pdb_frame %>% 
  filter(elety == "CA")

# convert to 1 letter aa codes
aa_vec <- pdb_frame$resid
aa_vec <- aa321(aa_vec)

# create output vector
dir_list <- read_tsv("dir_names_for_jobarray.txt",col_names=F)

dir_list <- do.call(paste0, dir_list)

# residue numbers table

res_nums <- str_extract(dir_list, "(?<=6vxx_)\\d+")

# replacement amino acid table

replacement_aa <- str_extract(dir_list, "(?<=PIKAA_).*")

# process and concatenate data --------------------------------------------

# read scanning output
### NOTE: try wrapper in function so missing folders will not halt loop but will
### print into output vector. Check for errors
# rosetta_reader <- function(code){
#   #set working directory
#   try({setwd(as.character(code))
# 
#   #read in scanning output
#   raw_fixbb <- read_table("score.sc",skip=1)
# 
#   #select total score column only
#   raw_fixbb <- raw_fixbb %>%
#     select(total_score)
# 
#   colnames(raw_fixbb) <- "Score"
# 
#   setwd("..")
#   },silent = F)
#   try(return(raw_fixbb),silent = F)
# }

rosetta_reader <- function(code){
  #set working directory
  setwd(code)

  #read in scanning output
  raw_fixbb <- read_table("score.sc",skip=1)
  
  raw_fixbb <- select(raw_fixbb, -c("SCORE:","description","X20"))

  setwd("..")
  return(raw_fixbb)
}

allres_fixbb <- as_tibble(do.call(rbind,lapply(dir_list, rosetta_reader)))

write_csv(allres_fixbb,"scores_fixbb_6vxx_softrep_minimise_ML.csv")


# read in wildtype (relaxed) structure score for soft and hard reps
setwd("/Users/user/OneDrive - The University of Manchester/SARS-CoV-2/Structure/Rosetta/6vxx_score")


# ### for hard rep
# hard_score <- read_table("hardrep_score.sc")
# #select total score column only
# hard_score <- hard_score %>%
#   select(score)
# 
# 
# allres_fixbb_output <- allres_fixbb_output %>%
#   mutate(Hard_WT_Score = hard_score$score)

### for soft rep
soft_score <- read_table("softrep_score.sc")
# select total score column only

soft_score <- soft_score %>%rename(total_score = score)

names_vec <- colnames(allres_fixbb)

soft_score <- soft_score %>%
  select(all_of(names_vec))

allres_fixbb_dif <- allres_fixbb %>% mutate(total_score = total_score - soft_score$total_score,
                                             dslf_ca_dih = dslf_ca_dih - soft_score$dslf_ca_dih,
                                             dslf_cs_ang = dslf_cs_ang - soft_score$dslf_cs_ang,
                                             dslf_ss_dih = dslf_ss_dih - soft_score$dslf_ss_dih,
                                             dslf_ss_dst = dslf_ss_dst - soft_score$dslf_ss_dst,
                                             fa_atr = fa_atr - soft_score$fa_atr,
                                             fa_dun = fa_dun - soft_score$fa_dun,
                                             fa_pair = fa_pair - soft_score$fa_pair,
                                             fa_rep = fa_rep - soft_score$fa_rep,
                                             fa_sol = fa_sol - soft_score$fa_sol,
                                             hbond_bb_sc = hbond_bb_sc - soft_score$hbond_bb_sc,
                                             hbond_lr_bb = hbond_lr_bb - soft_score$hbond_lr_bb,
                                             hbond_sc = hbond_sc - soft_score$hbond_sc,
                                             hbond_sr_bb = hbond_sr_bb - soft_score$hbond_sr_bb,
                                             p_aa_pp = p_aa_pp - soft_score$p_aa_pp,
                                             pro_close = pro_close - soft_score$pro_close,
                                             ref = ref - soft_score$ref)


allres_fixbb_output <- cbind(allres_fixbb_dif,res_nums,replacement_aa)


setwd("/Users/user/OneDrive - The University of Manchester/SARS-CoV-2/Structure/Rosetta/6vxx_full_run")


# print CSV of output
write_csv(allres_fixbb_output,"scores_fixbb_6vxx_softrep_min_dif_ML.csv")
