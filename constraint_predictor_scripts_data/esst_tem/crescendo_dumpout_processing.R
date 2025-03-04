library(tidyverse)
library(viridis)
library(writexl)
library(bio3d)

# read reference spike structure
ref_pdb <- read.pdb("6vxx_Repair.pdb")

# create replacement vectors with correct residue numbers
ref_pdb <- ref_pdb$atom

ref_pdb_A <- ref_pdb %>% filter(chain=="A",elety=="CA")
ref_pdb_A_unique <- ref_pdb_A$resno

ref_pdb_B <- ref_pdb %>% filter(chain=="B",elety=="CA")
ref_pdb_B_unique <- ref_pdb_B$resno

ref_pdb_C <- ref_pdb %>% filter(chain=="C",elety=="CA")
ref_pdb_C_unique <- ref_pdb_C$resno

# Read, process and plot mutational probability outputs--------
raw_probs_C <- read_table("6vxxC.probs", col_names = F)

colnames(raw_probs_C) <- c("resno", "res", "A", "C", "D", "E", "F", "G", "H", "I",
                         "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# pivot longer
processed_probs_C <- pivot_longer(raw_probs_C, -c("resno", "res"), names_to = "replacement_res", 
                          values_to = "probability")

processed_probs_C$resno <- as.numeric(processed_probs_C$resno)
processed_probs_C$probability <- as.numeric(processed_probs_C$probability)

# replace with correct residue numbers
processed_probs_C$resno <- rep(ref_pdb_C_unique, each = 20)

processed_probs_A$chain <- "A"
processed_probs_B$chain <- "B"
processed_probs_C$chain <- "C"

processed_probs_full <- bind_rows(processed_probs_A, processed_probs_B, processed_probs_C)
write_csv(processed_probs_full, "PROCESSED_crescendo_dumpouts_6vxx.csv")
