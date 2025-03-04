library(tidyverse)
library(writexl)
library(stringr)
library(viridis)
library(bio3d)


# MAC paths
setwd("/Users/user/OneDrive - The University of Manchester/SARS-CoV-2/BLAST/DELTABLAST")

# read pdb file 
raw_pdb <- read.pdb("6vxx")

pdb_frame <- raw_pdb$atom

pdb_frame <- as_tibble(pdb_frame)

# filter for Calpha atoms only
pdb_frame <- pdb_frame %>% 
  filter(elety == "CA")

# read sequence from the FASTA file used for DELTABLAST run
spike_seq <- read.fasta("rcsb_pdb_6VXX.fasta")
spike_seq <- spike_seq$ali

# make tibble and name
spike_seq <- as_tibble(spike_seq)

spike_seq <- spike_seq %>% 
  pivot_longer(cols = everything(), names_to = "Resno",values_to = "WTAA")

spike_seq$WTAA <-aa123(spike_seq$WTAA)
# create list of residues present in PDB for each chain
# FASTA is starting on different residue number than PDB here so we need to renumber
spike_seq$Resno <- -18:1262


# read PSSM file from DELTA BLAST
spike_pssm_raw <- read_file("dblast_6VXX_pssm.out")

# extract likelihood scores
spike_pssm <- str_match(spike_pssm_raw, "(?s).*finalData\\s*\\{\\s*scores\\s*\\{(.*)\\},\\s*lambda\\s*\\{.*")

spike_pssm <- spike_pssm[2]
spike_pssm <- str_split(spike_pssm, ",\\n\\s*", simplify=T)

#clean first entry
spike_pssm[1] <- gsub(".*\\n\\s*(.+).*", "\\1", spike_pssm[1])

# convert to tibble and pivot
spike_pssm <- as_tibble(spike_pssm)
spike_pssm <- pivot_longer(spike_pssm, cols = everything(), names_to = NULL, values_to = "LogLikelihood")

# add amino acid identities to new column
aa_list <- c('-','A','B','C','D','E','F','G','H','I','K','L','M',
'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
 'O', 'J')

spike_pssm$ReplacementAA <- rep(aa_list, times = length(spike_pssm$LogLikelihood)/28)

spike_pssm <- spike_pssm %>%  
  filter(!ReplacementAA %in% c("B","Z","-", "X", "U", "*", "O", "J"))

spike_pssm$LogLikelihood <-as.numeric(spike_pssm$LogLikelihood)

# fill the WTAA and residue number columns from the spike_seq tibble already 
# checked against the PDB file
spike_pssm$WTAA <- rep(spike_seq$WTAA, each = 20)
spike_pssm$Resno <- rep(spike_seq$Resno, each = 20)

# write simple processed files as csv
write_csv(spike_pssm, "PROCESSED_DELTABLAST_full_6vxx_spike_pssm.csv")

# filter the PDB file to make 3 lists of residues present in each chain
A_chain_list <- pdb_frame %>% 
  filter(chain == "A")
A_chain_list <- A_chain_list$resno

B_chain_list <- pdb_frame %>% 
  filter(chain == "B")
B_chain_list <- B_chain_list$resno

C_chain_list <- pdb_frame %>% 
  filter(chain == "C")
C_chain_list <- C_chain_list$resno


# now filter the PSSM for residues present in each chain of the PDB file
spike_pssm_A <- spike_pssm %>% 
  filter(Resno %in% A_chain_list)

spike_pssm_B <- spike_pssm %>% 
  filter(Resno %in% B_chain_list)

spike_pssm_C <- spike_pssm %>% 
  filter(Resno %in% C_chain_list)

# NOTE: we are producing 3 separate PSSM files for easy joining to the ESST
# and FoldX analysis which is performed on each separate chain. DELTABLAST is
# NOT performed on chains separately (as they have identical sequences) so all
# values are identical across the 3 chains

# add chain value col and bind
spike_pssm_A$Chain <- "A"
spike_pssm_B$Chain <- "B"
spike_pssm_C$Chain <- "C"

spike_pssm <- bind_rows(spike_pssm_A, spike_pssm_B, spike_pssm_C)


# plot heatmap
p <- ggplot(spike_pssm, aes(x = Resno, y = ReplacementAA, fill = LogLikelihood)) +
  geom_tile() +
  theme_bw(base_size = 12) +
  xlab("Canonical Residue") +
  ylab("Replacement residue likelihood") +
  scale_x_continuous(breaks = seq(from = 0, to = 667, by=25)) +
  scale_fill_gradient2(low="black", mid="white", high = "red", midpoint = 0, name = "Likelihood")
p

# write processed files for each chain as csv
write_csv(spike_pssm, "PROCESSED_DELTABLAST_bychain_6vxx_spike_pssm.csv")




