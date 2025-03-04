require(Biostrings)
require(tidyverse)
require(bio3d)
require(writexl)

# script to read buildmodel 'Dif' output and clean up to give ddG

# build list of residue numbers from PDB
# read pdb file
raw_pdb <- read.pdb("6vxx_Repair.PDB")

pdb_frame <- raw_pdb$atom
rm(raw_pdb)

pdb_frame <- as_tibble(pdb_frame)

#remove NAG
pdb_frame <- pdb_frame %>% 
  filter(resid != "NAG")

# filter for Calpha atoms only
pdb_frame <- pdb_frame %>% 
  filter(elety == "CA")

# identify longest chain (i.e. highest coverage)
a_length <- nrow(pdb_frame %>% 
  filter(chain == "A"))

b_length <- nrow(pdb_frame %>% 
  filter(chain == "B"))

c_length <- nrow(pdb_frame %>% 
  filter(chain == "C"))

# filter for longest chain
pdb_frame <- pdb_frame %>% 
  filter(chain == c("A","B","C")[which.max(c(a_length,b_length,c_length))]
)

resno_vec <- pdb_frame$resno

# read in output and remove all columns except ddG (difference between WT and 
# corresponding mutant)
clean_and_write <- function(folder,aa_sub){
setwd(folder)

raw_dif <- read_tsv("Dif_6vxx_Repair.fxout", skip = 8)

# remove any columns that are all 0
clean_dif <- raw_dif[, colSums(raw_dif != 0) > 0]

# bind to residue number vector
clean_dif$Resno <- resno_vec
clean_dif$ReplacementAA <- as.character(aa_sub)

#remove unnecessary column
clean_dif <- select(clean_dif, -Pdb)



colnames(clean_dif) <- c("DeltaDeltaG","Backbone_Hbond","Sidechain_Hbond","Van_der_Waals","Electrostatics","Solvation_Polar",
                         "Solvation_Hydrophobic","Van_der_Waals_clashes","entropy_sidechain","entropy_mainchain","torsional_clash",
                         "backbone_clash","helix_dipole","disulfide","electrostatic_kon","energy_Ionisation","Resno","ReplacementAA")

# write output in individual folders
write_csv(clean_dif, paste(aa_sub,"_clean_ddG_ML.csv",sep=""))

setwd("..")
}

# list of files
file_list <- c("buildmodel_output_A","buildmodel_output_C","buildmodel_output_D","buildmodel_output_E","buildmodel_output_F","buildmodel_output_G","buildmodel_output_H",
  "buildmodel_output_I","buildmodel_output_K","buildmodel_output_L","buildmodel_output_M","buildmodel_output_N","buildmodel_output_P","buildmodel_output_Q",
  "buildmodel_output_R","buildmodel_output_S","buildmodel_output_T","buildmodel_output_V","buildmodel_output_W","buildmodel_output_Y")

aa_sub <- c("A","C","D","E","F","G","H",
              "I","K","L","M","N","P","Q",
              "R","S","T","V","W","Y")


mapply(FUN=clean_and_write,folder=file_list,aa_sub=aa_sub)

# set up tibble to join to
clean_combine <- tibble("DeltaDeltaG" = NA,"Backbone_Hbond" = NA,"Sidechain_Hbond" = NA,"Van_der_Waals" = NA,"Electrostatics" = NA,"Solvation_Polar" = NA,
                        "Solvation_Hydrophobic" = NA,"Van_der_Waals_clashes" = NA,"entropy_sidechain" = NA,"entropy_mainchain" = NA,"torsional_clash" = NA,
                        "backbone_clash" = NA,"helix_dipole" = NA,"disulfide" = NA,"electrostatic_kon" = NA,"energy_Ionisation" = NA,"Resno" = NA,"ReplacementAA" = NA)

# function to read in all files just written, combine and write into single master
read_and_combine <- function(folder,aa_sub){
  setwd(folder)
  
  clean_indiv <- read_csv(paste(aa_sub,"_clean_ddG_ML.csv",sep=""),
                          col_types=cols(DeltaDeltaG="?",Backbone_Hbond="?",Sidechain_Hbond="?",Van_der_Waals="?",Electrostatics="?",
                                         Solvation_Polar="?",Solvation_Hydrophobic="?",Van_der_Waals_clashes="?",entropy_sidechain="?",
                                         entropy_mainchain="?",torsional_clash="?",backbone_clash="?",helix_dipole="?",disulfide="?",
                                         electrostatic_kon="?",energy_Ionisation="?",Resno="?",ReplacementAA="c"))
  
  # join each individual file
  clean_combine <<- rbind(clean_combine,clean_indiv)
  
  setwd("..")
  
  return(clean_combine)
}


mapply(read_and_combine,folder=file_list,aa_sub=aa_sub)

clean_combine <- clean_combine[-1,]

write_csv(clean_combine, "buildmodel_output_master_ML.csv")
