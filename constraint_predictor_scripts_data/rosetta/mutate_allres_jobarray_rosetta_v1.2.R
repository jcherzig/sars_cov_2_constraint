require(Biostrings)
require(tidyverse)
require(bio3d)
require(writexl)

#function to read pdb file and create residue files for all individual residues

#WARNING: ensure correct wd is set before code is run 
##(folder containing relevant pdb file only)
  
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
  
  fixbb_output <- tibble(.rows=length(aa_vec))

  # create output vector
  for(i in 1:20){
  fixbb_output[i] <- paste(pdb_frame$resno," ",pdb_frame$chain," PIKAA ",(unique(aa_vec))[i],sep="")
  }
  
  # convert to strings
  fixbb_output_A <- do.call(paste0, fixbb_output[,1])
  fixbb_output_Y <- do.call(paste0, fixbb_output[,2])
  fixbb_output_T <- do.call(paste0, fixbb_output[,3])
  fixbb_output_N <- do.call(paste0, fixbb_output[,4])
  fixbb_output_S <- do.call(paste0, fixbb_output[,5])
  fixbb_output_F <- do.call(paste0, fixbb_output[,6])
  fixbb_output_R <- do.call(paste0, fixbb_output[,7])
  fixbb_output_G <- do.call(paste0, fixbb_output[,8])
  fixbb_output_V <- do.call(paste0, fixbb_output[,9])
  fixbb_output_P <- do.call(paste0, fixbb_output[,10])
  fixbb_output_D <- do.call(paste0, fixbb_output[,11])
  fixbb_output_K <- do.call(paste0, fixbb_output[,12])
  fixbb_output_L <- do.call(paste0, fixbb_output[,13])
  fixbb_output_H <- do.call(paste0, fixbb_output[,14])
  fixbb_output_Q <- do.call(paste0, fixbb_output[,15])
  fixbb_output_W <- do.call(paste0, fixbb_output[,16])
  fixbb_output_I <- do.call(paste0, fixbb_output[,17])
  fixbb_output_E <- do.call(paste0, fixbb_output[,18])
  fixbb_output_C <- do.call(paste0, fixbb_output[,19])
  fixbb_output_M <- do.call(paste0, fixbb_output[,20])
  
  
  # define utility function to bind vectors with uneven length (fill with NAs)
  bind_cols_fill <- function(df_list) {
    
    max_rows <- map_int(df_list, nrow) %>% max()
    
    map(df_list, function(df) {
      if(nrow(df) == max_rows) return(df)
      first <- names(df)[1] %>% sym()
      df %>% add_row(!!first := rep("", max_rows - nrow(df)))
    }) %>% bind_cols()
  }
  
  # function to separate each list of mutants by chain
  sep_by_col <-  function(x){
    
    #define empty vectors
    A_vec <-vector()
    B_vec <-vector()
    C_vec <-vector()
    
    # separate into vectors by chain
    for(i in 1:length(x)){  if(str_detect(x[i], "\\sA\\sPIK.")==T){
      A_vec[i] <- x[i]}
    }
    
    for(j in 1:length(x)){if(str_detect(x[j], "\\sB\\sPIK.")==T){
      B_vec[j] <- x[j]}
    }
    
    for(k in 1:length(x)){if(str_detect(x[k], "\\sC\\sPIK.")==T){
      C_vec[k] <- x[k]}
    }
    
    B_vec <- B_vec[!is.na(B_vec)]
    C_vec <- C_vec[!is.na(C_vec)]
    
    # bind to create tibble containing all mutations in rows by chain
    bind_cols_fill(list(tibble(A_vec),tibble(B_vec),tibble(C_vec)))
  }
  
  # apply the function to every substitution list
  scan_list <- list(fixbb_output_A,
                    fixbb_output_Y,
                    fixbb_output_T,
                    fixbb_output_N,
                    fixbb_output_S,
                    fixbb_output_F,
                    fixbb_output_R,
                    fixbb_output_G,
                    fixbb_output_V,
                    fixbb_output_P,
                    fixbb_output_D,
                    fixbb_output_K,
                    fixbb_output_L,
                    fixbb_output_H,
                    fixbb_output_Q,
                    fixbb_output_W,
                    fixbb_output_I,
                    fixbb_output_E,
                    fixbb_output_C,
                    fixbb_output_M)
  
  
  printme <- sapply(scan_list,sep_by_col)
  # convert to tibble and transpose so chains are columns and substitutions are rows
  printme_t <- as_tibble(printme)
  printme_t <- t.data.frame(printme_t)
  
  
  # create individual tibbles for every mutant
  fixbb_output_A <- tibble("A_vec"=unlist(printme[1,][1]),"B_vec"=unlist(printme[2,][1]),"C_vec"=unlist(printme[3,][1]))
  fixbb_output_Y <- tibble("A_vec"=unlist(printme[1,][2]),"B_vec"=unlist(printme[2,][2]),"C_vec"=unlist(printme[3,][2]))
  fixbb_output_T <- tibble("A_vec"=unlist(printme[1,][3]),"B_vec"=unlist(printme[2,][3]),"C_vec"=unlist(printme[3,][3]))
  fixbb_output_N <- tibble("A_vec"=unlist(printme[1,][4]),"B_vec"=unlist(printme[2,][4]),"C_vec"=unlist(printme[3,][4]))
  fixbb_output_S <- tibble("A_vec"=unlist(printme[1,][5]),"B_vec"=unlist(printme[2,][5]),"C_vec"=unlist(printme[3,][5]))
  fixbb_output_F <- tibble("A_vec"=unlist(printme[1,][6]),"B_vec"=unlist(printme[2,][6]),"C_vec"=unlist(printme[3,][6]))
  fixbb_output_R <- tibble("A_vec"=unlist(printme[1,][7]),"B_vec"=unlist(printme[2,][7]),"C_vec"=unlist(printme[3,][7]))
  fixbb_output_G <- tibble("A_vec"=unlist(printme[1,][8]),"B_vec"=unlist(printme[2,][8]),"C_vec"=unlist(printme[3,][8]))
  fixbb_output_V <- tibble("A_vec"=unlist(printme[1,][9]),"B_vec"=unlist(printme[2,][9]),"C_vec"=unlist(printme[3,][9]))
  fixbb_output_P <- tibble("A_vec"=unlist(printme[1,][10]),"B_vec"=unlist(printme[2,][10]),"C_vec"=unlist(printme[3,][10]))
  fixbb_output_D <- tibble("A_vec"=unlist(printme[1,][11]),"B_vec"=unlist(printme[2,][11]),"C_vec"=unlist(printme[3,][11]))
  fixbb_output_K <- tibble("A_vec"=unlist(printme[1,][12]),"B_vec"=unlist(printme[2,][12]),"C_vec"=unlist(printme[3,][12]))
  fixbb_output_L <- tibble("A_vec"=unlist(printme[1,][13]),"B_vec"=unlist(printme[2,][13]),"C_vec"=unlist(printme[3,][13]))
  fixbb_output_H <- tibble("A_vec"=unlist(printme[1,][14]),"B_vec"=unlist(printme[2,][14]),"C_vec"=unlist(printme[3,][14]))
  fixbb_output_Q <- tibble("A_vec"=unlist(printme[1,][15]),"B_vec"=unlist(printme[2,][15]),"C_vec"=unlist(printme[3,][15]))
  fixbb_output_W <- tibble("A_vec"=unlist(printme[1,][16]),"B_vec"=unlist(printme[2,][16]),"C_vec"=unlist(printme[3,][16]))
  fixbb_output_I <- tibble("A_vec"=unlist(printme[1,][17]),"B_vec"=unlist(printme[2,][17]),"C_vec"=unlist(printme[3,][17]))
  fixbb_output_E <- tibble("A_vec"=unlist(printme[1,][18]),"B_vec"=unlist(printme[2,][18]),"C_vec"=unlist(printme[3,][18]))
  fixbb_output_C <- tibble("A_vec"=unlist(printme[1,][19]),"B_vec"=unlist(printme[2,][19]),"C_vec"=unlist(printme[3,][19]))
  fixbb_output_M <- tibble("A_vec"=unlist(printme[1,][20]),"B_vec"=unlist(printme[2,][20]),"C_vec"=unlist(printme[3,][20]))
  
  
  # new list with new tibbles
  scan_list <- list(fixbb_output_A,
                    fixbb_output_Y,
                    fixbb_output_T,
                    fixbb_output_N,
                    fixbb_output_S,
                    fixbb_output_F,
                    fixbb_output_R,
                    fixbb_output_G,
                    fixbb_output_V,
                    fixbb_output_P,
                    fixbb_output_D,
                    fixbb_output_K,
                    fixbb_output_L,
                    fixbb_output_H,
                    fixbb_output_Q,
                    fixbb_output_W,
                    fixbb_output_I,
                    fixbb_output_E,
                    fixbb_output_C,
                    fixbb_output_M)
  
  
  ###############################
  
  # create list of nat amino acids
  aa_list <- c('A','C','D','E','F','G','H','I','K','L','M',
               'N','P','Q','R','S','T','V','W','Y')
  
  
  resno_max_vec <- pdb_frame %>% 
    group_by(chain) %>% 
    summarise(max(length(resno)))
  
# B chain longest (contains all residue positions) so we use that for naming
  
pdb_frame_B <- pdb_frame %>% 
  filter(chain == "B")
  
  # create naming output vectors
  for(i in 1:20){
  fixbb_output_names_A <- paste(pdb_frame_B$resno, "PIKAA","A", sep="_")
 }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_C <- paste(pdb_frame_B$resno, "PIKAA","C", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_D <- paste(pdb_frame_B$resno, "PIKAA","D", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_E <- paste(pdb_frame_B$resno, "PIKAA","E", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_F <- paste(pdb_frame_B$resno, "PIKAA","F", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_G <- paste(pdb_frame_B$resno, "PIKAA","G", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_H <- paste(pdb_frame_B$resno, "PIKAA","H", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_I <- paste(pdb_frame_B$resno, "PIKAA","I", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_K <- paste(pdb_frame_B$resno, "PIKAA","K", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_L <- paste(pdb_frame_B$resno, "PIKAA","L", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_M <- paste(pdb_frame_B$resno, "PIKAA","M", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_N <- paste(pdb_frame_B$resno, "PIKAA","N", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_P <- paste(pdb_frame_B$resno, "PIKAA","P", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_Q <- paste(pdb_frame_B$resno, "PIKAA","Q", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_R <- paste(pdb_frame_B$resno, "PIKAA","R", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_S <- paste(pdb_frame_B$resno, "PIKAA","S", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_T <- paste(pdb_frame_B$resno, "PIKAA","T", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_V <- paste(pdb_frame_B$resno, "PIKAA","V", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_W <- paste(pdb_frame_B$resno, "PIKAA","W", sep="_")
  }
  # create naming output vectors
  for(i in 1:20){
    fixbb_output_names_Y <- paste(pdb_frame_B$resno, "PIKAA","Y", sep="_")
  }
  
  
config_writer <- function(content,names){
  # print csv
  for(i in 1:nrow(fixbb_output_A)){
      cat(paste("NATAA
start
",content[i,1],"
",content[i,2],"
",content[i,3],sep=""),
          file=paste("6vxx_",
                     names[i],
                     ".cfg",
                     sep=""),
          sep="")
    }
}



# new list with new tibbles
name_list <- list(fixbb_output_names_A,
                  fixbb_output_names_Y,
                  fixbb_output_names_T,
                  fixbb_output_names_N,
                  fixbb_output_names_S,
                  fixbb_output_names_F,
                  fixbb_output_names_R,
                  fixbb_output_names_G,
                  fixbb_output_names_V,
                  fixbb_output_names_P,
                  fixbb_output_names_D,
                  fixbb_output_names_K,
                  fixbb_output_names_L,
                  fixbb_output_names_H,
                  fixbb_output_names_Q,
                  fixbb_output_names_W,
                  fixbb_output_names_I,
                  fixbb_output_names_E,
                  fixbb_output_names_C,
                  fixbb_output_names_M)

mapply(config_writer, content=scan_list,names=name_list)

#function to generate jobscript file and list of directories

#NOTE must run allres_mutator function first to generate cfg files

  #create jobscript to run each config file one by one
  filelist <- list.files(pattern="*.cfg")
  #remove the .cfg
  namelist <- str_sub(filelist,1,nchar(filelist)-4)
  
  full_list <- tibble(names_list = NA, number = NA, .rows=length(filelist))
  
  # create list of nat amino acids
  aa_list <- c('A','C','D','E','F','G','H','I','K','L','M',
               'N','P','Q','R','S','T','V','W','Y')
  
  for(i in 1:length(filelist)){
      full_list$names_list[i] <- paste("fixbb -s 6vxx_relax_soft.pdb -resfile ", filelist[i]," -nstruct 1 -ex1 -ex2 -use_input_sc -out:path:all ./",namelist[i]," -out:pdb_gz",sep="")
      full_list$number[i] <- i
      i+1
    }
  
  jobscript_out <- paste(full_list$names_list, collapse="\n")
  
  # print list of config files 1 per line for bash scripting
  cat(filelist,file="6vxx_config_filelist.txt",sep="
")
  
  # print list of directory names 1 per line for bash scripting
  cat(namelist,file="dir_names_for_jobarray.txt",sep="
")
  ###add below to fixbb command line to use soft rep weights
  # -score:weights soft_rep_design
  
  ###add below to fixbb command line to add minimization step
  # -minimize_sidechains
  
  cat(paste
      ("#!/bin/bash --login

#$ -cwd
#$ -t 1-",nrow(full_list),"

INFILE=`awk \"NR==$SGE_TASK_ID\" 6vxx_config_filelist.txt`
OUTDIR=`awk \"NR==$SGE_TASK_ID\" dir_names_for_jobarray.txt`

module load apps/binapps/rosetta/2020.08
fixbb -s 6vxx_relax_soft.pdb -resfile $INFILE -nstruct 1 -ex1 -ex2 -use_input_sc -score:weights soft_rep_design -minimize_sidechains -out:path:all $OUTDIR -out:pdb_gz"
        ,sep=""),
file="jobscript_6vxx_fixbb_AA_singles.job",sep="")

# print the list of directory names to pass to xargs
  
  cat(namelist,file="dir_names_for_fixbb_singles.txt")
