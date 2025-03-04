require(Biostrings)
require(tidyverse)
require(bio3d)
require(writexl)

#function to read pdb file and create buildmodel ready list of 
##all residues mutated

#WARNING: ensure correct wd is set before code is run 
##(folder containing relevant pdb file only)

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
  
  # convert to 1 letter aa codes
  aa_vec <- pdb_frame$resid
  aa_vec <- aa321(aa_vec)
  
  posscan_output <- tibble(.rows=nrow(pdb_frame))
  
  # create mutation list containing all chains and individual mutants (3x20)
  for(i in 1:20){
    posscan_output[i]<- paste(aa_vec, pdb_frame$chain, pdb_frame$resno,(unique(aa_vec))[i],sep="")
  }
  
  # convert to strings
  posscan_output_A <- do.call(paste0, posscan_output[,1])
  posscan_output_Y <- do.call(paste0, posscan_output[,2])
  posscan_output_T <- do.call(paste0, posscan_output[,3])
  posscan_output_N <- do.call(paste0, posscan_output[,4])
  posscan_output_S <- do.call(paste0, posscan_output[,5])
  posscan_output_F <- do.call(paste0, posscan_output[,6])
  posscan_output_R <- do.call(paste0, posscan_output[,7])
  posscan_output_G <- do.call(paste0, posscan_output[,8])
  posscan_output_V <- do.call(paste0, posscan_output[,9])
  posscan_output_P <- do.call(paste0, posscan_output[,10])
  posscan_output_D <- do.call(paste0, posscan_output[,11])
  posscan_output_K <- do.call(paste0, posscan_output[,12])
  posscan_output_L <- do.call(paste0, posscan_output[,13])
  posscan_output_H <- do.call(paste0, posscan_output[,14])
  posscan_output_Q <- do.call(paste0, posscan_output[,15])
  posscan_output_W <- do.call(paste0, posscan_output[,16])
  posscan_output_I <- do.call(paste0, posscan_output[,17])
  posscan_output_E <- do.call(paste0, posscan_output[,18])
  posscan_output_C <- do.call(paste0, posscan_output[,19])
  posscan_output_M <- do.call(paste0, posscan_output[,20])
  
  # define utility function to bind vectors with uneven length (fill with NAs)
  bind_cols_fill <- function(df_list) {
    
    max_rows <- map_int(df_list, nrow) %>% max()
    
    map(df_list, function(df) {
      if(nrow(df) == max_rows) return(df)
      first <- names(df)[1] %>% sym()
      df %>% add_row(!!first := rep(NA, max_rows - nrow(df)))
    }) %>% bind_cols()
  }

  # function to separate each list of mutants by chain
sep_by_col <-  function(x){
  
  #define empty vectors
A_vec <-vector()
B_vec <-vector()
C_vec <-vector()

# separate into vectors by chain
for(i in 1:length(x)){  if(str_detect(x[i], "^(.{1})A.")==T){
  A_vec[i] <- x[i]}
}

for(j in 1:length(x)){if(str_detect(x[j], "^(.{1})B.")==T){
    B_vec[j] <- x[j]}
}

for(k in 1:length(x)){if(str_detect(x[k], "^(.{1})C.")==T){
  C_vec[k] <- x[k]}
}

B_vec <- B_vec[!is.na(B_vec)]
C_vec <- C_vec[!is.na(C_vec)]

# bind to create tibble containing all mutations in rows by chain
bind_cols_fill(list(tibble(A_vec),tibble(B_vec),tibble(C_vec)))
}

# apply the function to every substitution list
scan_list <- list(posscan_output_A,
                  posscan_output_Y,
                  posscan_output_T,
                  posscan_output_N,
                  posscan_output_S,
                  posscan_output_F,
                  posscan_output_R,
                  posscan_output_G,
                  posscan_output_V,
                  posscan_output_P,
                  posscan_output_D,
                  posscan_output_K,
                  posscan_output_L,
                  posscan_output_H,
                  posscan_output_Q,
                  posscan_output_W,
                  posscan_output_I,
                  posscan_output_E,
                  posscan_output_C,
                  posscan_output_M)

printme <- sapply(scan_list,sep_by_col)
# convert to tibble and transpose so chains are columns and substitutions are rows
printme_t <- as_tibble(printme)
printme_t <- t.data.frame(printme_t)


# create individual tibbles for every mutant
posscan_output_A <- tibble("A_vec"=unlist(printme[1,][1]),"B_vec"=unlist(printme[2,][1]),"C_vec"=unlist(printme[3,][1]))
posscan_output_Y <- tibble("A_vec"=unlist(printme[1,][2]),"B_vec"=unlist(printme[2,][2]),"C_vec"=unlist(printme[3,][2]))
posscan_output_T <- tibble("A_vec"=unlist(printme[1,][3]),"B_vec"=unlist(printme[2,][3]),"C_vec"=unlist(printme[3,][3]))
posscan_output_N <- tibble("A_vec"=unlist(printme[1,][4]),"B_vec"=unlist(printme[2,][4]),"C_vec"=unlist(printme[3,][4]))
posscan_output_S <- tibble("A_vec"=unlist(printme[1,][5]),"B_vec"=unlist(printme[2,][5]),"C_vec"=unlist(printme[3,][5]))
posscan_output_F <- tibble("A_vec"=unlist(printme[1,][6]),"B_vec"=unlist(printme[2,][6]),"C_vec"=unlist(printme[3,][6]))
posscan_output_R <- tibble("A_vec"=unlist(printme[1,][7]),"B_vec"=unlist(printme[2,][7]),"C_vec"=unlist(printme[3,][7]))
posscan_output_G <- tibble("A_vec"=unlist(printme[1,][8]),"B_vec"=unlist(printme[2,][8]),"C_vec"=unlist(printme[3,][8]))
posscan_output_V <- tibble("A_vec"=unlist(printme[1,][9]),"B_vec"=unlist(printme[2,][9]),"C_vec"=unlist(printme[3,][9]))
posscan_output_P <- tibble("A_vec"=unlist(printme[1,][10]),"B_vec"=unlist(printme[2,][10]),"C_vec"=unlist(printme[3,][10]))
posscan_output_D <- tibble("A_vec"=unlist(printme[1,][11]),"B_vec"=unlist(printme[2,][11]),"C_vec"=unlist(printme[3,][11]))
posscan_output_K <- tibble("A_vec"=unlist(printme[1,][12]),"B_vec"=unlist(printme[2,][12]),"C_vec"=unlist(printme[3,][12]))
posscan_output_L <- tibble("A_vec"=unlist(printme[1,][13]),"B_vec"=unlist(printme[2,][13]),"C_vec"=unlist(printme[3,][13]))
posscan_output_H <- tibble("A_vec"=unlist(printme[1,][14]),"B_vec"=unlist(printme[2,][14]),"C_vec"=unlist(printme[3,][14]))
posscan_output_Q <- tibble("A_vec"=unlist(printme[1,][15]),"B_vec"=unlist(printme[2,][15]),"C_vec"=unlist(printme[3,][15]))
posscan_output_W <- tibble("A_vec"=unlist(printme[1,][16]),"B_vec"=unlist(printme[2,][16]),"C_vec"=unlist(printme[3,][16]))
posscan_output_I <- tibble("A_vec"=unlist(printme[1,][17]),"B_vec"=unlist(printme[2,][17]),"C_vec"=unlist(printme[3,][17]))
posscan_output_E <- tibble("A_vec"=unlist(printme[1,][18]),"B_vec"=unlist(printme[2,][18]),"C_vec"=unlist(printme[3,][18]))
posscan_output_C <- tibble("A_vec"=unlist(printme[1,][19]),"B_vec"=unlist(printme[2,][19]),"C_vec"=unlist(printme[3,][19]))
posscan_output_M <- tibble("A_vec"=unlist(printme[1,][20]),"B_vec"=unlist(printme[2,][20]),"C_vec"=unlist(printme[3,][20]))


# new list with new tibbles
scan_list <- list(posscan_output_A,
                  posscan_output_Y,
                  posscan_output_T,
                  posscan_output_N,
                  posscan_output_S,
                  posscan_output_F,
                  posscan_output_R,
                  posscan_output_G,
                  posscan_output_V,
                  posscan_output_P,
                  posscan_output_D,
                  posscan_output_K,
                  posscan_output_L,
                  posscan_output_H,
                  posscan_output_Q,
                  posscan_output_W,
                  posscan_output_I,
                  posscan_output_E,
                  posscan_output_C,
                  posscan_output_M)

# function to generate correct outputs for foldx individual list mutant file
build_output <- function(x){
  # unite into single column
x <- x %>% 
  unite(A_vec,B_vec,C_vec,col = "cmd",sep = ",")

# remove NAs from binding
x <- x %>% mutate(cmd = str_replace_all(cmd, "NA\\,",""))
x <- x %>% mutate(cmd = str_replace_all(cmd, "\\,NA",""))

# add end line ;
x$temp <- ";"

x <- x %>% 
  unite(cmd,temp,col="cmd",sep="")

x <- paste((x),sep="")


# various cleanups of the paste using regex
x <- str_replace(x,"^c","")
x <- str_replace(x,"\\)$","")


x <- str_replace_all(x,"\"+?","")
x <- str_replace_all(x,"\\(+?","")
x <- str_replace_all(x,"\\n","")
x <- str_replace_all(x,"\\s+","")
# remove end line commas and replace with a new line
x <- str_replace_all(x,"\\;\\,+?","\\;\\\n")

# print individual files for each substitution
  cat(paste((x),sep=""),file=paste("individual_list_",str_sub(x,-2,-2),".txt",sep=""),sep="
")
}

final_out_list <- sapply(scan_list, build_output)

# to build a single individual_list containing all mutations instead of separated
# by substitution, remove the cat command inside the function and run following 2 lines instead

# final_out_list <- sapply(scan_list, build_output)

# cat(paste((c(final_out_list[1],final_out_list[2],final_out_list[3],final_out_list[4],final_out_list[5],final_out_list[6],final_out_list[7],final_out_list[8],final_out_list[9],final_out_list[10],final_out_list[11],final_out_list[12],final_out_list[13],final_out_list[14],final_out_list[15],final_out_list[16],final_out_list[17],final_out_list[18],final_out_list[19],final_out_list[20])),sep=""),file="individual_list_allres.txt",sep="
# ")


letter_list <- c("A","Y","T","N","S","F","R","G","V","P","D","K","L","H","Q","W","I","E","C","M")

job_writer <- function(x){
  cat(paste
    ("#!/bin/bash --login

#$ -cwd
#$ -pe smp.pe 4

chmod u+x ./foldx_20231231
./foldx_20231231 --command=BuildModel --pdb=./6vxx_Repair.pdb --mutant-file=./individual_list_",x,".txt --out-pdb=false --output-dir=./buildmodel_output_",x
      ,sep=""),
    file=paste("foldx5_buildmodel_jobscript_",x,".job",sep=""),sep="")
}

lapply(letter_list,job_writer)

namelist <- c("buildmodel_output_A",
                  "buildmodel_output_Y",
                  "buildmodel_output_T",
                  "buildmodel_output_N",
                  "buildmodel_output_S",
                  "buildmodel_output_F",
                  "buildmodel_output_R",
                  "buildmodel_output_G",
                  "buildmodel_output_V",
                  "buildmodel_output_P",
                  "buildmodel_output_D",
                  "buildmodel_output_K",
                  "buildmodel_output_L",
                  "buildmodel_output_H",
                  "buildmodel_output_Q",
                  "buildmodel_output_W",
                  "buildmodel_output_I",
                  "buildmodel_output_E",
                  "buildmodel_output_C",
                  "buildmodel_output_M")

# print list of directory names 1 per line for bash scripting
cat(namelist,file="dir_names.txt",sep="
")
