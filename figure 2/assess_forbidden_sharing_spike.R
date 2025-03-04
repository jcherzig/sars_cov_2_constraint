library(readxl)
library(tidyverse)
library(viridis)
library(bio3d)
library(writexl)


shared_forbidden_list <- function(structure,date){
  
  tryCatch({setwd(paste0("/Users/user/Documents/Case\ study\ paper/3mo_loos_ensemble/case_study_3moloos_unaligned_50/spike_filtered_",date,"/Master_tables/",structure))
    # read in predictions made by ensemble model
    alltrain_count_probs <- read.csv("rf_50_ensemble_count_probs.csv")
    alltrain_count_probs$Structure <- structure
    alltrain_count_probs$Date <- date
    
    return(alltrain_count_probs)
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
shared_list <- cmapply(shared_forbidden_list,structure="6vxx",date_list)

# shared_list = shared_list[-which(sapply(shared_list, is.null))]

indiv_test_shared <- function(x){
  
  out_list <- list()
  
  for(i in 1:length(shared_list)){
    
    set1 <- shared_list[[x]] %>% 
      group_by(Group) %>% 
      group_split()
    
    set2 <- shared_list[[i]] %>% 
      group_by(Group) %>% 
      group_split()
    
    if(length(set1) == 4 & length(set2) == 4){
      
      test_vec1 <- set1[[1]]$Forbidden == set2[[1]]$Forbidden
      test_vec2 <- set1[[2]]$Forbidden == set2[[2]]$Forbidden
      test_vec3 <- set1[[3]]$Forbidden == set2[[3]]$Forbidden
      test_vec4 <- set1[[4]]$Forbidden == set2[[4]]$Forbidden
      
      sum1 <- sum(test_vec1, na.rm=TRUE)
      sum2 <- sum(test_vec2, na.rm=TRUE)
      sum3 <- sum(test_vec3, na.rm=TRUE)
      sum4 <- sum(test_vec4, na.rm=TRUE)
      
      sum1 <- sum1/nrow(set1[[1]])
      sum2 <- sum2/nrow(set1[[2]])
      sum3 <- sum3/nrow(set1[[3]])
      sum4 <- sum4/nrow(set1[[4]])
      
      
      output <- tibble(shared_prop = c(sum1,sum2,sum3,sum4),
                       group = c(unique(set1[[1]]$Group),unique(set1[[2]]$Group),unique(set1[[3]]$Group),unique(set1[[4]]$Group)),
                       date = c(unique(set1[[1]]$Date),unique(set1[[2]]$Date),unique(set1[[3]]$Date),unique(set1[[4]]$Date)),
                       comparison = c(unique(set2[[1]]$Date),unique(set2[[2]]$Date),unique(set2[[3]]$Date),unique(set2[[4]]$Date)))
      
    } else if(length(set1) == 4 & length(set2) == 3) {
      test_vec1 <- set1[[2]]$Forbidden == set2[[1]]$Forbidden
      test_vec2 <- set1[[3]]$Forbidden == set2[[2]]$Forbidden
      test_vec3 <- set1[[4]]$Forbidden == set2[[3]]$Forbidden
      
      sum1 <- sum(test_vec1, na.rm=TRUE)
      sum2 <- sum(test_vec2, na.rm=TRUE)
      sum3 <- sum(test_vec3, na.rm=TRUE)
      
      sum1 <- sum1/nrow(set1[[2]])
      sum2 <- sum2/nrow(set1[[3]])
      sum3 <- sum3/nrow(set1[[4]])
      
      output <- tibble(shared_prop = c(sum1,sum2,sum3),
                       group = c(unique(set1[[2]]$Group),unique(set1[[3]]$Group),unique(set1[[4]]$Group)),
                       date = c(unique(set1[[2]]$Date),unique(set1[[3]]$Date),unique(set1[[4]]$Date)),
                       comparison = c(unique(set2[[1]]$Date),unique(set2[[2]]$Date),unique(set2[[3]]$Date)))
      
    } else if(length(set1) == 3 & length(set2) == 4) {
      test_vec1 <- set1[[1]]$Forbidden == set2[[2]]$Forbidden
      test_vec2 <- set1[[2]]$Forbidden == set2[[3]]$Forbidden
      test_vec3 <- set1[[3]]$Forbidden == set2[[4]]$Forbidden
      
      sum1 <- sum(test_vec1, na.rm=TRUE)
      sum2 <- sum(test_vec2, na.rm=TRUE)
      sum3 <- sum(test_vec3, na.rm=TRUE)
      
      sum1 <- sum1/nrow(set1[[1]])
      sum2 <- sum2/nrow(set1[[2]])
      sum3 <- sum3/nrow(set1[[3]])
      
      output <- tibble(shared_prop = c(sum1,sum2,sum3),
                       group = c(unique(set1[[1]]$Group),unique(set1[[2]]$Group),unique(set1[[3]]$Group)),
                       date = c(unique(set1[[1]]$Date),unique(set1[[2]]$Date),unique(set1[[3]]$Date)),
                       comparison = c(unique(set2[[2]]$Date),unique(set2[[3]]$Date),unique(set2[[4]]$Date)))
      
    } else if (length(set1) == 3 & length(set2) == 3){
      test_vec1 <- set1[[1]]$Forbidden == set2[[1]]$Forbidden
      test_vec2 <- set1[[2]]$Forbidden == set2[[2]]$Forbidden
      test_vec3 <- set1[[3]]$Forbidden == set2[[3]]$Forbidden
      
      sum1 <- sum(test_vec1, na.rm=TRUE)
      sum2 <- sum(test_vec2, na.rm=TRUE)
      sum3 <- sum(test_vec3, na.rm=TRUE)
      
      sum1 <- sum1/nrow(set1[[1]])
      sum2 <- sum2/nrow(set1[[2]])
      sum3 <- sum3/nrow(set1[[3]])
      
      output <- tibble(shared_prop = c(sum1,sum2,sum3),
                       group = c(unique(set1[[1]]$Group),unique(set1[[2]]$Group),unique(set1[[3]]$Group)),
                       date = c(unique(set1[[1]]$Date),unique(set1[[2]]$Date),unique(set1[[3]]$Date)),
                       comparison = c(unique(set2[[1]]$Date),unique(set2[[2]]$Date),unique(set2[[3]]$Date)))
    }
    out_list[[i]] <- output
  }
  out_tbl <- do.call(rbind,out_list)
  return(out_tbl)
}


x <- 1:length(shared_list)

full_out <- lapply(x,indiv_test_shared)

full_out_tbl <- do.call(rbind,full_out)

full_out_tbl <- full_out_tbl %>% 
  mutate(date = case_when(date == "DEC_20" ~ "12/20",
                          date == "JAN_21" ~ "01/21",
                          date == "FEB_21" ~ "02/21",
                          date == "MAR_21" ~ "03/21",
                          date == "APR_21" ~ "04/21",
                          date == "MAY_21" ~ "05/21",
                          date == "JUN_21" ~ "06/21",
                          date == "JUL_21" ~ "07/21",
                          date == "AUG_21" ~ "08/21",
                          date == "SEP_21" ~ "09/21",
                          date == "OCT_21" ~ "10/21",
                          date == "NOV_21" ~ "11/21",
                          date == "DEC_21" ~ "12/21",
                          date == "JAN_22" ~ "01/22",
                          date == "FEB_22" ~ "02/22",
                          date == "MAR_22" ~ "03/22",
                          date == "APR_22" ~ "04/22",
                          date == "MAY_22" ~ "05/22",
                          date == "JUN_22" ~ "06/22",
                          date == "JUL_22" ~ "07/22",
                          date == "AUG_22" ~ "08/22",
                          date == "SEP_22" ~ "09/22",
                          date == "OCT_22" ~ "10/22",
                          date == "NOV_22" ~ "11/22",
                          date == "DEC_22" ~ "12/22"
  ))

full_out_tbl <- full_out_tbl %>% 
  mutate(comparison = case_when(comparison == "DEC_20" ~ "12/20",
                          comparison == "JAN_21" ~ "01/21",
                          comparison == "FEB_21" ~ "02/21",
                          comparison == "MAR_21" ~ "03/21",
                          comparison == "APR_21" ~ "04/21",
                          comparison == "MAY_21" ~ "05/21",
                          comparison == "JUN_21" ~ "06/21",
                          comparison == "JUL_21" ~ "07/21",
                          comparison == "AUG_21" ~ "08/21",
                          comparison == "SEP_21" ~ "09/21",
                          comparison == "OCT_21" ~ "10/21",
                          comparison == "NOV_21" ~ "11/21",
                          comparison == "DEC_21" ~ "12/21",
                          comparison == "JAN_22" ~ "01/22",
                          comparison == "FEB_22" ~ "02/22",
                          comparison == "MAR_22" ~ "03/22",
                          comparison == "APR_22" ~ "04/22",
                          comparison == "MAY_22" ~ "05/22",
                          comparison == "JUN_22" ~ "06/22",
                          comparison == "JUL_22" ~ "07/22",
                          comparison == "AUG_22" ~ "08/22",
                          comparison == "SEP_22" ~ "09/22",
                          comparison == "OCT_22" ~ "10/22",
                          comparison == "NOV_22" ~ "11/22",
                          comparison == "DEC_22" ~ "12/22"
  ))

full_out_tbl$date <- factor(full_out_tbl$date,levels=c("12/20","01/21","02/21","03/21","04/21","05/21","06/21","07/21","08/21",
                                                          "09/21","10/21","11/21","12/21","01/22","02/22","03/22","04/22",
                                                          "05/22","06/22","07/22","08/22","09/22","10/22","11/22","12/22"))

full_out_tbl$comparison <- factor(full_out_tbl$comparison,levels=c("12/20","01/21","02/21","03/21","04/21","05/21","06/21","07/21","08/21",
                                                          "09/21","10/21","11/21","12/21","01/22","02/22","03/22","04/22",
                                                          "05/22","06/22","07/22","08/22","09/22","10/22","11/22","12/22"))

full_out_tbl <- full_out_tbl %>% 
  mutate(group = case_when(group == "1month" ~ "1-month",
                           group == "3month" ~ "3-month",
                           group == "6month" ~ "6-month",
                           group == "12month" ~ "12-month"))

full_out_tbl$group <- factor(full_out_tbl$group,levels=c("1-month","3-month","6-month","12-month"))

full_out_tbl_print <- full_out_tbl %>% 
  pivot_wider(id_cols = c(group,comparison),names_from=date,values_from=shared_prop)

indiv_groups <- full_out_tbl_print %>% 
  group_by(group) %>% 
  group_split()

setwd("/Users/user/Documents/Paper\ planning/results+figures/Figure_4")

write_xlsx(indiv_groups[[1]],"forbidden_comp_6vxx_12month.xlsx")
write_xlsx(indiv_groups[[2]],"forbidden_comp_6vxx_1month.xlsx")
write_xlsx(indiv_groups[[3]],"forbidden_comp_6vxx_3month.xlsx")
write_xlsx(indiv_groups[[4]],"forbidden_comp_6vxx_6month.xlsx")


p <- ggplot(full_out_tbl %>% filter(group=="6-month" | group=="1-month"),aes(x=date,y=comparison,fill=shared_prop)) +
  geom_tile()+
  scale_fill_viridis(name = "Shared predictions")+
  xlab("Date")+
  ylab("Comparison") +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  facet_wrap(~group)
p  

ggsave("6vxx_shared_forbidden_heatmap_1_6month.png", width = 3500, height = 1600, unit="px",dpi = 300)
