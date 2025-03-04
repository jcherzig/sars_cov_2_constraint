library(tidyverse)
library(viridis)
library(writexl)
library(readxl)
library(scales)

# make a list of all VOC subsittutions from Carabelli et al.

VOC_subs <- c("T19I","T19R","A67V","V83A","T95I","G142D","H146Q","Q183E",
              "V213E","V213G","G339D","G339H","R356T","L368I","S371F","S371L",
              "S373P","S375F","T376A","D405N","R408S","K417N","N440K","K444T",
              "V445P","G446S","L452R","N460K","S477N","T478K","E484A","F486V",
              "F490S","Q493R","G496S","Q498R","N501Y","Y505H","T547K","A570D",
              "D614G","H655Y","N679K","P681H","P681R","T716I","N764K","D796Y",
              "N856K","D950N","Q954H","N969K","L981F","S982A","D1118H")

VOCs_through_time <- function(date){
  
  setwd(paste0("/Users/user/Documents/Case study paper/3mo_loos_ensemble/case_study_3moloos_unaligned_50/spike_filtered_",date,"/Master_tables/6vxx"))
  
  # check if this date includes 12 month test set
  test_sets <- list.files(pattern="*TEST_month12_*")
  if(length(test_sets)==0){
    full_tests<-F
  } else {full_tests<-T}
  
  
  if(full_tests==T){
  # read in 1 month test set
  # read in master table of scores
  master_conservation_1 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month1_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_1$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_1_50 <- master_conservation_1 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  master_conservation_1_50$Period <- "1 month"
  
  rm(master_conservation_1)
  
  # read in 3 month test set
  # read in master table of scores
  master_conservation_3 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month3_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_3$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_3_50 <- master_conservation_3 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  master_conservation_3_50$Period <- "3 month"
  
  rm(master_conservation_3)
  
  # read in 6 month test set
  # read in master table of scores
  master_conservation_6 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month6_",date,".xlsx"))
  
  # calculate quantiles of observed probability to defin classes later
  quantiles <- quantile(master_conservation_6$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_6_50 <- master_conservation_6 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
  
  master_conservation_6_50$Period <- "6 month"
  
  rm(master_conservation_6)
  
  # read in 12 month test set
  # read in master table of scores
  master_conservation_12 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month12_",date,".xlsx"))
  
  # calculate quantiles of observed probability to define classes later
  quantiles <- quantile(master_conservation_12$Observed_probability)
  
  mean_classifier <- quantiles[3]
  
  master_conservation_12_50 <- master_conservation_12 %>% 
    mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted"))
  
  master_conservation_12_50$Period <- "12 month"
  
  rm(master_conservation_12)
  
  
  # add them all into one tibble
  master_conservation_full <- rbind(master_conservation_1_50,master_conservation_3_50,master_conservation_6_50,master_conservation_12_50)
  
  master_conservation_full <- master_conservation_full %>% 
    filter(Chain == "A")
  
  # filter for only VOCs 
  master_conservation_full <- master_conservation_full %>% 
    filter(Replacement %in% VOC_subs) %>% 
    dplyr::select(Replacement,Observed_class,Period)

  master_conservation_full$Date <- date
  
  return(master_conservation_full)  
  } else {
    # read in 1 month test set
    # read in master table of scores
    master_conservation_1 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month1_",date,".xlsx"))
    
    # calculate quantiles of observed probability to defin classes later
    quantiles <- quantile(master_conservation_1$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_1_50 <- master_conservation_1 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    master_conservation_1_50$Period <- "1 month"
    
    rm(master_conservation_1)
    
    # read in 3 month test set
    # read in master table of scores
    master_conservation_3 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month3_",date,".xlsx"))
    
    # calculate quantiles of observed probability to defin classes later
    quantiles <- quantile(master_conservation_3$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_3_50 <- master_conservation_3 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    master_conservation_3_50$Period <- "3 month"
    
    rm(master_conservation_3)
    
    # read in 6 month test set
    # read in master table of scores
    master_conservation_6 <- read_xlsx(paste0("6vxx_spike_master_scores_ML_TEST_month6_",date,".xlsx"))
    
    # calculate quantiles of observed probability to define classes later
    quantiles <- quantile(master_conservation_6$Observed_probability)
    
    mean_classifier <- quantiles[3]
    
    master_conservation_6_50 <- master_conservation_6 %>% 
      mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 
    
    master_conservation_6_50$Period <- "6 month"
    
    rm(master_conservation_6)

    # add them all into one tibble
    master_conservation_full <- rbind(master_conservation_1_50,master_conservation_3_50,master_conservation_6_50)
    
    master_conservation_full <- master_conservation_full %>% 
      filter(Chain == "A")
    
    # filter for only VOCs 
    master_conservation_full <- master_conservation_full %>% 
      filter(Replacement %in% VOC_subs) %>% 
      dplyr::select(Replacement,Observed_class,Period)
    
    master_conservation_full$Date <- date
    
    return(master_conservation_full) 
  }
}


date_list <- c("DEC_20","JAN_21","FEB_21","MAR_21","APR_21","MAY_21","JUN_21","JUL_21","AUG_21","SEP_21","OCT_21","NOV_21","DEC_21",
               "JAN_22","FEB_22","MAR_22","APR_22","MAY_22","JUN_22","JUL_22","AUG_22","SEP_22","OCT_22","NOV_22","DEC_22")

all_out <- do.call(rbind,lapply(date_list,VOCs_through_time))

all_out <- all_out %>% 
  mutate(Period = case_when(Period == "1 month" ~ "1-month",
                            Period == "3 month" ~ "3-month",
                            Period == "6 month" ~ "6-month",
                            Period == "12 month" ~ "12-month"))

all_out$Period <- factor(all_out$Period,levels=c("1-month","3-month","6-month","12-month"))

all_out$Date <- factor(all_out$Date,levels=c("DEC_20","JAN_21","FEB_21","MAR_21","APR_21","MAY_21","JUN_21","JUL_21","AUG_21","SEP_21","OCT_21","NOV_21","DEC_21",
                                             "JAN_22","FEB_22","MAR_22","APR_22","MAY_22","JUN_22","JUL_22","AUG_22","SEP_22","OCT_22","NOV_22","DEC_22"))

all_out$Replacement <- factor(all_out$Replacement,levels=rev(c("T19I","T19R","A67V","V83A","T95I","G142D","H146Q","Q183E",
                                                           "V213E","V213G","G339D","G339H","R356T","L368I","S371F","S371L",
                                                           "S373P","S375F","T376A","D405N","R408S","K417N","N440K","K444T",
                                                           "V445P","G446S","L452R","N460K","S477N","T478K","E484A","F486V",
                                                           "F490S","Q493R","G496S","Q498R","N501Y","Y505H","T547K","A570D",
                                                           "D614G","H655Y","N679K","P681H","P681R","T716I","N764K","D796Y",
                                                           "N856K","D950N","Q954H","N969K","L981F","S982A","D1118H")))

all_out <- all_out %>% 
  mutate(month_nums = case_when(Date == "DEC_20" ~ 0,
                                Date == "JAN_21" ~ 1,
                                Date == "FEB_21" ~ 2,
                                Date == "MAR_21" ~ 3,
                                Date == "APR_21" ~ 4,
                                Date == "MAY_21" ~ 5,
                                Date == "JUN_21" ~ 6,
                                Date == "JUL_21" ~ 7,
                                Date == "AUG_21" ~ 8,
                                Date == "SEP_21" ~ 9,
                                Date == "OCT_21" ~ 10,
                                Date == "NOV_21" ~ 11,
                                Date == "DEC_21" ~ 12,
                                Date == "JAN_22" ~ 13,
                                Date == "FEB_22" ~ 14,
                                Date == "MAR_22" ~ 15,
                                Date == "APR_22" ~ 16,
                                Date == "MAY_22" ~ 17,
                                Date == "JUN_22" ~ 18,
                                Date == "JUL_22" ~ 19,
                                Date == "AUG_22" ~ 20,
                                Date == "SEP_22" ~ 21,
                                Date == "OCT_22" ~ 22,
                                Date == "NOV_22" ~ 23,
                                Date == "DEC_22" ~ 24
  ))

all_out <- all_out %>% 
  mutate(Date = case_when(Date == "DEC_20" ~ "12/20",
                          Date == "JAN_21" ~ "01/21",
                          Date == "FEB_21" ~ "02/21",
                          Date == "MAR_21" ~ "03/21",
                          Date == "APR_21" ~ "04/21",
                          Date == "MAY_21" ~ "05/21",
                          Date == "JUN_21" ~ "06/21",
                          Date == "JUL_21" ~ "07/21",
                          Date == "AUG_21" ~ "08/21",
                          Date == "SEP_21" ~ "09/21",
                          Date == "OCT_21" ~ "10/21",
                          Date == "NOV_21" ~ "11/21",
                          Date == "DEC_21" ~ "12/21",
                          Date == "JAN_22" ~ "01/22",
                          Date == "FEB_22" ~ "02/22",
                          Date == "MAR_22" ~ "03/22",
                          Date == "APR_22" ~ "04/22",
                          Date == "MAY_22" ~ "05/22",
                          Date == "JUN_22" ~ "06/22",
                          Date == "JUL_22" ~ "07/22",
                          Date == "AUG_22" ~ "08/22",
                          Date == "SEP_22" ~ "09/22",
                          Date == "OCT_22" ~ "10/22",
                          Date == "NOV_22" ~ "11/22",
                          Date == "DEC_22" ~ "12/22"
  ))

axis_labs <- c("12/20","03/21","06/21","09/21","12/21","03/22","06/22","09/22","12/22")
axis_points <- c(0,3,6,9,12,15,18,21,24)

p <- ggplot(all_out,aes(x = month_nums,y=Replacement,fill=Observed_class)) +
  geom_tile(colour="black")+
  theme_bw(base_size=12)+
  xlab("Training cutoff date") +
  ylab("VOC substitution")+
  scale_x_continuous(labels=axis_labs,breaks=axis_points) +
  scale_fill_manual(values=c("firebrick3","steelblue3"))+
  labs(fill="Observed class")+
  facet_wrap(~Period)
p

ggsave("VOC_classification_heatmap.png",height=4500,width=3000,units="px")
