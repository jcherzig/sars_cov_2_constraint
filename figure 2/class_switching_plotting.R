library(readxl)
library(tidyverse)
library(viridis)
library(bio3d)
library(writexl)

# setwd to print results
setwd("/Users/user/Documents/Paper\ planning/results+figures/class_switch_data")

filelist_wt <- list.files(pattern="^6vxx*")
filelist_alpha <- list.files(pattern="^7lws*")
filelist_delta <- list.files(pattern="^7v7n*")
filelist_omicron <- list.files(pattern="^7wp9*")

wt_cs <- do.call(rbind,lapply(filelist_wt,read_xlsx))
alpha_cs <- do.call(rbind,lapply(filelist_alpha,read_xlsx))
delta_cs <- do.call(rbind,lapply(filelist_delta,read_xlsx))
omicron_cs <- do.call(rbind,lapply(filelist_omicron,read_xlsx))

wt_cs$Variant <- "WT"
alpha_cs$Variant <- "Alpha"
delta_cs$Variant <- "Delta"
omicron_cs$Variant <- "Omicron"

all_cs <- rbind(wt_cs,alpha_cs,delta_cs,omicron_cs)

# first calculate the predicted class
all_cs <- all_cs %>% 
  mutate(Predicted_class = case_when(Forbidden>0.5 ~ "Forbidden",Forbidden<0.5 ~ "Permitted",Forbidden==0.5 ~ "Split")) 

all_cs$Resno <- as.numeric(all_cs$Resno)

all_cs <- all_cs %>% 
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

all_cs$Date <- factor(all_cs$Date,levels=c("12/20","01/21","02/21","03/21","04/21","05/21","06/21","07/21","08/21","09/21","10/21","11/21","12/21",
                                           "01/22","02/22","03/22","04/22","05/22","06/22","07/22","08/22","09/22","10/22","11/22","12/22"))




# filter for only incorrect predictions
all_cs_incorrect <- all_cs %>% 
  filter(Observed_class != Predicted_class)

# filter for only correct predictions
all_cs_correct <- all_cs %>% 
  filter(Observed_class == Predicted_class)

# all_cs$Variant <- factor(all_cs$Variant,levels=c("WT","Alpha","Delta","Omicron"))

palette=viridis(n=6,option="G")

p <- ggplot(all_cs_incorrect %>% filter(Group=="3month"),aes(x=Resno,colour=Variant,fill=Variant))+
  geom_density(alpha=0.7) +  
  geom_vline(xintercept = 331, linetype=2) +
  geom_vline(xintercept = 528, linetype=2) +
  xlab("Position") +
  ylab("Density") +
  scale_colour_manual(values=palette[2:5]) +
  scale_fill_manual(values=palette[2:5]) +
  # scale_y_continuous(labels = label_comma())+
  scale_y_continuous(limits=c(0,0.0020),breaks=c(0,0.0010))+
  # scale_x_log10()+
  theme_bw(base_size=12) +
  facet_wrap(~Date,ncol = 1)+
  theme(strip.background = element_blank(),
        strip.text=element_blank())
p

ggsave("incorrect_cs_preds_density_plot_3month.png",width=2000,height=4000,units="px")

p <- ggplot(all_cs_correct %>% filter(Group=="3month"),aes(x=Resno,colour=Variant,fill=Variant))+
  geom_density(alpha=0.7) +  
  geom_vline(xintercept = 331, linetype=2) +
  geom_vline(xintercept = 528, linetype=2) +
  xlab("Position") +
  ylab("Density") +
  scale_colour_manual(values=palette[2:5]) +
  scale_fill_manual(values=palette[2:5]) +
  # scale_y_continuous(labels = label_comma())+
  scale_y_continuous(limits=c(0,0.0020),breaks=c(0,0.0010))+
  # scale_x_log10()+
  theme_bw(base_size=12) +
  facet_wrap(~Date,ncol = 1) +
  theme(strip.background = element_blank(),
        strip.text=element_blank())
p

ggsave("correct_cs_preds_density_plot_3month.png",width=2000,height=4000,units="px")

