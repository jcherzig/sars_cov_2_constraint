library(tidyverse)
library(viridis)
library(writexl)
library(readxl)
library(scales)

setwd("/Users/user/Documents/Paper\ planning/results+figures/Figure_1")

# make a list of all VOC subsittutions from Carabelli et al.

VOC_subs <- c("T19I","T19R","A67V","V83A","T95I","G142D","H146Q","Q183E",
              "V213E","V213G","G339D","G339H","R356T","L368I","S371F","S371L",
              "S373P","S375F","T376A","D405N","R408S","K417N","N440K","K444T",
              "V445P","G446S","L452R","N460K","S477N","T478K","E484A","F486V",
              "F490S","Q493R","G496S","Q498R","N501Y","Y505H","T547K","A570D",
              "D614G","H655Y","N679K","P681H","P681R","T716I","N764K","D796Y",
              "N856K","D950N","Q954H","N969K","L981F","S982A","D1118H")

# MAC path
setwd("./data_sets")


# read in WT master table of scores
WT_master <- read_xlsx("6vxx_spike_master_scores_ML_DEC_20.xlsx")

# create arbitrary all positive foldx and rosetta scores, maintaining rank order
WT_master$ddG_arb <- WT_master$DeltaDeltaG + 
  abs(floor(min(WT_master$DeltaDeltaG,na.rm=T)))

WT_master$dR_arb <- abs(WT_master$total_score - 
                            (max(WT_master$total_score,na.rm=T)))

# classify according to mean frequency
quantiles <- quantile(WT_master$Observed_probability)

mean_classifier <- quantiles[3]

WT_master <- WT_master %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

WT_master <- WT_master %>% 
  filter(Chain == "A")

WT_VOC_subs_tbl <- WT_master %>% 
  filter(Replacement %in% VOC_subs)

# read in alpha master table of scores
alpha_master <- read_xlsx("7lws_spike_master_scores_ML_JUN_21.xlsx")

# create arbitrary all positive foldx and rosetta scores, maintaining rank order
alpha_master$ddG_arb <- alpha_master$DeltaDeltaG + 
  abs(floor(min(alpha_master$DeltaDeltaG,na.rm=T)))

alpha_master$dR_arb <- abs(alpha_master$total_score - 
                          (max(alpha_master$total_score,na.rm=T)))

# classify according to mean frequency
quantiles <- quantile(alpha_master$Observed_probability)

mean_classifier <- quantiles[3]

alpha_master <- alpha_master %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

alpha_master <- alpha_master %>% 
  filter(Chain == "A")

alpha_VOC_subs_tbl <- alpha_master %>% 
  filter(Replacement %in% VOC_subs)


# read in delta master table of scores
delta_master <- read_xlsx("7v7n_spike_master_scores_ML_DEC_21.xlsx")

# create arbitrary all positive foldx and rosetta scores, maintaining rank order
delta_master$ddG_arb <- delta_master$DeltaDeltaG + 
  abs(floor(min(delta_master$DeltaDeltaG,na.rm=T)))

delta_master$dR_arb <- abs(delta_master$total_score - 
                          (max(delta_master$total_score,na.rm=T)))

# classify according to mean frequency
quantiles <- quantile(delta_master$Observed_probability)

mean_classifier <- quantiles[3]

delta_master <- delta_master %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

delta_master <- delta_master %>% 
  filter(Chain == "A")

delta_VOC_subs_tbl <- delta_master %>% 
  filter(Replacement %in% VOC_subs)


# read in omicron master table of scores
omicron_master <- read_xlsx("7wp9_spike_master_scores_ML_JUN_22.xlsx")

# create arbitrary all positive foldx and rosetta scores, maintaining rank order
omicron_master$ddG_arb <- omicron_master$DeltaDeltaG + 
  abs(floor(min(omicron_master$DeltaDeltaG,na.rm=T)))

omicron_master$dR_arb <- abs(omicron_master$total_score - 
                          (max(omicron_master$total_score,na.rm=T)))

# classify according to mean frequency
quantiles <- quantile(omicron_master$Observed_probability)

mean_classifier <- quantiles[3]

omicron_master <- omicron_master %>% 
  mutate(Observed_class = case_when(Observed_probability <=mean_classifier ~ "Forbidden", Observed_probability>mean_classifier ~ "Permitted")) 

omicron_master <- omicron_master %>% 
  filter(Chain == "A")

omicron_VOC_subs_tbl <- omicron_master %>% 
  filter(Replacement %in% VOC_subs)

WT_master$Structure <- "WT"
WT_VOC_subs_tbl$Structure <- "WT"

alpha_master$Structure <- "Alpha"
alpha_VOC_subs_tbl$Structure <- "Alpha"

delta_master$Structure <- "Delta"
delta_VOC_subs_tbl$Structure <- "Delta"

omicron_master$Structure <- "Omicron"
omicron_VOC_subs_tbl$Structure <- "Omicron"


full_master <- rbind(WT_master,alpha_master,delta_master,omicron_master)
VOC_subs_tbl_master <- rbind(WT_VOC_subs_tbl,alpha_VOC_subs_tbl,delta_VOC_subs_tbl,omicron_VOC_subs_tbl)

full_master$Structure <- factor(full_master$Structure,levels=c("WT","Alpha","Delta","Omicron"))
VOC_subs_tbl_master$Structure <- factor(VOC_subs_tbl_master$Structure,levels=c("WT","Alpha","Delta","Omicron"))
#

### these plots show the proportional substitution frequency rather than the binary classification result ###
# # log likelihood
# p <- ggplot(full_master,aes(x=Resno,y=LogLikelihood)) +
#   geom_point(aes(col=log10(Observed_probability)),alpha=0.2,position = position_jitter(width = 0.2, height = 0.2),size=0.9)+
#   geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=LogLikelihood,fill=log10(Observed_probability)),position = position_jitter(width = 0.2, height = 0.2),size=1.5,shape=21) +
#   theme_bw(base_size=12)+
#   scale_color_distiller(palette="Blues",direction=1)+
#   scale_fill_distiller(palette="Blues",direction=1)+
#   scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
#   xlab("Position") +
#   ylab("PSSM LogLikelihood") +
#   guides(fill = "none",color="none") +
#   facet_wrap(~Structure)
# p
# 
# ggsave("PSSM_over_protein_with_VOCs.png",width=1080,height=1080,units="px")
# 
# 
# p <- ggplot(full_master,aes(x=Resno,y=ddG_arb)) +
#   geom_point(aes(col=log10(Observed_probability)),alpha=0.2,size=0.9)+
#   geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=ddG_arb,fill=log10(Observed_probability)),size=1.5,shape=21) +
#   theme_bw(base_size=12)+
#   scale_color_distiller(palette="Blues",direction=1)+
#   scale_fill_distiller(palette="Blues",direction=1)+
#   scale_y_log10(labels = label_comma())+
#   scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
#   xlab("Position") +
#   ylab("FoldX \u0394\u0394G") +
#   guides(fill = "none",color="none") +
#   facet_wrap(~Structure)
# p
# 
# ggsave("foldx_over_protein_with_VOCs.png",width=1080,height=1080,units="px")
# 
# p <- ggplot(full_master,aes(x=Resno,y=dR_arb)) +
#   geom_point(aes(col=log10(Observed_probability)),alpha=0.2,size=0.9)+
#   geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=dR_arb,fill=log10(Observed_probability)),size=1.5,shape=21) +
#   theme_bw(base_size=12)+
#   scale_color_distiller(palette="Blues",direction=1)+
#   scale_fill_distiller(palette="Blues",direction=1)+
#   # scale_y_log10(labels = label_comma())+
#   scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
#   xlab("Position") +
#   ylab("Rosetta \u0394R") +
#   guides(fill = "none",color="none") +
#   facet_wrap(~Structure)
# p
# 
# ggsave("rosetta_over_protein_with_VOCs.png",width=1110,height=1080,units="px")
# 
# p <- ggplot(full_master,aes(x=Resno,y=ESST_probability)) +
#   geom_point(aes(col=log10(Observed_probability)),alpha=0.2,size=0.9)+
#   geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=ESST_probability,fill=log10(Observed_probability)),size=1.5,shape=21) +
#   theme_bw(base_size=12)+
#   scale_color_distiller(palette="Blues",direction=1)+
#   scale_fill_distiller(palette="Blues",direction=1)+
#   scale_y_log10(labels = label_comma())+
#   scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
#   xlab("Position") +
#   ylab("ESST probability") +
#   guides(fill = "none",color="none") +
#   facet_wrap(~Structure)
# p
# 
# ggsave("esst_over_protein_with_VOCs.png",width=1160,height=1080,units="px")
# 
# 
# # statistical analysis to compare VOC substitutions to all substitutions and only permitted substitutions
# wilcox.test(master_conservation$ddG_buildmodel_arb,VOC_subs_tbl$ddG_buildmodel_arb,digits.rank = 4,p.adjust.method = "BH",alternative="greater")
# wilcox.test(master_conservation$LogLikelihood,VOC_subs_tbl$LogLikelihood,digits.rank = 4,p.adjust.method = "BH",alternative="less")
# wilcox.test(master_conservation$ESST_probability,VOC_subs_tbl$ESST_probability,digits.rank = 4,p.adjust.method = "BH",alternative="less")
# wilcox.test(master_conservation$deltaRosetta_soft_min_arb,VOC_subs_tbl$deltaRosetta_soft_min_arb,digits.rank = 4,p.adjust.method = "BH",alternative="less")
# 
# wilcox.test(master_conservation_permitted$ddG_buildmodel_arb,VOC_subs_tbl$ddG_buildmodel_arb,digits.rank = 4,p.adjust.method = "BH")
# wilcox.test(master_conservation_permitted$LogLikelihood,VOC_subs_tbl$LogLikelihood,digits.rank = 4,p.adjust.method = "BH")
# wilcox.test(master_conservation_permitted$ESST_probability,VOC_subs_tbl$ESST_probability,digits.rank = 4,p.adjust.method = "BH")
# wilcox.test(master_conservation_permitted$deltaRosetta_soft_min_arb,VOC_subs_tbl$deltaRosetta_soft_min_arb,digits.rank = 4,p.adjust.method = "BH")



# log likelihood
p <- ggplot(full_master,aes(x=Resno,y=LogLikelihood)) +
  geom_point(aes(col=Observed_class),alpha=0.2,position = position_jitter(width = 0.2, height = 0.2),size=0.9)+
  geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=LogLikelihood,fill=Observed_class),position = position_jitter(width = 0.2, height = 0.2),size=1.7,shape=21) +
  theme_bw(base_size=12)+
  scale_colour_manual(values=c("firebrick2","steelblue2"))+
  scale_fill_manual(values=c("firebrick2","steelblue2"))+
  scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
  xlab("Position") +
  ylab("PSSM LogLikelihood") +
  # guides(fill=guide_legend(title="Observed class"),col=guide_legend(title="Observed class"))+
  facet_wrap(~Structure,ncol=1)
p

ggsave("PSSM_over_protein_class_wVOC.png",width=1080,height=1440,units="px")

p <- ggplot(full_master,aes(x=Resno,y=ddG_arb)) +
  geom_point(aes(col=Observed_class),alpha=0.2,size=0.9)+
  geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=ddG_arb,fill=Observed_class),size=1.7,shape=21) +
  theme_bw(base_size=12)+
  scale_colour_manual(values=c("firebrick2","steelblue2"))+
  scale_fill_manual(values=c("firebrick2","steelblue2"))+
  scale_y_log10(labels = label_comma())+
  scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
  xlab("Position") +
  ylab("FoldX \u0394\u0394G") +
  guides(fill = "none",color="none") +
  facet_wrap(~Structure,ncol=1)
p

ggsave("foldx_over_protein_class_wVOC.png",width=1080,height=1440,units="px")


p <- ggplot(full_master,aes(x=Resno,y=dR_arb)) +
  geom_point(aes(col=Observed_class),alpha=0.2,size=0.9)+
  geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=dR_arb,fill=Observed_class),size=1.7,shape=21) +
  theme_bw(base_size=12)+
  scale_colour_manual(values=c("firebrick2","steelblue2"))+
  scale_fill_manual(values=c("firebrick2","steelblue2"))+
  # scale_y_log10(labels = label_comma())+
  scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
  xlab("Position") +
  ylab("Rosetta \u0394R") +
  guides(fill = "none",color="none") +
  facet_wrap(~Structure,ncol=1)
p

ggsave("rosetta_over_protein_class_wVOC.png",width=1080,height=1440,units="px")


p <- ggplot(full_master,aes(x=Resno,y=ESST_probability)) +
  geom_point(aes(col=Observed_class),alpha=0.2,size=0.9)+
  geom_point(data=VOC_subs_tbl_master,aes(x=Resno,y=ESST_probability,fill=Observed_class),size=1.7,shape=21) +
  theme_bw(base_size=12)+
  scale_colour_manual(values=c("firebrick2","steelblue2"))+
  scale_fill_manual(values=c("firebrick2","steelblue2"))+
  # scale_y_log10(labels = label_number())+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100),labels=paste0(c(0.01,0.1,1,10,100)))+
  scale_x_continuous(limits = c(0,1200),breaks=c(0,200,400,600,800,1000,1200))+
  xlab("Position") +
  ylab("ESST probability") +
  guides(fill = "none",color="none") +
  facet_wrap(~Structure,ncol=1)
p

ggsave("esst_over_protein_class_wVOC.png",width=1080,height=1440,units="px")
