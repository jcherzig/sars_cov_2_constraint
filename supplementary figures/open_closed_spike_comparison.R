library(tidyverse)
library(viridis)
library(writexl)
library(bio3d)
library(fishualize)

# read in ESST results for 6xm0 and 6xm5
ESST_closed <- read.csv("PROCESSED_crescendo_dumpouts_6xm5A.csv")
ESST_open <- read.csv("PROCESSED_crescendo_dumpouts_6xm0AB.csv")

ESST <- rbind(ESST_closed,ESST_open)

 # read in FoldX results for 6XM0 (open) and 6XM5 (closed)

 raw_6xm5 <- read.csv("processed_foldx_posscan_6xm5.csv")

 raw_6xm0 <- read.csv("processed_foldx_posscan_6xm0.csv")

 raw_6xm5 <- raw_6xm5 %>%
   filter(Chain == "A" | Chain == "B")

 raw_6xm0 <- raw_6xm0 %>%
   filter(Chain == c("A"))

 raw_6xm5 <- raw_6xm5 %>%
   mutate(Structure = case_when(Chain == "A" ~ "Open conformation; 'down' chain",
                                Chain == "B" ~ "Open conformation; 'up' chain"))

raw_6xm0$Structure <- "Closed conformation"

foldx <- rbind(raw_6xm0,raw_6xm5)

# remove one missing residue from the 6XM5 run
foldx <- foldx %>% 
  filter(Resno != 431)

colnames(ESST) <- c("Resno","WTAA","ReplacementAA","Probability","Structure","Chain")

foldx$WTAA <- aa321(foldx$WTAA)

# look at the individual methods against themselves in different conf
# foldx first
foldx_wide <- foldx %>% 
  pivot_wider(id_cols=-Chain,names_from=Structure,values_from=DeltaG)

foldx_wide <- foldx_wide %>%
  mutate(Region = case_when(Resno >= 14 & Resno <= 306 ~ "N terminal",
                            Resno >= 331 & Resno <= 528 ~ "RBD",
                            Resno >= 529 & Resno <= 686 ~ "C terminal",
                            Resno >= 687 ~ "S2",
                            .default="Indeterminate"))

foldx_wide <- foldx_wide %>% 
  filter(WTAA != "H")

foldx_wide$Region <- factor(foldx_wide$Region, levels=c("N terminal","RBD","C terminal","S2","Indeterminate"))


p <- ggplot(foldx_wide,aes(x=`Open conformation; 'down' chain`,y=`Open conformation; 'up' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16)+
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("Open conformation; 'down' chain") +
  # ylab("Open conformation; 'up' chain")
p
ggsave("open_down_vs_open_up_foldx.png",width=1080,height=1080,units="px")

p <- ggplot(foldx_wide,aes(x=`Closed conformation`,y=`Open conformation; 'down' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16)+
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("Closed conformation") +
  # ylab("Open conformation; 'down' chain")
p
ggsave("closed_vs_open_down_foldx.png",width=1080,height=1080,units="px")


p <- ggplot(foldx_wide,aes(x=`Closed conformation`,y=`Open conformation; 'up' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16)+
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("Closed conformation") +
  # ylab("Open conformation; 'up' chain")
p
ggsave("closed_vs_open_up_foldx.png",width=1080,height=1080,units="px")


# then ESST
ESST_wide <- ESST %>% 
  pivot_wider(id_cols=-Chain,names_from=Structure,values_from=Probability)


ESST_wide <- ESST_wide %>%
  mutate(Region = case_when(Resno >= 14 & Resno <= 306 ~ "N terminal",
                            Resno >= 331 & Resno <= 528 ~ "RBD",
                            Resno >= 529 & Resno <= 686 ~ "C terminal",
                            Resno >= 687 ~ "S2",
                            .default="Indeterminate"))

ESST_wide$Region <- factor(ESST_wide$Region, levels=c("N terminal","RBD","C terminal","S2","Indeterminate"))

p <- ggplot(ESST_wide,aes(x=`Open conformation; 'down' chain`,y=`Open conformation; 'up' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16)+
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("Open conformation; 'down' chain") +
  # ylab("Open conformation; 'up' chain")
p

ggsave("open_down_vs_open_up_ESST.png",width=1080,height=1080,units="px")


p <- ggplot(ESST_wide,aes(x=`Closed conformation`,y=`Open conformation; 'down' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16) +
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("Closed conformation") +
  # ylab("Open conformation; 'down' chain")
p

ggsave("closed_vs_open_down_ESST.png",width=1080,height=1080,units="px")


p <- ggplot(ESST_wide,aes(x=`Closed conformation`,y=`Open conformation; 'up' chain`,col=Region))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm") +
  geom_abline(intercept=0,linetype="dashed")+
  scale_color_fish_d(option = "Scarus_quoyi",guide="none") +
  theme_bw(base_size = 16)+
  xlab(NULL) +
  ylab(NULL)
  # scale_x_log10()+
  # xlab("ESST probability") +
  # ylab("ESST probability")
p

ggsave("closed_vs_open_up_ESST.png",width=1080,height=1080,units="px")


# add a dummy var for region
ESST_wide <- ESST_wide %>% 
  mutate(region_var = case_when(Region == "N terminal"~0,
                                Region == "RBD"~1,
                                Region == "C terminal"~2,
                                Region == "S2"~3,
                                Region == "Indeterminate"~4))

# means significance tests between up and down conformation chains
foldx_wide <- foldx_wide %>% 
  dplyr::rename("Closed" = "Closed conformation", "Open_down"="Open conformation; 'down' chain","Open_up"="Open conformation; 'up' chain")

ESST_wide <- ESST_wide %>% 
  dplyr::rename("Closed" = "Closed conformation", "Open_down"="Open conformation; 'down' chain","Open_up"="Open conformation; 'up' chain")


foldx_closed_open_down <- wilcox.test(foldx_wide$`Closed conformation`,foldx_wide$`Open conformation; 'down' chain`,paired=T,p.adjust.method = "BH")
foldx_closed_open_up <- wilcox.test(foldx_wide$`Closed conformation`,foldx_wide$`Open conformation; 'up' chain`,paired=T,p.adjust.method = "BH")
foldx_open_down_open_up <- wilcox.test(foldx_wide$`Open conformation; 'down' chain`,foldx_wide$`Open conformation; 'up' chain`,paired=T,p.adjust.method = "BH")

ESST_closed_open_down <- wilcox.test(ESST_wide$`Closed conformation`,ESST_wide$`Open conformation; 'down' chain`,paired=T,p.adjust.method = "BH")
ESST_closed_open_up <- wilcox.test(ESST_wide$`Closed conformation`,ESST_wide$`Open conformation; 'up' chain`,paired=T,p.adjust.method = "BH")
ESST_open_down_open_up <- wilcox.test(ESST_wide$`Open conformation; 'down' chain`,ESST_wide$`Open conformation; 'up' chain`,paired=T,p.adjust.method = "BH")
