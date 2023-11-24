#############################################################################################
# Analysis of qPCR Data
#############################################################################################

# Date : 21 November 2023
# Author : Ari Fina Bintarti

# Install packages
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("ggtext")
install.packages("sjmisc") 
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages('mvtnorm', dep = TRUE)
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("lmerTest")
install.packages("multcomp")
install.packages("ape")
install.packages("devtools", dependencies = TRUE)
BiocManager::install("phyloseq")
install.packages("datarium")
install.packages("rstatix")
library(devtools)
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggtext)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)
library(ape)
library(phyloseq)
library(datarium)
library(rstatix)

###########################################################################
# 1. Response variable: AOA COPY NUMBER
###########################################################################
# 1a. Analyses of Bulk Soil
setwd('D:/Fina/INRAE_Project/microservices/')
qPCR <- read.csv("qPCR_results.csv")
qPCR.BS <- qPCR[1:120,]
str(qPCR.BS)
qPCR.BS$nr <- as.factor(qPCR.BS$nr)
qPCR.BS$plot <- as.factor(qPCR.BS$plot)
qPCR.BS$block <- as.factor(qPCR.BS$block)
qPCR.BS$irrigation <- as.factor(qPCR.BS$irrigation)
qPCR.BS$fertilization <- as.factor(qPCR.BS$fertilization)
qPCR.BS$type <- as.factor(qPCR.BS$type)
qPCR.BS$sampling.date <- as.factor(qPCR.BS$sampling.date)
qPCR.BS$var3 <- as.factor(qPCR.BS$var3)
qPCR.BS$x <- as.factor(qPCR.BS$x)
qPCR.BS$sampling.date <- factor(qPCR.BS$sampling.date, levels = c("28/04/2022", "01/06/2022", "05/07/2022", "20/07/2022", "13/09/2022"),
                                labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))


qPCR.BS$AOA_logDWSoil <- log10(qPCR.BS$AOA_nbc_per_g_DW_soil)
qPCR.BS$AOB_logDWSoil <- log10(qPCR.BS$AOB_nbc_per_g_DW_soil)
qPCR.BS$AOA_logngDNA <- log10(qPCR.BS$AOA_nbc_per_ngDNA)
qPCR.BS$AOB_logngDNA <- log10(qPCR.BS$AOB_nbc_per_ngDNA)


aoa.BS.copies <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOA_nbc_per_g_DW_soil, type = "mean_sd")


aoa.BS.copies.plot <- ggboxplot(
  qPCR.BS, x = "irrigation", y = "AOA_nbc_per_g_DW_soil",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aoa.BS.copies.plot

# check assumption (outliers)
aoa.BS.copies.out <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOA_nbc_per_g_DW_soil) # no extreme outliers
# Saphiro-Wilk for normality
aoa.BS.copies.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOA_nbc_per_g_DW_soil)
ggqqplot(qPCR.BS, "AOA_nbc_per_g_DW_soil", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.BS.copies.Lave <- qPCR.BS %>%
  group_by(sampling.date) %>%
  levene_test(AOA_nbc_per_g_DW_soil ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.BS.copies.aov <- anova_test(
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.BS.copies.aov)
# TWO-WAY
res.aov2 <- aov(AOA_nbc_per_g_DW_soil ~ irrigation + fertilization, data = qPCR.BS)
summary(res.aov2)
# THREE WAY
res.aov3 <- qPCR.BS %>% anova_test(AOA_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date)
res.aov3

res.aov3.ngDNA <- qPCR.BS %>% anova_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date)
res.aov3.ngDNA

aoa.BS.ngDNA.aov <- anova_test(
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.BS.ngDNA.aov)
# model
aoa.BS.copies.mod <- lm(qPCR.BS$AOA_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aoa.BS.copies.mod)
# pairwise comparisons
aoa.test.pw <- qPCR.BS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_g_DW_soil ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.BS.copies.mod)
aoa.test.pw
# Pairwise comparisons
library(emmeans)
pwc.aoa.DWS <- qPCR.BS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_g_DW_soil ~ irrigation, p.adjust.method = "BH") %>%
  select(-df, -statistic, -p) # Remove details
pwc.aoa.DWS # either grouped by fertilization or both by fertilization and per sampling date, drought is not significant

##subset
aoa.cop.M <- qPCR.BS[which(qPCR.BS$fertilization == "M"),]
aoa.cop.M
aoa.cop.D <- qPCR.BS[which(qPCR.BS$fertilization == "D"),]
aoa.cop.D
aoa.cop.K <- qPCR.BS[which(qPCR.BS$fertilization == "K"),]
aoa.cop.K
aoa.cop.M.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M)
summary(aoa.cop.M.aov)
TukeyHSD(aoa.cop.M.aov, conf.level=.95)

aoa.cop.D.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.D)
summary(aoa.cop.D.aov)
TukeyHSD(aoa.cop.D.aov, conf.level=.95)

aoa.cop.K.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.K)
summary(aoa.cop.K.aov)
TukeyHSD(aoa.cop.K.aov, conf.level=.95)

aoa.cop.M1 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "04-28-22" ),]
aoa.cop.M1.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M1)
summary(aoa.cop.M1.aov)
TukeyHSD(aoa.cop.M1.aov, conf.level=.95)

aoa.cop.M2 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "06-01-22" ),]
aoa.cop.M2.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M2)
summary(aoa.cop.M2.aov)
TukeyHSD(aoa.cop.M2.aov, conf.level=.95)

aoa.cop.M3 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "07-05-22" ),]
aoa.cop.M3.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M3)
summary(aoa.cop.M3.aov)
TukeyHSD(aoa.cop.M3.aov, conf.level=.95)

aoa.cop.M4 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "07-20-22" ),]
aoa.cop.M4.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M4)
summary(aoa.cop.M4.aov)
TukeyHSD(aoa.cop.M4.aov, conf.level=.95)

aoa.cop.M5 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "09-13-22" ),]
aoa.cop.M5.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M5)
summary(aoa.cop.M5.aov)
TukeyHSD(aoa.cop.M5.aov, conf.level=.95)

# AOA copies between irrigations
qPCR.BS$x
qPCR.BS.ed <- qPCR.BS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

# plot
aoa.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  labs(y="AOA copy number per g dry soil", fill='Treatment')+
  scale_fill_manual(values = c("#009E73","#B6E3D7","#FF618C","#FFBBCD","#E69F00","#F4D591"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
aoa.cop.pwc.plot

aoa.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote(~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
aoa.cop.pwc.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.Bulk.tiff",
       aoa.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

# line chart AOA Copies per gram DWS Bulk Soil 

qPCR.BS.ed


# tidy up the data frame
qPCR.BS.sum <- qPCR.BS.ed %>%
  group_by(irrigation, fertilization,x, sampling.date,var3) %>%
  summarize(AOA_nbc_per_g_DW_soil=mean(AOA_nbc_per_g_DW_soil),
            sd.aoa.DWS=sd(AOA_nbc_per_g_DW_soil, na.rm = TRUE),
            AOA_nbc_per_ngDNA=mean(AOA_nbc_per_ngDNA),
            sd.aoa.ngDNA=sd(AOA_nbc_per_ngDNA, na.rm = TRUE),
            AOB_nbc_per_g_DW_soil=mean(AOB_nbc_per_g_DW_soil),
            sd.aob.DWS=sd(AOB_nbc_per_g_DW_soil, na.rm = TRUE),
            AOB_nbc_per_ngDNA=mean(AOB_nbc_per_ngDNA),
            sd.aob.ngDNA=sd(AOB_nbc_per_ngDNA, na.rm = TRUE))
qPCR.BS.sum
str(qPCR.BS.sum)

qPCR.BS.ed$lm_pred_val <- predict(aoa.BS.copies.mod,newdata = qPCR.BS.ed,
                              interval =  "confidence"
) %>% as.data.frame()


ggplot(qPCR.BS.ed, aes(x = sampling.date, y = AOA_nbc_per_g_DW_soil, linetype=irrigation)) +
  geom_jitter(position = position_jitter(0.2),
               color = "darkgray") + 
  geom_line(linewidth=1.15,aes(group = x,col=fertilization), data = qPCR.BS.sum) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"))+
  geom_errorbar(aes(ymin = AOA_nbc_per_g_DW_soil-sd.aoa.DWS, ymax = AOA_nbc_per_g_DW_soil+sd.aoa.DWS),
                data = qPCR.BS.sum, width = 0.2) +
  geom_point(data = qPCR.BS.sum, size = 2)+
  #geom_ribbon(aes(ymin = AOA_nbc_per_g_DW_soil+lm_pred_val$lwr, ymax = AOA_nbc_per_g_DW_soil-lm_pred_val$upr), alpha = 0.1)+
  ylab(bquote(~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("A. Bulk Soil")









#_______________________________________________________________________________
# AOA abundance per ng DNA

set.seed(13)
aoa.BS.copies.ngDNA.mod <- lm(qPCR.BS$AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aoa.BS.copies.ngDNA.mod)

aoa.ngDNA.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote(~italic(amoA)~'gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")
aoa.ngDNA.pwc.plot

#install.packages("patchwork")
#library(patchwork)
aoa.copies.bulk.plot <- aoa.cop.pwc.plot+aoa.ngDNA.pwc.plot + plot_layout(nrow = 2)
aoa.copies.bulk.plot

#_______________________________________________________________________________


# AOA abundance per ng DNA - RHIZOSPHERE


setwd('D:/Fina/INRAE_Project/microservices/')
qPCR <- read.csv("qPCR_results.csv")
qPCR.RS <- qPCR[121:192,]
str(qPCR.RS)
qPCR.RS$nr <- as.factor(qPCR.RS$nr)
qPCR.RS$plot <- as.factor(qPCR.RS$plot)
qPCR.RS$block <- as.factor(qPCR.RS$block)
qPCR.RS$irrigation <- as.factor(qPCR.RS$irrigation)
qPCR.RS$fertilization <- as.factor(qPCR.RS$fertilization)
qPCR.RS$type <- as.factor(qPCR.RS$type)
qPCR.RS$sampling.date <- as.factor(qPCR.RS$sampling.date)
qPCR.RS$var3 <- as.factor(qPCR.RS$var3)
qPCR.RS$x <- as.factor(qPCR.RS$x)
qPCR.RS$sampling.date <- factor(qPCR.RS$sampling.date, levels = c("28/04/2022", "01/06/2022", "05/07/2022"),
                                labels = c("04-28-22", "06-01-22", "07-05-22"))

##subset
aoa.rh.M <- qPCR.RS[which(qPCR.RS$fertilization == "M"),]
aoa.rh.M
aoa.rh.D <- qPCR.RS[which(qPCR.RS$fertilization == "D"),]
aoa.rh.D
aoa.rh.K <- qPCR.RS[which(qPCR.RS$fertilization == "K"),]
aoa.rh.K
aoa.rh.D1 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.D2 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.D3 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "07-05-22" ),]
aoa.rh.K1 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.K2 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.K3 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "07-05-22" ),]
aoa.rh.M1 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.M2 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.M3 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "07-05-22" ),]
# check stats
aoa.RS.copies <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOA_nbc_per_ngDNA, type = "mean_sd")

aoa.RS.copies.plot <- ggboxplot(
  qPCR.RS, x = "irrigation", y = "AOA_nbc_per_ngDNA",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aoa.RS.copies.plot

# check assumption (outliers)
aoa.RS.copies.out <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOA_nbc_per_ngDNA) # no extreme outliers
# Saphiro-Wilk for normality
aoa.RS.copies.SW <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOA_nbc_per_ngDNA)
ggqqplot(qPCR.RS, "AOA_nbc_per_ngDNA", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.RS.copies.Lave <- qPCR.RS %>%
  group_by(sampling.date) %>%
  levene_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.RS.copies.aov <- anova_test(
  data = qPCR.RS, dv = AOA_nbc_per_ngDNA, wid = nr ,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.RS.copies.aov)

# TWO-WAY
res.aov2.rh <- aov(AOA_nbc_per_ngDNA ~ irrigation + fertilization, data = qPCR.RS)
summary(res.aov2.rh)
# THREE WAY
res.aov3.rh <- qPCR.RS %>% anova_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date)
res.aov3.rh
# Model Fit
set.seed(13)
aoa.RS.copies.lm.mod <- lm(qPCR.RS$AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.RS)
anova(aoa.RS.copies.lm.mod) # similar to type II anova in three-way 
# Pairwise comparisons
aoa.rh.K %>%
  pairwise_t_test(AOA_nbc_per_ngDNA ~ irrigation, p.adjust.method = "BH")
# Pairwise comparisons
library(emmeans)
pwc.aoa.ngDNA.rhizo <- qPCR.RS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_ngDNA ~ irrigation, p.adjust.method = "BH") %>%
  select(-df, -statistic, -p) # Remove details
pwc.aoa.ngDNA.rhizo # significant if only i grouped by fertilization not significant per sampling date


qPCR.RS$x
qPCR.RS.ed <- qPCR.RS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")
# plot
aoa.ngDNA.rhizo.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ ggtitle("B. Rhizosphere")
aoa.ngDNA.rhizo.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.Rhizo.tiff",
       aoa.ngDNA.rhizo.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


aoa.copies.all.plot <- aoa.cop.pwc.plot+aoa.ngDNA.rhizo.plot + plot_layout(nrow = 2)
aoa.copies.all.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.ALL.tiff",
       aoa.copies.all.plot, device = "tiff",
       width = 12, height =11, 
       units= "in", dpi = 600)


############################################################################################################################
# 2. Response variable: AOB COPY NUMBER
############################################################################################################################
# 1. Bulk Soil

aob.BS.copies <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOB_nbc_per_g_DW_soil, type = "mean_sd")


aob.BS.copies.plot <- ggboxplot(
  qPCR.BS, x = "irrigation", y = "AOB_nbc_per_g_DW_soil",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aob.BS.copies.plot

# check assumption (outliers)
aob.BS.copies.out <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOB_nbc_per_g_DW_soil) # no extreme outliers
# Saphiro-Wilk for normality
aob.BS.copies.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOB_nbc_per_g_DW_soil)
ggqqplot(qPCR.BS, "AOB_nbc_per_g_DW_soil", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.BS.copies.Lave <- qPCR.BS %>%
  group_by(sampling.date) %>%
  levene_test(AOB_nbc_per_g_DW_soil ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.BS.copies.aov <- anova_test(
  data = qPCR.BS, type=3, dv = AOB_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.BS.copies.aov)

aob.bulk.DW <- aov(AOB_nbc_per_g_DW_soil~irrigation, qPCR.BS)
anova(aob.bulk.DW)




# Model Fit
set.seed(13)
aob.BS.copies.mod <- lm(qPCR.BS$AOB_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aob.BS.copies.mod)

set.seed(13)
aob.BS.copies.mod.lmer <- lmerTest::lmer(qPCR.BS$AOB_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date+(1|plot), data=qPCR.BS)
anova(aob.BS.copies.mod.lmer)

# plot
aob.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
aob.cop.pwc.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.Bulk.tiff",
       aob.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

#______________________________________________________________________________
# 2. Rhizosphere

# check assumption (outliers)
aob.RS.copies.out <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOB_nbc_per_ngDNA) # no extreme outliers
# Saphiro-Wilk for normality
aob.RS.copies.SW <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOB_nbc_per_ngDNA)
ggqqplot(qPCR.RS, "AOB_nbc_per_ngDNA", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.RS.copies.Lave <- qPCR.RS %>%
  group_by(sampling.date) %>%
  levene_test(AOB_nbc_per_ngDNA ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.RS.copies.aov <- anova_test(
  data = qPCR.RS, dv = AOB_nbc_per_ngDNA, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.RS.copies.aov)
# Model Fit


set.seed(13)
aob.RS.copies.lm.mod <- lm(qPCR.RS$AOB_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.RS)
anova(aob.RS.copies.lm.mod)


aob.RS.copies.mod <- lmerTest::lmer(qPCR.RS$AOB_nbc_per_ngDNA ~ irrigation*fertilization + (1|plot), data=qPCR.RS)
anova(aob.RS.copies.mod)


aob.test.pw.rhizo <- qPCR.RS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOB_nbc_per_ngDNA ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
aob.test.pw.rhizo

# plot
aob.ngDNA.rhizo.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOB_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ ggtitle("B. Rhizosphere")
aob.ngDNA.rhizo.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.Rhizo.tiff",
       aob.ngDNA.rhizo.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


aob.copies.all.plot <- aob.cop.pwc.plot+aob.ngDNA.rhizo.plot + plot_layout(nrow = 2)
aob.copies.all.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.ALL.tiff",
       aob.copies.all.plot, device = "tiff",
       width = 12, height =11, 
       units= "in", dpi = 600)

