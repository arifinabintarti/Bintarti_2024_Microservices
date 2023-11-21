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
install.packages("car")
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
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = block,
  within = fertilization, between = irrigation)
get_anova_table(aoa.BS.copies.aov)
############################################################################################################
# Model Fit
set.seed(13)
#install.packages("lmerTest", type = "source")
library(lmerTest)
aoa.BS.copies.mod <- lme4::lmer(qPCR.BS$AOA_nbc_per_g_DW_soil ~ irrigation*fertilization +(fertilization:irrigation)+(1|block), data=qPCR.BS)
car::Anova(aoa.BS.copies.mod)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.rich.bulk <- aoa.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.bulk.mod)
# 2. between irrigation:
aoa.emm.rich.irri.bulk <- aoa.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.bulk.mod)

######################################################################












































