################################################################################
# Microservices Project: Soil Properties Statistical Analysis
################################################################################
# Date : 23 August 2023
# Author : Ari Fina Bintarti

library(datarium)
library(rstatix)

# set working directory
setwd('D:/Fina/INRAE_Project/microservices/')
wd <- print(getwd())
# load the metadata
soilprop <- read.csv('Soil_SamplingPlotData_DOK2022_edited.csv')
soilprop$Date <- factor(soilprop$Date, levels = c("28/04/2022", "01/06/2022", "05/07/2022", "20/07/2022", "13/09/2022"),
                          labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))
soilprop$Treatment <- factor(soilprop$Treatment, levels = c("D", "K", "M"),
                                    labels = c("Biodynamic (D)", "Conventional (K)", "Mineral fertilized (M)"))
###########################################################################
# 1. Response variable: GWC (Gravimetric water content in water content (g) / dry soil (g))
str(soilprop)
gwc.sum <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  get_summary_stats(GWC, type = "mean_sd")
gwc.sum.plot <- ggboxplot(
  soilprop, x = "Type", y = "GWC",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
gwc.sum.plot
# check assumption (outliers)
gwc.out <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  identify_outliers(GWC) # no extreme outliers
# Saphiro-Wilk for normality
gwc.SW <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  shapiro_test(GWC)
ggqqplot(soilprop, "GWC", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
gwc.Lave <- soilprop %>%
  group_by(Date) %>%
  levene_test(GWC ~ Type*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
gwc.3aov <- anova_test(
  data = soilprop, type=3, dv = GWC, wid = PlotID,
  within = Date, between = c(Type, Treatment))
get_anova_table(gwc.3aov)
# Test Method 3 
gwc.lme <- lme(GWC ~ Type*Treatment*Date, random=~1 | PlotID, method="ML", data=soilprop)
anova(gwc.lme, type = "III") # similar results as the three way mixed ANOVA above!!!
gwc.aov <- aov(soilprop$GWC ~ Type*Treatment*Date + Error(PlotID/Irrigation*Treatment) , data=aob.meta.bulk)
#summary(model1)

############################################################################################################
# Model Fit
set.seed(13)
rich.bulk.mod <- lmerTest::lmer(aob.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk)
anova(rich.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.rich.bulk <- aob.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.bulk.mod)
# 2. between irrigation:
emm.rich.irri.bulk <- aob.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.bulk.mod)
#############################################################################################################






























