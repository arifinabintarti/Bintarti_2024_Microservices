###########################################################################
#################### ALPHA DIVERSITY ANALYSIS: COM ########################
###########################################################################

# Author: Ari Fina Bintarti
# Date: 02/08/2023

#install.packages("datarium")
#install.packages("rstatix")
library(datarium)
library(rstatix)
###########################################################################
# 1. Response variable: Richness
###########################################################################
# 1a. Analyses of Bulk Soil
com.meta.bulk <- com.meta.df[1:118,]
str(com.meta.bulk)
com.meta.bulk.sum.rich <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
com.bulk.sum.rich.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.rich.plot
# check assumption (outliers)
com.bulk.rich.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
com.bulk.rich.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(com.meta.bulk, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.rich.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
com.bulk.rich.aov <- anova_test(
  data = com.meta.bulk, type=2, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.bulk.rich.aov)

# Test Method 3 
#model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="REML", data=com.meta.bulk)
#anova(model3)
#model1 <- aov(com.meta.bulk$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment) , data=com.meta.bulk)
#summary(model1)

############################################################################################################
# Model Fit
set.seed(13)
com.rich.bulk.mod <- lmerTest::lmer(com.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk)
anova(com.rich.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.rich.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.bulk.mod)
# 2. between irrigation:
com.emm.rich.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.bulk.mod)

######################################################################
# 1b. Analyses of rhizosphere Soil

com.meta.rh <- com.meta.df[119:190,]
str(com.meta.rh)
com.meta.rh.sum.rich <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
com.rh.sum.rich.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.rich.plot
# check assumption (outliers)
com.rh.rich.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.rich.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(com.meta.rh, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.rich.Lave <- com.meta.rh %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.rh.rich.aov <- anova_test(
  data = com.meta.rh, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.rich.aov)
# other model
# model1 <- aov(com.meta.rh$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment), data=com.meta.rh)
# summary(model1)
# Test Method 3
# model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="ML", data=com.meta.rh)
# anova(model3)

# Model Fit
set.seed(13)
com.rich.rhizo.mod <- lmerTest::lmer(com.meta.rh$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.rh)
anova(com.rich.rhizo.mod, type = 2)
# pairwise comparisons
# 1. between fertilization treatment:
com.emm.rich.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.rhizo.mod)
# 2. between irrigation:
com.emm.rich.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.rhizo.mod)


###########################################################################
# 2. Response variable: Shannon
###########################################################################
# 2a. Analyses of Bulk Soil
com.meta.bulk.sum.sha <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
com.bulk.sum.sha.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.sha.plot
# check assumption (outliers)
com.bulk.sha.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
com.bulk.sha.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(com.meta.bulk, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.sha.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.bulk.sha.aov <- anova_test(
  data = com.meta.bulk, type=2, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.bulk.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.sha.bulk.mod <- lmerTest::lmer(com.meta.bulk$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk)
anova(com.sha.bulk.mod, type = 3)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.sha.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.bulk.mod)
# 2. between irrigation:
com.emm.sha.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.bulk.mod)
#############################################################################################################

##################################################################
# 2b. Analyses of Rhizosphere
com.meta.rh.sum.sha <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
com.rh.sum.sha.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.sha.plot
# check assumption (outliers)
com.rh.sha.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.sha.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(com.meta.rh, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.sha.Lave <- com.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(com.meta.rh)
com.rh.sha.aov <- anova_test(
  data = com.meta.rh, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.sha.rh.mod <- lmerTest::lmer(com.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.rh)
anova(com.sha.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.sha.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.rh.mod)
# 2. between irrigation:
com.emm.sha.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.rh.mod)
#############################################################################################################

###########################################################################
# 3. Response variable: Inverse Simpson
###########################################################################
# 4a. Analyses of Bulk Soil
com.meta.bulk.invsum.simp <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
com.bulk.sum.invsimp.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.invsimp.plot
# check assumption (outliers)
com.bulk.invsimp.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
com.bulk.invsimp.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(com.meta.bulk, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.invsimp.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.bulk.invsimp.aov <- anova_test(
  data = com.meta.bulk, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.bulk.invsimp.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.invsimp.bulk.mod <- lmerTest::lmer(com.meta.bulk$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk)
anova(com.invsimp.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.invsimp.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.bulk.mod)
# 2. between irrigation:
com.emm.invsimp.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.bulk.mod)
#############################################################################################################

##################################################################
# 3b. Analyses of Rhizosphere
com.meta.rh.sum.invsimp <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
com.rh.sum.invsimp.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.invsimp.plot
# check assumption (outliers)
com.rh.invsimp.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.invsimp.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(com.meta.rh, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.invsimp.Lave <- com.meta.rh %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.rh.invsimp.aov <- anova_test(
  data = com.min.meta.rh, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.invsimp.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.invsimp.rh.mod <- lmerTest::lmer(com.meta.rh$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.rh)
anova(com.invsimp.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.invsimp.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.rh.mod)
# 2. between irrigation:
com.emm.invsimp.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.rh.mod)
#############################################################################################################













