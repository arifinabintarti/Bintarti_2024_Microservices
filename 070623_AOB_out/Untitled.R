##### ALPHA DIVERSITY ANALYSIS: MICROSERVICES ######

install.packages("datarium")
install.packages("rstatix")
library(datarium)
library(rstatix)
###########################################################################
# 1. Response variable: Richness
###########################################################################
# 1a. Analyses of Bulk Soil
str(aob.min.meta.bulk)
aob.min.meta.bulk.sum.rich <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aob.bulk.sum.rich.plot <- ggboxplot(
  aob.min.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.rich.plot
# check assumption (outliers)
aob.bulk.rich.out <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.rich.SW <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aob.min.meta.bulk, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.rich.Lave <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.min.meta$SampleID <- as.factor(aob.min.meta$SampleID)
aob.min.meta$PlotID <- as.factor(aob.min.meta$PlotID)
aob.min.meta.bulk$SampleID <- as.factor(aob.min.meta.bulk$SampleID)
aob.min.meta.bulk$PlotID <- as.factor(aob.min.meta.bulk$PlotID)
aob.bulk.rich.aov <- anova_test(
  data = aob.min.meta.bulk, type=3, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.rich.aov)
# RESULT: Significant effect of treatment but not irrigation within date on the AOB richness
# There are no significant three-way interaction on the AOB richnesss
# ANOVA table (type I tests)
# Effect      DFn    DFd      F        p p<.05       ges
# Irrigation  1     18     0.124      7.28e-01       0.002
# Treatment   2     18    18.130      4.87e-05     * 0.404
# Date        4     72     6.705      2.00e-03     * 0.198
# Check other methods:
# Test Method 1
model1 <- aov(aob.min.meta.bulk$Richness ~ Irrigation*Treatment*Date + Error(PlotID/(Date)) + Irrigation*Treatment, data=aob.min.meta.bulk)
summary(model1)
# Test Method 2
install.packages("lmerTest")
library(lmerTest)
str(aob.min.meta.bulk)
aob.min.meta.bulk$Date <- as.factor(aob.min.meta.bulk$Date)
aob.min.meta.bulk$var <- as.factor(aob.min.meta.bulk$var)
aob.min.meta.bulk$var2 <- as.factor(aob.min.meta.bulk$var2)
model2 <- lmerTest::lmer(aob.min.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.bulk)
anova(model2, type = 2)
# Test Method 3
model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="ML", data=aob.min.meta.bulk)
anova(model3)
# Test Method 4
aob.bulk.rich.aov <- anova_test(
  data = aob.min.meta.bulk, type=1, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.rich.aov)
# two-way interaction
aob.bulk.rich.aov.two.way <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = Richness, wid = PlotID, between = c(Irrigation, Treatment))
aob.bulk.rich.aov.two.way # there are no significant two-way interactions
# There are only main effect significance in specific date
# simple simple test for main effect
aob.rich.trt.effect <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Richness, wid = PlotID, between = Treatment)
aob.rich.trt.effect
# simple simple multiple paiwaise comparisons or posthoc test
# Fit pairwise comparisons
aob.min.meta.bulk$Treatment <- factor(aob.min.meta.bulk$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
aob.min.meta.bulk$Date  <- as.Date(aob.min.meta.bulk$Date , "%m/%d/%Y")
aob.rich.pwc <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Richness ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.rich.pwc
######################################################################
# 1b. Analyses of rhizosphere Soil
aob.min.meta.rh <- aob.min.meta[121:192,]
str(aob.min.meta.rh)
aob.min.meta.rh$SampleID <- as.factor(aob.min.meta.rh$SampleID)
aob.min.meta.rh$PlotID <- as.factor(aob.min.meta.rh$PlotID)
aob.min.meta.rh$Date <- as.factor(aob.min.meta.rh$Date)
aob.min.meta.rh.sum.rich <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aob.rh.sum.rich.plot <- ggboxplot(
  aob.min.meta.rh, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.rich.plot
# check assumption (outliers)
aob.rh.rich.out <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.rich.SW <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aob.min.meta.rh, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.rich.Lave <- aob.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.min.meta.rh$SampleID <- as.factor(aob.min.meta.rh$SampleID)
aob.min.meta.rh$PlotID <- as.factor(aob.min.meta.rh$PlotID)
aob.rh.rich.aov <- anova_test(
  data = aob.min.meta.rh, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.rich.aov)
# RESULT: Significant effect of treatment but not irrigation within date on the AOB richness
# There are no significant three-way interaction on the AOB richnesss
# ANOVA table (type I/II tests)
# Effect      DFn    DFd      F       p p<.05     ges
# Irrigation  1     18     0.133      0.72        0.002
# Treatment   2     18    13.291      0.000285  * 0.404
# Date        4     72     1.657      0.205       0.198
# Check other methods:
# Test Method 1
model1 <- aov(aob.min.meta.rh$Richness ~ Irrigation*Treatment*Date + Error(PlotID/(Date)) + Irrigation*Treatment, data=aob.min.meta.rh)
summary(model1)
# Test Method 2
#install.packages("lmerTest")
library(lmerTest)
str(aob.min.meta.rh)
model2 <- lmerTest::lmer(aob.min.meta.rh$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.rh)
anova(model2)
# Test Method 3
model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="ML", data=aob.min.meta.rh)
anova(model3)
# two-way interaction
aob.rh.rich.aov.two.way <- aob.min.meta.rh %>%
  group_by(Date) %>%
  anova_test(dv = Richness, wid = PlotID, between = c(Irrigation, Treatment))
aob.rh.rich.aov.two.way # there are no significant two-way interactions
# There are only main effect significance in specific date
# simple simple test for main effect
aob.rich.trt.effect.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Richness, wid = PlotID, between = Treatment)
aob.rich.trt.effect.rh
# simple simple multiple paiwaise comparisons or posthoc test
# Fit pairwise comparisons
aob.min.meta.rh$Treatment <- factor(aob.min.meta.rh$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
aob.min.meta.rh$Date  <- as.Date(aob.min.meta.rh$Date , "%m/%d/%Y")
aob.rich.pwc.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Richness ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.rich.pwc.rh

# pairwise comparisons for all data
aob.min.meta$Treatment <- factor(aob.min.meta$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
aob.min.meta$Date  <- as.Date(aob.min.meta$Date , "%m/%d/%Y")
aob.min.meta$Type <- factor(aob.min.meta$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere"))
aob.rich.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(Richness ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)

###########################################################################
# 2. Response variable: Shannon
###########################################################################
# 2a. Analyses of Bulk Soil
str(aob.min.meta.bulk)
aob.min.meta.bulk.sum.sha <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aob.bulk.sum.sha.plot <- ggboxplot(
  aob.min.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.sha.plot
# check assumption (outliers)
aob.bulk.sha.out <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.sha.SW <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(aob.min.meta.bulk, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.sha.Lave <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.bulk)
aob.bulk.sha.aov <- anova_test(
  data = aob.min.meta.bulk, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.sha.aov)
# other method with lmer
model2 <- lmerTest::lmer(aob.min.meta.bulk$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.bulk)
anova(model2)
# two-way interaction
aob.bulk.sha.aov.two.way <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = Shannon, wid = PlotID, between = c(Irrigation, Treatment))
aob.bulk.sha.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.sha.trt.effect.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Shannon, wid = PlotID, between = Treatment)
aob.sha.trt.effect.bulk
# simple simple multiple paiwaise comparisons or posthoc test
# Fit pairwise comparisons
str(aob.min.meta.bulk)
aob.sha.pwc.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Shannon ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.sha.pwc.bulk
##################################################################
# 2b. Analyses of Rhizosphere
str(aob.min.meta.rh)
aob.min.meta.rh.sum.sha <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aob.rh.sum.sha.plot <- ggboxplot(
  aob.min.meta.rh, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.sha.plot
# check assumption (outliers)
aob.rh.sha.out <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.sha.SW <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(aob.min.meta.rh, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.sha.Lave <- aob.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.rh)
aob.rh.sha.aov <- anova_test(
  data = aob.min.meta.rh, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.sha.aov)
# other method with lmer
model2 <- lmerTest::lmer(aob.min.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.rh)
anova(model2)
# two-way interaction
aob.rh.sha.aov.two.way <- aob.min.meta.rh %>%
  group_by(Date) %>%
  anova_test(dv = Shannon, wid = PlotID, between = c(Irrigation, Treatment))
aob.rh.sha.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.sha.trt.effect.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Shannon, wid = PlotID, between = Treatment)
aob.sha.trt.effect.rh
# simple simple multiple pairwaise comparisons or posthoc test
# Fit pairwise comparisons
aob.sha.pwc.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Shannon ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.sha.pwc.rh
###########################################################################
# 3. Response variable: Simpson
###########################################################################
# 3a. Analyses of Bulk Soil
str(aob.min.meta.bulk)
aob.min.meta.bulk.sum.simp <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Simpson, type = "mean_sd")
aob.bulk.sum.simp.plot <- ggboxplot(
  aob.min.meta.bulk, x = "Irrigation", y = "Simpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.simp.plot
# check assumption (outliers)
aob.bulk.simp.out <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Simpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.simp.SW <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Simpson)
ggqqplot(aob.min.meta.bulk, "Simpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.simp.Lave <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Simpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.bulk)
aob.bulk.simp.aov <- anova_test(
  data = aob.min.meta.bulk, dv = Simpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.simp.aov)
# other method with lmer
model2 <- lmerTest::lmer(aob.min.meta.bulk$Simpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.bulk)
anova(model2, type =2)
# two-way interaction
aob.bulk.simp.aov.two.way <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = InvSimpson, wid = PlotID, between = c(Irrigation, Treatment))
aob.bulk.simp.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.simp.trt.effect.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Simpson, wid = PlotID, between = Treatment)
aob.simp.trt.effect.bulk
# simple simple multiple paiwaise comparisons or posthoc test
# Fit pairwise comparisons
str(aob.min.meta.bulk)
aob.simp.pwc.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Simpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.simp.pwc.bulk
##################################################################
# 3b. Analyses of Rhizosphere
str(aob.min.meta.rh)
aob.min.meta.rh.sum.simp <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Simpson, type = "mean_sd")
aob.rh.sum.simp.plot <- ggboxplot(
  aob.min.meta.rh, x = "Irrigation", y = "Simpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.simp.plot
# check assumption (outliers)
aob.rh.simp.out <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Simpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.simp.SW <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Simpson)
ggqqplot(aob.min.meta.rh, "Simpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.simp.Lave <- aob.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Simpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.rh)
aob.rh.simp.aov <- anova_test(
  data = aob.min.meta.rh, type=3, dv = Simpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.simp.aov)
# other method with lmer
model2 <- lmerTest::lmer(aob.min.meta.rh$Simpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.rh)
anova(model2)
# two-way interaction
aob.rh.simp.aov.two.way <- aob.min.meta.rh %>%
  group_by(Date) %>%
  anova_test(dv = Simpson, wid = PlotID, between = c(Irrigation, Treatment))
aob.rh.simp.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.simp.trt.effect.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = Simpson, wid = PlotID, between = Treatment)
aob.simp.trt.effect.rh
# simple simple multiple pairwaise comparisons or posthoc test
# Fit pairwise comparisons
aob.simp.pwc.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(Simpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.simp.pwc.rh

###########################################################################
# 4. Response variable: Inverse Simpson
###########################################################################
# 4a. Analyses of Bulk Soil
str(aob.min.meta.bulk)
aob.min.meta.bulk.invsum.simp <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aob.bulk.sum.invsimp.plot <- ggboxplot(
  aob.min.meta.bulk, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.invsimp.plot
# check assumption (outliers)
aob.bulk.invsimp.out <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.invsimp.SW <- aob.min.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aob.min.meta.bulk, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.invsimp.Lave <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.bulk)
aob.bulk.invsimp.aov <- anova_test(
  data = aob.min.meta.bulk, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.invsimp.aov)
# other method with lmer
model2 <- lmerTest::lmer(aob.min.meta.bulk$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.bulk)
anova(model2)
# two-way interaction
aob.bulk.invsimp.aov.two.way <- aob.min.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = InvSimpson, wid = PlotID, between = c(Irrigation, Treatment))
aob.bulk.invsimp.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.invsimp.trt.effect.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = InvSimpson, wid = PlotID, between = Treatment)
aob.invsimp.trt.effect.bulk
# simple simple multiple paiwaise comparisons or posthoc test
# Fit pairwise comparisons
str(aob.min.meta.bulk)
aob.invsimp.pwc.bulk <- aob.min.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(InvSimpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.invsimp.pwc.bulk
##################################################################
# 4b. Analyses of Rhizosphere
str(aob.min.meta.rh)
aob.min.meta.rh.sum.invsimp <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aob.rh.sum.invsimp.plot <- ggboxplot(
  aob.min.meta.rh, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.invsimp.plot
# check assumption (outliers)
aob.rh.invsimp.out <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.invsimp.SW <- aob.min.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aob.min.meta.rh, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.invsimp.Lave <- aob.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.min.meta.rh)
aob.rh.invsimp.aov <- anova_test(
  data = aob.min.meta.rh, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.invsimp.aov)
# other method with lmer
aob.min.meta.rh$var2 <- as.factor(aob.min.meta.rh$var2)
model2 <- lmerTest::lmer(aob.min.meta.rh$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.min.meta.rh)
anova(model2)
# two-way interaction
aob.rh.invsimp.aov.two.way <- aob.min.meta.rh %>%
  group_by(Date) %>%
  anova_test(dv = InvSimpson, wid = PlotID, between = c(Irrigation, Treatment))
aob.rh.invsimp.aov.two.way # there are no significant two-way interactions
# There are only main effect significance of Treatment in specific date
# simple simple test for main effect
aob.invsimp.trt.effect.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  anova_test(dv = InvSimpson, wid = PlotID, between = Treatment)
aob.invsimp.trt.effect.rh
# simple simple multiple pairwaise comparisons or posthoc test
# Fit pairwise comparisons
aob.invsimp.pwc.rh <- aob.min.meta.rh %>%
  group_by(Date, Irrigation) %>%
  pairwise_t_test(InvSimpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif) # Remove details
aob.invsimp.pwc.rh









