#############################################################################################
# Analysis of amoA gene of AOA Illumina MiSeq Data 
#############################################################################################

# Date : 08 February 2024
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
# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA.ASV-analysis')
#setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA.ASV-analysis')
wd <- print(getwd())
# load the asv table
aoa.asv <- read.table('annotated.AOA.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(aoa.asv) # 646  192
sort(rowSums(aoa.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/')
#setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/')
aoa.tax <- read.table("besthit.diamond.output.curateddb.AOA.ASVs.edited.csv", sep = ';', header=T)
dim(aoa.tax) # 646
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA-rooted-tree/')
AOA_rooted_tree <- ape::read.tree("tree.nwk")
AOA_rooted_tree

## make a phyloseq object of the asv table, taxonomy table, metadata
# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aoa.asv))
aoa.asv.ord  <- aoa.asv[ ,re_order]
aoa.asv.physeq <- otu_table(aoa.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aoa.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aoa.asv.physeq) <- paste0("S", sample_names(aoa.asv.physeq))
# phyloseq object of the taxonomy table
aoa.tax <- column_to_rownames(aoa.tax, var = "ASVid")
#row.names(aoa.tax) <- aoa.tax$ASVid
aoa.tax.physeq <- tax_table(as.matrix(aoa.tax)) # taxonomy table
# phyloseq object of the metadata
meta_micro$Date <- factor(meta_micro$Date, levels = c("4/28/22", "6/1/22", "7/5/22", "7/20/22", "9/13/22"),
                          labels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"))
rownames(meta_micro) <- sample_names(aoa.asv.physeq)
aoa.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aoa.meta.physeq)
# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA-rooted-tree/')
AOA_rooted_tree <- ape::read.tree("tree.nwk")
# make phyloseq object
aoa.physeq <- merge_phyloseq(aoa.asv.physeq,aoa.tax.physeq,aoa.meta.physeq,AOA_rooted_tree)
aoa.physeq
sample_data(aoa.physeq)$SampleID <- paste0("S", sample_data(aoa.physeq)$SampleID)
sample_data(aoa.physeq)
#remove the last sampling date
aoa.physeq.no.last = subset_samples(aoa.physeq, Date != "Sept 13th")
View(sample_data(aoa.physeq.no.last))
aoa.physeq.no.last # 646 samples BS & RS
otu_table(aoa.physeq.no.last)
sort(taxa_sums(aoa.physeq.no.last), decreasing =F) # contain ASVs that don't exist in at least 1 sample
aoa.physeq.no.last1 <- prune_taxa(taxa_sums(aoa.physeq.no.last)>0, aoa.physeq.no.last)
aoa.physeq.no.last1 # 607 taxa & 168 samples BS & RS
# rarefy to minimum sequencing depth (minimum reads = 3832 reads)
sort(colSums(otu_table(aoa.physeq.no.last1), na.rm = FALSE, dims = 1), decreasing = F)
set.seed(13)
aoa.rare.nolast <- rarefy_even_depth(aoa.physeq.no.last1, sample.size = min(sample_sums(aoa.physeq.no.last1)),
  rngseed = 13, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
aoa.rare.nolast
sort(sample_sums(aoa.rare.nolast), decreasing = F) # 51 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed
sort(rowSums(otu_table(aoa.rare.nolast), na.rm = FALSE, dims = 1), decreasing = F)
# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 
#### AOA alpha diversity using rarefied data ###
colSums(otu_table(aoa.rare.nolast))
aoa.asv.rare.df <- as.data.frame(otu_table(aoa.rare.nolast))
dim(aoa.asv.rare.df) #  556 ASVs 168 Samples
aoa.asv.rare.df_pa <- 1*(aoa.asv.rare.df>0)
aoa.s2 <- specnumber(aoa.asv.rare.df, MARGIN = 2) # richness
aoa.richness2 <- as.data.frame(aoa.s2) 
aoa.h2 <- diversity(t(aoa.asv.rare.df), index = 'shannon') # Shannon index
aoa.shannon2 <- as.data.frame(aoa.h2)
aoa.d2 <- diversity(t(aoa.asv.rare.df), index = 'simpson') # Simpson index
aoa.simpson2 <- as.data.frame(aoa.d2)
aoa.ainv.d2 <- diversity(t(aoa.asv.rare.df), index = 'invsimpson')
# edit env mapping data of AOA
aoa.meta.df2 <- data.frame(meta_micro)
head(aoa.meta.df2)
dim(aoa.meta.df2) #192 Samples
# filter out S11 from the metadata and the last sampling date
aoa.meta.df.sub2 <- aoa.meta.df2 %>% filter(Date != "Sept 13th")
dim(aoa.meta.df.sub2) # 168 Samples
aoa.meta.df.sub2$Date <- factor(aoa.meta.df.sub2$Date)
# adding alpha diversity calculation
aoa.meta.df.sub2$AOA_Richness <- aoa.s2
aoa.meta.df.sub2$AOA_Shannon <- aoa.h2
aoa.meta.df.sub2$AOA_Simpson <- aoa.d2
aoa.meta.df.sub2$AOA_InvSimpson <- aoa.ainv.d2
# tidy up
str(aoa.meta.df.sub2)
aoa.meta.df.sub2$Type <- factor(aoa.meta.df.sub2$Type, levels = c("BS", "RS"),
                               labels = c("Bulk_Soil", "Rhizosphere"))
aoa.meta.df.sub2$Treatment <- factor(aoa.meta.df.sub2$Treatment, levels = c("D", "K", "M"),
                                    labels = c("BIODYN", "CONFYM", "CONMIN"))
aoa.meta.df.sub2$SampleID<-as.factor(aoa.meta.df.sub2$SampleID)
aoa.meta.df.sub2$PlotID<-as.factor(aoa.meta.df.sub2$PlotID)
aoa.meta.df.sub2$Irrigation<-as.factor(aoa.meta.df.sub2$Irrigation)
aoa.meta.df.sub2$Block<-as.factor(aoa.meta.df.sub2$Block)
aoa.meta.df.sub2$x<-as.factor(aoa.meta.df.sub2$x)
aoa.meta.df.sub2$rep<-as.factor(aoa.meta.df.sub2$rep)
aoa.meta.df.sub2$rep2<-as.factor(aoa.meta.df.sub2$rep2)
aoa.meta.df.sub2$var2<-as.factor(aoa.meta.df.sub2$var2)
aoa.meta.df.sub2$var3<-as.factor(aoa.meta.df.sub2$var3)
aoa.meta.df.sub2[sapply(aoa.meta.df.sub2, is.character)] <- 
 lapply(aoa.meta.df.sub2[sapply(aoa.meta.df.sub2, is.character)], as.numeric)
aoa.meta.df.sub2[sapply(aoa.meta.df.sub2, is.integer)] <- 
 lapply(aoa.meta.df.sub2[sapply(aoa.meta.df.sub2, is.integer)], as.numeric)

###########################################################################
# 1. Response variable: Alpha Diversity
###########################################################################
install.packages("datarium")
install.packages("rstatix")
library(datarium)
library(rstatix)

# 1a. Analyses of Bulk Soil
aoa.meta.bulk2 <- aoa.meta.df.sub2[1:96,]
str(aoa.meta.bulk2)
dim(aoa.meta.bulk2)
aoa.meta.bulk2$SampleID<-factor(aoa.meta.bulk2$SampleID)
aoa.meta.bulk2$PlotID<-factor(aoa.meta.bulk2$PlotID)
aoa.meta.bulk2$Irrigation<-factor(aoa.meta.bulk2$Irrigation)
aoa.meta.bulk2$Block<-factor(aoa.meta.bulk2$Block)
aoa.meta.bulk2$x<-factor(aoa.meta.bulk2$x)
aoa.meta.bulk2$rep<-factor(aoa.meta.bulk2$rep)
aoa.meta.bulk2$rep2<-factor(aoa.meta.bulk2$rep2)
aoa.meta.bulk2$var2<-factor(aoa.meta.bulk2$var2)
aoa.meta.bulk2$var3<-factor(aoa.meta.bulk2$var3)
aoa.meta.bulk2$AOA_sqrtRich <- sqrt(aoa.meta.bulk2$AOA_Richness)


aoa.BS.sum.rich <- aoa.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(AOA_Richness, type = "mean_sd")
View(aoa.BS.sum.rich)
aoa.BS.rich.plot <- ggboxplot(
  aoa.meta.bulk2, x = "Irrigation", y = "AOA_Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.BS.rich.plot
# check assumption (outliers)
aoa.BS.rich.out <- aoa.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(AOA_Richness) # no extreme outliers
View(aoa.BS.rich.out)
# Saphiro-Wilk for normality
aoa.BS.rich.SW <- aoa.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(AOA_Richness)
View(aoa.BS.rich.SW)
ggqqplot(aoa.meta.bulk2, "AOA_Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.BS.rich.Lave <- aoa.meta.bulk2 %>%
  group_by(Date) %>%
  levene_test(AOA_Richness ~ Irrigation*Treatment)
View(aoa.BS.rich.Lave)
# interaction plot
with(aoa.meta.bulk2, interaction.plot(x.factor = Treatment, 
                               trace.factor = Irrigation, 
                               response = AOA_Richness))
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.meta.bulk2$SampleID <- factor(aoa.meta.bulk2$SampleID)
aoa.BS.rich.aov <- anova_test(
  data = aoa.meta.bulk2, type=3, dv = AOA_Richness, wid = SampleID,
  between = c(Irrigation, Treatment, Date))
get_anova_table(aoa.BS.rich.aov)

# Three-Way Repeated-Measures ANOVA
aoa.bulk.rich.aov2 <- anova_test(
  data = aoa.meta.bulk2, type=2,dv = AOA_Richness, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.rich.aov2)


# with afex package
library(afex)
install.packages("ez")
aoa.BS.rich.afex <- aov_ez("PlotID", "AOA_Richness", aoa.meta.bulk2, 
             within = "Date",
             between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
aoa.BS.rich.afex
# Test Method 3 - Richness
aoa.BS.rich.mod1 <- lme(AOA_Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=aoa.meta.bulk2, na.action = na.omit)
anova(aoa.BS.rich.mod1) 
aoa.BS.rich.mod2 <- lmerTest::lmer(aoa.meta.bulk2$AOA_Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.bulk2)
anova(aoa.BS.rich.mod2)
VarCorr(aoa.BS.rich.mod2)
aoa.BS.rich.mod3 <- lme4::lmer(aoa.meta.bulk2$AOA_Richness ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                   contrasts = list(irrigation="contr.sum",fertilization="contr.sum",sampling.date="contr.sum"),data=aoa.meta.bulk2)#control = lmerControl(check.nobs.vs.nRE="ignore"))
car::Anova(aoa.BS.rich.mod3, test="F", type="III") 
summary(aoa.BS.rich.mod3)
VarCorr(aoa.BS.rich.mod3)

anova(aoa.BS.rich.mod3,aoa.BS.rich.mod2)
shapiro.test(resid(aoa.BS.rich.mod3)) # not normal
plot(simulateResiduals(aoa.BS.rich.mod3)) # okay

#confidence interval
confint(aoa.BS.rich.mod3, oldNames = FALSE)

# Test Method 3 - Shannon Index

# Three-Way Repeated-Measures ANOVA
aoa.bulk.sha.aov2 <- anova_test(
  data = aoa.meta.bulk2, type=2,dv = AOA_Shannon, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.sha.aov2)


aoa.BS.sha.mod1 <- lme(AOA_Shannon ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=aoa.meta.bulk2, na.action = na.omit)
anova(aoa.BS.sha.mod1) 
aoa.BS.sha.mod2 <- lmerTest::lmer(aoa.meta.bulk2$AOA_Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.bulk2, na.action = na.omit)
anova(aoa.BS.sha.mod2)
aoa.BS.sha.mod3 <- lmerTest::lmer(aoa.meta.bulk2$AOA_Shannon ~ Irrigation*Treatment*Date+(1|Block:Date), data=aoa.meta.bulk2, na.action=na.omit)
anova(aoa.BS.sha.mod3)

library(DHARMa)
shapiro.test(resid(aoa.BS.sha.mod3)) # not normal
plot(simulateResiduals(aoa.BS.sha.mod3)) # okay

# with afex package
library(afex)
aoa.BS.sha.afex <- aov_ez("rep2", "AOA_Shannon", aoa.meta.bulk2, 
             within = c("Irrigation","Treatment", "Date"),
             #between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
aoa.BS.sha.afex

#BETA DIVERSITY ANALYSIS #

# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# Bray-Curtis - Bulk Soil :
aoa.asv.rare.df <- as.data.frame(otu_table(aoa.rare.nolast))
aoa.asv4.BS <- aoa.asv.rare.df[,1:96]
colnames(aoa.asv4.BS)
aoa.asv4.BS1 <- aoa.asv4.BS[rowSums(aoa.asv4.BS)>0,]
sort(rowSums(aoa.asv4.BS1, na.rm = FALSE, dims = 1), decreasing = FALSE)
setwd('D:/Fina/INRAE_Project/microservices/070623_AOA_out/')
#write.csv(aoa.asv.bulk1, file = "aoa.bulk1.otutab.csv")
aoa4.BS_dist_bc <- vegdist(t(aoa.asv4.BS1), method = "bray")

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
# Bray-Curtis - Bulk Soil:
aoa4.BS_pcoa_bc <- cmdscale(aoa4.BS_dist_bc, eig=T)

# 3. scores of PC1 and PC2
# Bray-Curtis - Bulk Soil:
ax1.scores4.aoa.BS <- aoa4.BS_pcoa_bc$points[,1]
ax2.scores4.aoa.BS <- aoa4.BS_pcoa_bc$points[,2] 

# 4. calculate percent variance explained, then add to plot
aoa.meta.bulk2 <- aoa.meta.df.sub2[1:96,]
# Bray-curtis - Bulk Soil:
ax1.aoa.BS4 <- aoa4.BS_pcoa_bc$eig[1]/sum(aoa4.BS_pcoa_bc$eig)
ax2.aoa.BS4 <- aoa4.BS_pcoa_bc$eig[2]/sum(aoa4.BS_pcoa_bc$eig)
aoa.map.pcoa.BS4 <- cbind(aoa.meta.bulk2,ax1.scores4.aoa.BS,ax2.scores4.aoa.BS)
#################################################################################################

# RHIZOSPHERE

# Bray-Curtis - Rhizosphere :
aoa.asv.rare.df <- as.data.frame(otu_table(aoa.rare.nolast))
colnames(aoa.asv.rare.df)
aoa.asv4.RS <- aoa.asv.rare.df[,97:168]
colnames(aoa.asv4.RS)
aoa.asv4.RS1 <- aoa.asv4.RS[rowSums(aoa.asv4.RS)>0,]
sort(rowSums(aoa.asv4.RS1, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(aoa.asv4.RS1) #847 72
aoa4.RS_dist_bc <- vegdist(t(aoa.asv4.RS1), method = "bray")

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
# Bray-Curtis - Rhizosphere :
aoa4.RS_pcoa_bc <- cmdscale(aoa4.RS_dist_bc, eig=T)

# 3. scores of PC1 and PC2
# Bray-Curtis - Rhizosphere :
ax1.scores4.aoa.RS <- aoa4.RS_pcoa_bc$points[,1]
ax2.scores4.aoa.RS <- aoa4.RS_pcoa_bc$points[,2] 

# 4. calculate percent variance explained, then add to plot
aoa.meta.rh2 <- aoa.meta.df.sub2[97:168,]
view(aoa.meta.rh2)
# Bray-curtis - Rhizosphere :
ax1.aoa.RS4 <- aoa4.RS_pcoa_bc$eig[1]/sum(aoa4.RS_pcoa_bc$eig)
ax2.aoa.RS4 <- aoa4.RS_pcoa_bc$eig[2]/sum(aoa4.RS_pcoa_bc$eig)
aoa.map.pcoa.RS4 <- cbind(aoa.meta.rh2,ax1.scores4.aoa.RS,ax2.scores4.aoa.RS)

###############################################################################
# Adonis
# 1. Using adonis2 package with defined perm to restrict the permutation 
str(aoa.meta.bulk2)
set.seed(13)
aoa.adonis.BS4.bc <- adonis2(aoa4.BS_dist_bc ~ Irrigation, strata=aoa.meta.bulk2$Block, data=aoa.meta.bulk2, 
                              permutations = 999) # not significant
aoa.adonis.BS4.bc

# similar with below:
perm1.aoa = how(within = Within(type="free"), 
            plots = Plots(type = "none"),
            blocks = aoa.meta.bulk2$Block,
            nperm = 999,
            observed = TRUE)
set.seed(13)
aoa.adonis.BS.bc.perm1 <- adonis2(aoa4.BS_dist_bc ~ Irrigation, data=aoa.meta.bulk2, 
                                    permutations = perm1.aoa)
aoa.adonis.BS.bc.perm1 

# another way to use how()
# Since our intent is to focus on the variation among treatments, 
# we need to restrict the permutations so that plots are permuted within each block, but plots are not permuted across blocks.
# these two ways are equivalent:
CTRL.t1.aoa <- how(within = Within(type = "free"),
               plots = Plots(type = "none"),
               blocks = aoa.meta.bulk2$Block,
               nperm = 999,
               observed = TRUE)
# and
CTRL.t2.aoa <- how(within = Within(type = "free"),
               plots = Plots(strata = aoa.meta.bulk2$Block, type = "none"),
               nperm = 999,
               observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(13)
aoa.adonis.BS.bc.CTRL.t2 <- adonis2(aoa4.BS_dist_bc ~ Irrigation, data=aoa.meta.bulk2, 
                                      permutations = CTRL.t2.aoa)
aoa.adonis.BS.bc.CTRL.t2
# 2. Using ANOSIM package and define the strata
set.seed(13)
aoa.BS.bc.anosim <- anosim(aoa4.BS_dist_bc,
                        grouping = aoa.meta.bulk2$Irrigation, permutations = 999, strata = aoa.meta.bulk2$Block)
summary(aoa.BS.bc.anosim) # almost 0.057
# test the permanova for whole plot
set.seed(13)
aoa.adonis.BS4 <- adonis2(aoa4.BS_dist_bc ~ Irrigation*Treatment*Date, data=aoa.meta.bulk2, 
                           permutation=999) # only treatment is significant
aoa.adonis.BS4



