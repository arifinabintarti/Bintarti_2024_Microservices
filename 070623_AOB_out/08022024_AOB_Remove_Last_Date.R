#############################################################################################
# Analysis of amoA gene of AOB Illumina MiSeq Data 
#############################################################################################

# Date : 07 February 2024
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
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.ASV-analysis')
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/AOB.ASV-analysis')
wd <- print(getwd())
# load the asv table
aob.asv <- read.table('annotated.AOB.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out')
#write.csv(aob.asv, file = "aob.asv.csv")
aob.asv
dim(aob.asv)# 1338  192
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out')
aob.tax <- read.csv("besthit.diamond.output.curateddb.AOB.ASVs.csv")
dim(aob.tax) # 1338
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices')
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.Phylogenetic-analysis/')
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/AOB.Phylogenetic-analysis')
aob.tre <- ape::read.tree("tree.AOB.nwk")

## make a phyloseq object of the asv table, taxonomy table, metadata
# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aob.asv))
aob.asv.ord  <- aob.asv[ ,re_order]
aob.asv.physeq = otu_table(aob.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aob.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aob.asv.physeq) <- paste0("S", sample_names(aob.asv.physeq))
# phyloseq object of the taxonomy table
aob.tax <- column_to_rownames(aob.tax, var = "ASVid")
aob.tax.physeq = tax_table(as.matrix(aob.tax)) # taxonomy table
# phyloseq object of the metadata
str(meta_micro)
meta_micro$Date <- factor(meta_micro$Date, levels = c("4/28/22", "06/01/2022", "07/05/2022", "7/20/22", "9/13/22"),
                          labels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"))
rownames(meta_micro) <- sample_names(aob.asv.physeq)
aob.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aob.meta.physeq)
# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/AOB-rooted-tree/')
AOB_rooted_tree <- ape::read.tree("tree.nwk")
AOB_rooted_tree
# make phyloseq object
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq,AOB_rooted_tree)
aob.physeq #192 samples BS & RS
aob.asv.ord <- as.data.frame(otu_table(aob.physeq))
aob.asv.ord
#write.csv(aob.asv.ord, file = "aob.asv.ord.csv")
sample_data(aob.physeq)$SampleID <- paste0("S", sample_data(aob.physeq)$SampleID)
sample_data(aob.physeq)
#remove the last sampling date
aob.physeq.no.last = subset_samples(aob.physeq, Date != "Sept 13th")
View(sample_data(aob.physeq.no.last))
aob.physeq.no.last #168 samples BS & RS
otu_table(aob.physeq.no.last)
sort(taxa_sums(aob.physeq.no.last), decreasing =F) # contain ASVs that don't exist in at least 1 sample
aob.physeq.no.last1 <- prune_taxa(taxa_sums(aob.physeq.no.last)>0, aob.physeq.no.last)
aob.physeq.no.last1 # 1262 taxa & 168 samples BS & RS

# rarefaction to 1282 reads (remove sample with reads < 1000)
# ASV Table
sort(colSums(otu_table(aob.physeq.no.last1), na.rm = FALSE, dims = 1), decreasing = F)
set.seed(333)
aob.rare.nolast <- rarefy_even_depth(aob.physeq.no.last1, sample.size = 1282,
                                       rngseed = 333, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.nolast # 1 samples removed (S11), 100 ASVs were removed
sort(rowSums(otu_table(aob.rare.nolast), na.rm = FALSE, dims = 1), decreasing = F)

# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 
#### AOB alpha diversity using rarefied data ###
colSums(otu_table(aob.rare.nolast))
aob.asv.rare.df <- as.data.frame(otu_table(aob.rare.nolast))
dim(aob.asv.rare.df) # 1162 ASVs
aob.asv.rare.df_pa <- 1*(aob.asv.rare.df>0)
aob.s2 <- specnumber(aob.asv.rare.df, MARGIN = 2) # richness
aob.richness2 <- as.data.frame(aob.s2) 
aob.h2 <- diversity(t(aob.asv.rare.df), index = 'shannon') # Shannon index
aob.shannon2 <- as.data.frame(aob.h2)
aob.d2 <- diversity(t(aob.asv.rare.df), index = 'simpson') # Simpson index
aob.simpson2 <- as.data.frame(aob.d2)
aob.inv.d2 <- diversity(t(aob.asv.rare.df), index = 'invsimpson')

# edit env mapping data of AOB
aob.meta.df2 <- data.frame(meta_micro)
head(aob.meta.df2)
dim(aob.meta.df2)
# filter out S11 from the metadata and the last sampling date
aob.meta.df.sub2 <- aob.meta.df2 %>% filter(SampleID != 11,
                                          Date != "Sept 13th")
str(aob.meta.df.sub2)
aob.meta.df.sub2$Date <- factor(aob.meta.df.sub2$Date)
# adding alpha diversity calculation
aob.meta.df.sub2$AOB_Richness <- aob.s2
aob.meta.df.sub2$AOB_Shannon <- aob.h2
aob.meta.df.sub2$AOB_Simpson <- aob.d2
aob.meta.df.sub2$AOB_InvSimpson <- aob.inv.d2
# tidy up
str(aob.meta.df.sub2)
aob.meta.df.sub2$Type <- factor(aob.meta.df.sub2$Type, levels = c("BS", "RS"),
                               labels = c("Bulk_Soil", "Rhizosphere"))
aob.meta.df.sub2$Treatment <- factor(aob.meta.df.sub2$Treatment, levels = c("D", "K", "M"),
                                    labels = c("BIODYN", "CONFYM", "CONMIN"))
aob.meta.df.sub2$SampleID<-as.factor(aob.meta.df.sub2$SampleID)
aob.meta.df.sub2$PlotID<-as.factor(aob.meta.df.sub2$PlotID)
aob.meta.df.sub2$Irrigation<-as.factor(aob.meta.df.sub2$Irrigation)
aob.meta.df.sub2$Block<-as.factor(aob.meta.df.sub2$Block)
aob.meta.df.sub2$x<-as.factor(aob.meta.df.sub2$x)
aob.meta.df.sub2$var<-as.factor(aob.meta.df.sub2$var)
aob.meta.df.sub2$var2<-as.factor(aob.meta.df.sub2$var2)
aob.meta.df.sub2$var3<-as.factor(aob.meta.df.sub2$var3)
aob.meta.df.sub2[sapply(aob.meta.df.sub2, is.character)] <- 
 lapply(aob.meta.df.sub2[sapply(aob.meta.df.sub2, is.character)], as.numeric)
aob.meta.df.sub2[sapply(aob.meta.df.sub2, is.integer)] <- 
 lapply(aob.meta.df.sub2[sapply(aob.meta.df.sub2, is.integer)], as.numeric)

###########################################################################
# 1. Response variable: Alpha Diversity
###########################################################################
install.packages("datarium")
install.packages("rstatix")
library(datarium)
library(rstatix)

# 1a. Analyses of Bulk Soil
aob.meta.bulk2 <- aob.meta.df.sub2[1:95,]
str(aob.meta.bulk2)
aob.BS.sum.rich <- aob.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(AOB_Richness, type = "mean_sd")
aob.BS.sum.rich
aob.BS.rich.plot <- ggboxplot(
  aob.meta.bulk2, x = "Irrigation", y = "AOB_Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.BS.rich.plot
# check assumption (outliers)
aob.BS.rich.out <- aob.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(AOB_Richness) # no extreme outliers
View(aob.BS.rich.out)
# Saphiro-Wilk for normality
aob.BS.rich.SW <- aob.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(AOB_Richness)
View(aob.BS.rich.SW)
ggqqplot(aob.meta.bulk2, "AOB_Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.BS.rich.Lave <- aob.meta.bulk2 %>%
  group_by(Date) %>%
  levene_test(AOB_Richness ~ Irrigation*Treatment)
View(aob.BS.rich.Lave)

# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.BS.rich.aov <- anova_test(
  data = aob.meta.bulk2, type=3, dv = AOB_Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.BS.rich.aov)
# with afex package
library(afex)
aob.BS.rich.afex <- aov_ez("PlotID", "AOB_Richness", aob.meta.bulk2, 
             within = "Date",
             between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
aob.BS.rich.afex
# Test Method 3 - Richness
aob.BS.rich.mod1 <- lme(AOB_Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=aob.meta.bulk2, na.action = na.omit)
anova(aob.BS.rich.mod1) 
aob.BS.rich.mod2 <- lmerTest::lmer(aob.meta.bulk2$AOB_Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk2, na.action = na.omit)
anova(aob.BS.rich.mod2)
aob.BS.rich.mod3 <- lmerTest::lmer(aob.meta.bulk2$AOB_Richness ~ Irrigation*Treatment*Date+(1|Block:Date), data=aob.meta.bulk2, na.action=na.omit)
anova(aob.BS.rich.mod3)

# Test Method 3 - Shannon Index
aob.BS.sha.mod1 <- lme(AOB_Shannon ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=aob.meta.bulk2, na.action = na.omit)
anova(aob.BS.sha.mod1) 
aob.BS.sha.mod2 <- lmerTest::lmer(aob.meta.bulk2$AOB_Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk2, na.action = na.omit)
anova(aob.BS.sha.mod2)
aob.BS.sha.mod3 <- lmerTest::lmer(aob.meta.bulk2$AOB_Shannon ~ Irrigation*Treatment*Date+(1|Block:Date), data=aob.meta.bulk2, na.action=na.omit)
anova(aob.BS.sha.mod3)
# with afex package
library(afex)
aob.BS.sha.afex <- aov_ez("PlotID", "AOB_Shannon", aob.meta.bulk2, 
             within = "Date",
             between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
aob.BS.sha.afex

