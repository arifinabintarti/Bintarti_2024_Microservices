#############################################################################################
# Analysis of amoA gene of COMAMMOX Illumina MiSeq Data 
#############################################################################################

# Date : 09 February 2024
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
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM.ASV-analysis')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM.ASV-analysis')
wd <- print(getwd())
# load the asv table
com.asv <- read.table('annotated.COM.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(com.asv) # 686 192
# remove the bad sample (sample # 26) from the OTU table
com.asv.sub <- com.asv[, -which(names(com.asv) == "26" )]
sort(rowSums(com.asv.sub, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(com.asv.sub)
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/')
com.tax <- read.table("besthit.diamond.output.curateddb.COM.ASVs.edited.csv", sep = ';', header=T)
dim(com.tax) # 680 
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# remove the bad sample (sample # 26) from the metadata
meta_micro_sub <- meta_micro[-26,]
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM-rooted-tree/')
COM_rooted_tree <- ape::read.tree("tree.nwk")
## make a phyloseq object of the asv table, taxonomy table, metadata
# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro_sub), colnames(com.asv.sub))
com.asv.ord  <- com.asv.sub[ ,re_order]
com.asv.physeq = otu_table(com.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(com.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(com.asv.physeq) <- paste0("S", sample_names(com.asv.physeq))
# phyloseq object of the taxonomy table
com.tax <- column_to_rownames(com.tax, var = "ASVid")
#row.names(com.tax) <- com.tax$ASVid
com.tax.physeq = tax_table(as.matrix(com.tax)) # taxonomy table
# phyloseq object of the metadata
meta_micro_sub$Date <- factor(meta_micro_sub$Date, levels = c("4/28/22", "06/01/2022", "07/05/2022", "7/20/22", "9/13/22"),
                          labels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"))
rownames(meta_micro_sub) <- sample_names(com.asv.physeq)
com.meta.physeq <- sample_data(meta_micro_sub)# meta data
sample_names(com.meta.physeq)
# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM-rooted-tree/')
COM_rooted_tree <- ape::read.tree("tree.nwk")
# make phyloseq object
com.physeq <- merge_phyloseq(com.asv.physeq,com.tax.physeq,com.meta.physeq,COM_rooted_tree)
com.physeq # 653 taxa 191 sample
sample_data(com.physeq)$SampleID <- paste0("S", sample_data(com.physeq)$SampleID)
sample_data(com.physeq)
#remove the last sampling date
com.physeq.no.last = subset_samples(com.physeq, Date != "Sept 13th")
View(sample_data(com.physeq.no.last))
com.physeq.no.last # 653 taxa 167 samples BS & RS
otu_table(com.physeq.no.last)
sort(taxa_sums(com.physeq.no.last), decreasing =F) # contain ASVs that don't exist in at least 1 sample
com.physeq.no.last1 <- prune_taxa(taxa_sums(com.physeq.no.last)>0, com.physeq.no.last)
com.physeq.no.last1 # 622 taxa & 167 samples BS & RS
# rarefy to minimum sequencing depth (minimum reads = 3832 reads)
sort(colSums(otu_table(com.physeq.no.last1), na.rm = FALSE, dims = 1), decreasing = F)
set.seed(333)
com.rare.nolast <- rarefy_even_depth(com.physeq.no.last1, sample.size = 5242,
  rngseed = 333, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
com.rare.nolast ## 20 OTUs were removed because they are no longer present in any sample after random subsampling
# 1 sample removed (S52)
sort(sample_sums(com.rare.nolast), decreasing = F) 
sort(rowSums(otu_table(com.rare.nolast), na.rm = FALSE, dims = 1), decreasing = F)
#### Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 
#### Comammox alpha diversity using rarefied data ###
colSums(otu_table(com.rare.nolast))
com.asv.rare.df <- as.data.frame(otu_table(com.rare.nolast))
dim(com.asv.rare.df) #  602 ASVs 166 Samples
com.asv.rare.df_pa <- 1*(com.asv.rare.df>0)
com.s2 <- specnumber(com.asv.rare.df, MARGIN = 2) # richness
com.richness2 <- as.data.frame(com.s2) 
com.h2 <- diversity(t(com.asv.rare.df), index = 'shannon') # Shannon index
com.shannon2 <- as.data.frame(com.h2)
com.d2 <- diversity(t(com.asv.rare.df), index = 'simpson') # Simpson index
com.simpson2 <- as.data.frame(com.d2)
com.inv.d2 <- diversity(t(com.asv.rare.df), index = 'invsimpson')
# edit env mapping data of Comammox
com.meta.df2 <- data.frame(meta_micro_sub)
head(com.meta.df2)
dim(com.meta.df2) #192 Samples
# filter out S52 from the metadata and the last sampling date
com.meta.df.sub2 <- com.meta.df2 %>% filter(SampleID != "52",Date != "Sept 13th")
dim(com.meta.df.sub2) # 166 Samples
com.meta.df.sub2$Date <- factor(com.meta.df.sub2$Date)
# adding alpha diversity calculation
com.meta.df.sub2$COM_Richness <- com.s2
com.meta.df.sub2$COM_Shannon <- com.h2
com.meta.df.sub2$COM_Simpson <- com.d2
com.meta.df.sub2$COM_InvSimpson <- com.inv.d2
# tidy up
str(com.meta.df.sub2)
com.meta.df.sub2$Type <- factor(com.meta.df.sub2$Type, levels = c("BS", "RS"),
                               labels = c("Bulk_Soil", "Rhizosphere"))
com.meta.df.sub2$Treatment <- factor(com.meta.df.sub2$Treatment, levels = c("D", "K", "M"),
                                    labels = c("BIODYN", "CONFYM", "CONMIN"))
com.meta.df.sub2$SampleID<-as.factor(com.meta.df.sub2$SampleID)
com.meta.df.sub2$PlotID<-as.factor(com.meta.df.sub2$PlotID)
com.meta.df.sub2$Irrigation<-as.factor(com.meta.df.sub2$Irrigation)
com.meta.df.sub2$Block<-as.factor(com.meta.df.sub2$Block)
com.meta.df.sub2$x<-as.factor(com.meta.df.sub2$x)
com.meta.df.sub2$var<-as.factor(com.meta.df.sub2$var)
com.meta.df.sub2$var2<-as.factor(com.meta.df.sub2$var2)
com.meta.df.sub2$var3<-as.factor(com.meta.df.sub2$var3)
com.meta.df.sub2[sapply(com.meta.df.sub2, is.character)] <- 
 lapply(com.meta.df.sub2[sapply(com.meta.df.sub2, is.character)], as.numeric)
com.meta.df.sub2[sapply(com.meta.df.sub2, is.integer)] <- 
 lapply(com.meta.df.sub2[sapply(com.meta.df.sub2, is.integer)], as.numeric)

###########################################################################
# 1. Response variable: Alpha Diversity
###########################################################################
install.packages("datarium")
install.packages("rstatix")
library(datarium)
library(rstatix)

# 1a. Analyses of Bulk Soil
com.meta.bulk2 <- com.meta.df.sub2[1:94,]
str(com.meta.bulk2)
dim(com.meta.bulk2)
com.BS.sum.rich <- com.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(COM_Richness, type = "mean_sd")
View(com.BS.sum.rich)
com.BS.rich.plot <- ggboxplot(
  com.meta.bulk2, x = "Irrigation", y = "COM_Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.BS.rich.plot
# check assumption (outliers)
com.BS.rich.out <- com.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(COM_Richness) # no extreme outliers
View(com.BS.rich.out)
# Saphiro-Wilk for normality
com.BS.rich.SW <- com.meta.bulk2 %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(COM_Richness)
View(com.BS.rich.SW)
ggqqplot(com.meta.bulk2, "COM_Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.BS.rich.Lave <- com.meta.bulk2 %>%
  group_by(Date) %>%
  levene_test(COM_Richness ~ Irrigation*Treatment)
View(com.BS.rich.Lave)
# interaction plot
with(com.meta.bulk2, interaction.plot(x.factor = Treatment, 
                               trace.factor = Irrigation, 
                               response = COM_Richness))
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
com.BS.rich.aov <- anova_test(
  data = com.meta.bulk2, type=3, dv = COM_Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.BS.rich.aov)

# with afex package
library(afex)
com.BS.rich.afex <- aov_ez("PlotID", "COM_Richness", com.meta.bulk2, 
             within = "Date",
             between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
com.BS.rich.afex
# Test Method 3 - Richness
com.BS.rich.mod1 <- lme(COM_Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=com.meta.bulk2, na.action = na.omit)
anova(com.BS.rich.mod1) 
com.BS.rich.mod2 <- lmerTest::lmer(com.meta.bulk2$COM_Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk2)
anova(com.BS.rich.mod2)
com.BS.rich.mod3 <- lmerTest::lmer(com.meta.bulk2$COM_Richness ~ Irrigation*Treatment*Date+(1|Block:Date),data=com.meta.bulk2)#control = lmerControl(check.nobs.vs.nRE="ignore"))
anova(com.BS.rich.mod3, type = 3)
# Assumption
shapiro.test(resid(com.BS.rich.mod3)) # normal
plot(simulateResiduals(com.BS.rich.mod3)) # okay
# Test Method 3 - Shannon Index
com.BS.sha.mod1 <- lme(COM_Shannon ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=com.meta.bulk2, na.action = na.omit)
anova(com.BS.sha.mod1) 
com.BS.sha.mod2 <- lmerTest::lmer(com.meta.bulk2$COM_Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk2, na.action = na.omit)
anova(com.BS.sha.mod2)
com.BS.sha.mod3 <- lmerTest::lmer(com.meta.bulk2$COM_Shannon ~ Irrigation*Treatment*Date+(Date|PlotID), data=com.meta.bulk2, na.action=na.omit)
anova(com.BS.sha.mod3)
# Assumption
library(DHARMa)
shapiro.test(resid(com.BS.sha.mod3)) # not normal
plot(simulateResiduals(com.BS.sha.mod3)) # okay
library(afex)
com.BS.sha.afex <- aov_ez("PlotID", "COM_Shannon", com.meta.bulk2, 
             within = "Date",
             between = c("Irrigation","Treatment"),
             type = 3,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
com.BS.sha.afex



