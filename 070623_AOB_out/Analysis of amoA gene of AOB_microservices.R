#################################### Analysis of amoA gene of AOB Illumina MiSeq Data #####################################
##
# Date : 30 May 2023
# Author : Ari Fina BINTARTI

# Install packages
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
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
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
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

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.ASV-analysis')
wd <- print(getwd())
# load the asv table
aob.asv <- read.table('annotated.AOB.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(aob.asv) # 1338  192
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
aob.tax <- read.csv("besthit.diamond.output.curateddb.AOB.ASVs.csv")
dim(aob.tax) # 1338
# load the metadata
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.Phylogenetic-analysis/')
aob.tre <- ape::read.tree("tree.AOB.nwk")

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(aob.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")
library(phyloseq)

# make a phyloseq asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(aob.meta), colnames(aob.asv))
aob.asv.ord  <- aob.asv[ ,re_order]
aob.asv.physeq = otu_table(aob.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aob.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aob.asv.physeq) <- paste0("S", sample_names(aob.asv.physeq))

# phyloseq object of the taxonomy table
aob.tax <- column_to_rownames(aob.tax, var = "ASVid")
aob.tax.physeq = tax_table(as.matrix(aob.tax)) # taxonomy table
 
# phyloseq object of the metadata
column_to_rownames(aob.meta, var = "SampleID")

rownames(aob.meta) <- aob.meta$SampleID
aob.meta.physeq <- sample_data(aob.meta)# meta data

aob.meta$SampleID <- as.factor(aob.meta$SampleID)
aob.meta.df <- aob.meta %>% select(-SampleID) %>% as.data.frame
rownames(aob.meta.df) <- aob.meta$SampleID

sample.data <- sample_data(aob.meta.df)





# make phyloseq object
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq, aob.meta.physeq)




row.names(aob.meta.physeq)










# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 

# 1. AOB alpha diversity for each sample 
aob.asv_pa <- 1*(aob.asv>0)
aob.s <- specnumber(aob.asv, MARGIN = 2) # richness
aob.richness <- as.data.frame(aob.s) 
aob.h <- diversity(t(aob.asv), index = 'shannon') # Shannon index
aob.shannon <- as.data.frame(aob.h)
aob.pielou <- aob.h/log(aob.s) # Pielou's evenness
aob.evenness <- as.data.frame(aob.pielou)
aob.meta <- data.frame(meta_micro) # make data frame of the map data
aob.meta$Richness <- aob.s
aob.meta$Shannon <- aob.h
aob.meta$Pielou <- aob.pielou

# Statistical Analyses: Alpha Diversity

# 1. Bulk Soil: 
# compare alpha diversity between control and shelter, and among managements systems within treatment
aob.meta.bulk <- aob.meta[1:120,] # select only bulk soil sample
aob.meta.bulk$Irrigation<-as.factor(aob.meta.bulk$Irrigation)
aob.meta.bulk$Treatment<-as.factor(aob.meta.bulk$Treatment)
# Richness: fit nested ANOVA
set.seed(13)
aob.nest.rich <- aov(aob.meta.bulk$Richness ~ aob.meta.bulk$Treatment / aob.meta.bulk$Irrigation)
summary(aob.nest.rich) # (p-value Irri & Treat = 0.73 & 0.95, F-val Irri & Treat = 0.42 & 0.047), no significant effect of irrigation and management systems to AOB richness
# double check the effect of irrigation on AOB richness
set.seed(13)
aob.rich.trt <- lm(aob.meta.bulk$Richness ~ Treatment, data=aob.meta.bulk, na.action=na.exclude)
aob.rich.trt
drop1(aob.rich.trt,~.,test="F") # not significant, similar p-val and F-val
# Shannon: fit nested ANOVA
set.seed(13)
aob.nest.sha <- aov(aob.meta.bulk$Shannon ~ aob.meta.bulk$Treatment / aob.meta.bulk$Irrigation)
summary(aob.nest.sha) # (p-value Irri & Treat = 0.35 & 0.95, F-val Irri & Treat = 1.1 & 0.048), no significant effect of irrigation and management systems to AOB shannon diversity index
# Evenness: fit nested ANOVA
set.seed(13)
aob.nest.pie <- aov(aob.meta.bulk$Pielou ~ aob.meta.bulk$Treatment / aob.meta.bulk$Irrigation)
summary(aob.nest.pie) # (p-value Irri & Treat = 0.78 & 0.71, F-val Irri & Treat = 0.35 & 0.33), no significant effect of irrigation and management systems to AOB evenness.

# 2. Rhizosphere Soil: 
# compare alpha diversity between control and shelter, and among managements systems within treatment
aob.meta.rhizo <- aob.meta[121:192,] # select only rhizosphere soil sample
aob.meta.rhizo$Irrigation<-as.factor(aob.meta.rhizo$Irrigation)
aob.meta.rhizo$Treatment<-as.factor(aob.meta.rhizo$Treatment)
# Richness: fit nested ANOVA
set.seed(13)
aob.nest.rich.rz <- aov(aob.meta.rhizo$Richness ~ aob.meta.rhizo$Treatment / aob.meta.rhizo$Irrigation)
summary(aob.nest.rich.rz) # (p-value Irri & Treat = 0.21 & 0.66, F-val Irri & Treat = 1.58 & 0.59), no significant effect of irrigation and management systems to AOB richness
# Shannon: fit nested ANOVA
set.seed(13)
aob.nest.sha.rz <- aov(aob.meta.rhizo$Shannon ~ aob.meta.rhizo$Treatment / aob.meta.rhizo$Irrigation)
summary(aob.nest.sha.rz) # (p-value Irri & Treat = 0.27 & 0.12, F-val Irri & Treat = 1.2 & 1.85), no significant effect of irrigation and management systems to AOB shannon diversity index
# Evenness: fit nested ANOVA
set.seed(13)
aob.nest.pie.rz <- aov(aob.meta.rhizo$Pielou ~ aob.meta.rhizo$Irrigation / aob.meta.rhizo$Treatment)
summary(aob.nest.pie.rz) # (p-value Irri & Treat = 0.38 & 0.29, F-val Irri & Treat = 0.75 & 1.25), no significant effect of irrigation and management systems to AOB evenness.

### Conclusion: Irrigation and management systems have no effect on the AOB Alpha diversity ###






