#################################### Analysis of amoA gene of AOA Illumina MiSeq Data #####################################
##
# Date : 08 June 2023
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

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA.ASV-analysis')
wd <- print(getwd())
# load the asv table
aoa.asv <- read.table('annotated.AOA.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(aoa.asv) # 646  192
sort(colSums(aoa.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/')
aoa.tax <- read.csv("besthit.diamond.output.curateddb.AOA.ASVs.csv")
dim(aoa.tax) # 646
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA.Phylogenetic-analysis/')
aoa.tre <- ape::read.tree("tree.AOA.nwk")

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(aoa.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")
library(phyloseq)

## make a phyloseq object of the asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aoa.asv))
aoa.asv.ord  <- aoa.asv[ ,re_order]
aoa.asv.physeq = otu_table(aoa.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aoa.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aoa.asv.physeq) <- paste0("S", sample_names(aoa.asv.physeq))

# phyloseq object of the taxonomy table
aoa.tax <- column_to_rownames(aoa.tax, var = "ASVid")
aoa.tax.physeq = tax_table(as.matrix(aoa.tax)) # taxonomy table
 
# phyloseq object of the metadata
rownames(meta_micro) <- sample_names(aoa.asv.physeq)
aoa.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aoa.meta.physeq)

# make phyloseq object
aoa.physeq <- merge_phyloseq(aoa.asv.physeq,aoa.tax.physeq,aoa.meta.physeq)
aoa.physeq

# run the ggrare function attached in the file "generating_rarecurve.r"
aoa.rare <- ggrare(aoa.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aoa.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aoa.rare <- aoa.rare + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(b) AOA", )+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "right",
        legend.title = element_text(size=15),
        legend.text = element_text(size = 13),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Number of ASVs")+xlab("Reads")
 
plot.aoa.rare

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')

ggsave("AOA_rarecurve.jpg",
       plot.aoa.rare, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

sort(sample_sums(aoa.physeq), decreasing = F)

# subset samples with sample sum less than 5000 reads
aoa.physeq5k <- prune_samples(sample_sums(aoa.physeq) < 5000, aoa.physeq)
aoa.physeq5k
sort(sample_sums(aoa.physeq5k), decreasing = F)

# rarefy to minimum sequencing depth
set.seed(13)
aoa.rare.min.physeq <- rarefy_even_depth(aoa.physeq, sample.size = min(sample_sums(aoa.physeq)),
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min), decreasing = F) # 70 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed

# run the ggrare function attached in the file "generating_rarecurve.r"
aoa.rare.min <- ggrare(aoa.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aoa.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aoa.rare.min <- aoa.rare.min + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(b) AOA", )+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "right",
        legend.title = element_text(size=15),
        legend.text = element_text(size = 13),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Number of ASVs")+xlab("Reads")
 
plot.aoa.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')

ggsave("AOA_rarecurve_min.jpg",
       plot.aoa.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)


