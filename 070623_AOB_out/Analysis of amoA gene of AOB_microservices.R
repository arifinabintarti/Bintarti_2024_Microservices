#############################################################################################
# Analysis of amoA gene of AOB Illumina MiSeq Data 
#############################################################################################

# Date : 30 May 2023
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
setwd('/Users/arifinabintarti/Documents/France/microservices/')
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
meta_micro$Date <- factor(meta_micro$Date, levels = c("4/28/22", "6/1/22", "7/5/22", "7/20/22", "9/13/22"),
                  labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))
rownames(meta_micro) <- sample_names(aob.asv.physeq)
aob.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aob.meta.physeq)

# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB-rooted-tree/')
AOB_rooted_tree <- ape::read.tree("tree.nwk")

# make phyloseq object
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq,AOB_rooted_tree)
aob.physeq
sample_data(aob.physeq)$SampleID <- paste0("S", sample_data(aob.physeq)$SampleID)
sample_data(aob.physeq)

# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
aob.rare <- ggrare(aob.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare <- aob.rare + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
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
 
plot.aob.rare

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve.jpg",
       plot.aob.rare, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

sort(sample_sums(aob.physeq), decreasing = F)


# rarefy to minimum sequencing depth
set.seed(13)
aob.rare.min.physeq <- rarefy_even_depth(aob.physeq, sample.size = min(sample_sums(aob.physeq)),
  rngseed =13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min.physeq), decreasing = F) # 172 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed

# run the ggrare function attached in the file "generating_rarecurve.r"
#aob.rare.min <- ggrare(aob.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)
test <- ggrare(aob.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare.min <- aob.rare.min + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
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
 
plot.aob.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve_min.jpg",
       plot.aob.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

# rarefy to 2k of sequencing depth
set.seed(13)
aob.rare.2k.physeq <- rarefy_even_depth(aob.physeq, sample.size = 2000,
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.2k.physeq
sort(sample_sums(aob.rare.2k.physeq), decreasing = F) # 102 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # 5 samples removed
# subset samples with sample sum less than 2000 reads
aob.physeq2k <- prune_samples(sample_sums(aob.physeq) < 2000, aob.physeq)
sort(sample_sums(aob.physeq2k), decreasing = F)

######################################################################################################

# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 

#### AOB alpha diversity using rarefied data to the minimum sequencing depth (reads=858) ###

colSums(otu_table(aob.rare.min.physeq))
aob.asv.min <- as.data.frame(otu_table(aob.rare.min.physeq))
dim(aob.asv.min) # 1166 ASVs
aob.asv.min_pa <- 1*(aob.asv.min>0)
aob.min.s <- specnumber(aob.asv.min, MARGIN = 2) # richness
aob.min.richness <- as.data.frame(aob.min.s) 
aob.min.h <- diversity(t(aob.asv.min), index = 'shannon') # Shannon index
aob.min.shannon <- as.data.frame(aob.min.h)
aob.min.pielou <- aob.min.h/log(aob.min.s) # Pielou's evenness
aob.min.evenness <- as.data.frame(aob.min.pielou)
aob.min.d <- diversity(t(aob.asv.min), index = 'simpson') # Simpson index
aob.min.simpson <- as.data.frame(aob.min.d)
aob.min.inv.d <- diversity(t(aob.asv.min), index = 'invsimpson')
#aob.min.meta <- data.frame(meta_micro) # make data frame of the map data
#aob.min.meta$Richness <- aob.min.s
#aob.min.meta$Shannon <- aob.min.h
#aob.min.meta$Pielou <- aob.min.pielou
#aob.min.meta$Simpson <- aob.min.d
#aob.min.meta$InvSimpson <- aob.min.inv.d
# Statistical Analyses: Alpha Diversity
#aob.min.meta$SampleID <- as.factor(aob.min.meta$SampleID)
#aob.min.meta$PlotID <- as.factor(aob.min.meta$PlotID)
#aob.min.meta$Irrigation<-as.factor(aob.min.meta$Irrigation)
#aob.min.meta$Treatment<-as.factor(aob.min.meta$Treatment)
#aob.min.meta$Type<-as.factor(aob.min.meta$Type)

# 1. Richness

# Line plot of AOB richness 
aob.min.meta.df <- data.frame(meta_micro)
aob.min.meta.df$Richness <- aob.min.s
aob.min.meta.df$Shannon <- aob.min.h
aob.min.meta.df$Pielou <- aob.min.pielou
aob.min.meta.df$Simpson <- aob.min.d
aob.min.meta.df$InvSimpson <- aob.min.inv.d
#aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
str(aob.min.meta.df)
aob.min.meta.df$Date <- factor(aob.min.meta.df$Date, levels = c("4/28/22", "6/1/22", "7/5/22", "7/20/22", "9/13/22"),
                  labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))
aob.min.meta.df$Type <- factor(aob.min.meta.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere"))
aob.min.meta.df$Treatment <- factor(aob.min.meta.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
aob.min.meta.df$SampleID<-as.factor(aob.min.meta.df$SampleID)
aob.min.meta.df$PlotID<-as.factor(aob.min.meta$PlotID)
aob.min.meta.df$Irrigation<-as.factor(aob.min.meta$Irrigation)
# tidy up the data frame
aob.min.meta.df.tidy <- aob.min.meta.df %>%
                             group_by(Irrigation, Treatment, Date,  Type, var2,var3) %>%
                             summarize(Mean.Rich=mean(Richness),
                                       Mean.Sha=mean(Shannon),
                                       Mean.Simp=mean(Simpson),
                                       Mean.invsimp=mean(InvSimpson))
str(aob.min.meta.df.tidy)
#setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
#write.csv(aob.min.meta.df.tidy, file = "aob.min.meta.df.tidy2.csv")
#aob.min.meta.df.tidy.ed <- read.csv("aob.min.meta.df.tidy.csv")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
display_carto_pal(7, "Vivid")
carto_pal(n = NULL, 'Vivid')
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
library(ggnewscale)
aob.min.meta.df$Date <- factor(aob.min.meta.df$Date, levels = unique(aob.min.meta.df$Date))
aob.min.rich.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
                             geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Richness")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
aob.min.rich.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich3.eps",
       aob.min.rich.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

dev.off()

# 2. Shannon

#Line plot of AOB Shannon
library(scales)
aob.min.sha.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
                             geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Shannon")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             
aob.min.sha.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha2.eps",
       aob.min.sha.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

# 3. Simpson

#Line plot of AOB Simpson
aob.min.simp.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Simp, linetype=Irrigation))+
                             geom_line(linewidth=1, aes(group = var2, col=Treatment))+
                             #ylim(0.9,0.97)+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             #expand_limits(y = 0.9)+
                             theme_bw() +
                             scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Simpson")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             #new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             #scale_color_manual(values = color.trt)
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             #scale_colour_manual(values=c('black', 'red'))
                             #scale_y_continuous(labels = label_number(accuracy = 0.01))
                             
aob.min.simp.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_simp.eps",
       aob.min.simp.plot, device = cairo_ps,
       width = 8, height = 7, 
       units= "in", dpi = 600)

#. 4. Inverse Simpson

#Line plot of AOB Inverse Simpson
aob.min.invsimp.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
                             geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Inverse Simpson")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             
aob.min.invsimp.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp.eps",
       aob.min.invsimp.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#############################################################################

# 1. Richness

# Box plot of AOB Richness

#install.packages("ggpattern")
#install.packages("sf")
library(sf)
library(ggpattern)
#soiltype <- c("Bulk Soil", "Rhizosphere")
#names(soiltype) <- c("BS", "RS")
#aob.min.meta.df$Type <- factor(aob.min.meta.df$Type, levels = c("BS", "RS"),
                  #labels = c("Bulk Soil", "Rhizosphere"))
#aob.min.meta.df$Treatment <- factor(aob.min.meta.df$Treatment, levels = c("D", "K", "M"),
                  #labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
#aob.min.rich.boxplot <- ggplot(aob.min.meta.df, aes(x=Treatment, y=Richness, pattern = Irrigation, fill = Treatment)) +
               #geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               #theme_bw() +
               #labs(x="Treatment", y="AOB Richness")+
               #expand_limits(y = 0)+
               #labs(pattern="Irrigation")+
               #scale_colour_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               #scale_pattern_manual(values = c(Rainout = "stripe", Control = "none"))+
               #facet_grid(Type~Date)+
               #geom_boxplot_pattern(position = position_dodge(preserve = "single"), 
                     #color = "black", pattern_fill = "black", 
                     #pattern_angle = 45, pattern_density = 0.05, 
                     #pattern_spacing = 0.018, 
                     #pattern_key_scale_factor = 0.6)+
                     #theme(legend.title = element_text(size=15, face='bold'),
                     #legend.text = element_text(size=15),
                     #strip.text = element_text(size=18),
                     #axis.text.y = element_text(size = 18),
                     #axis.title.y = element_text(size=18,face="bold"),
                     #axis.text.x = element_blank(),
                     #axis.title.x =element_blank(),
                     #axis.ticks.x = element_blank(),
                     #plot.background = element_blank(),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank())+
                     #guides(pattern = guide_legend(override.aes = list(fill = "white")),
                     #fill = guide_legend(override.aes = list(pattern = "none")))
#aob.min.rich.boxplot
   
#setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
#ggsave("AOB_min_rich_box.eps",
       #aob.min.rich.boxplot, device = "eps",
       #width = 12.5, height = 7, 
       #units= "in", dpi = 600)

str(aob.min.meta.df)
aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%y")
# Richness: plotting the significance across treatment
aob.min.rich.pwc.plot <- ggplot(aob.min.meta.df, aes(x=Irrigation, y=Richness)) +
               geom_boxplot(aes(fill=Treatment))+
               #stat_summary(fun.y=mean, geom="point", aes(group=Treatment), position=position_dodge(0.75), 
               #color="red", size=1.5)+
               #geom_point(aob.min.meta.df.tidy, mapping=aes(y = Mean.Rich, x = Irrigation), shape = 20, size = 3, color = 'red')+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Richness")+
               #ylim(0,130)+
               expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               facet_grid(Type~ Date,scales="free_x")+
                     theme(legend.title = element_text(size=15, face='bold'),
                     legend.text = element_text(size=15),
                     strip.text = element_text(size=18),
                     axis.text = element_text(size = 18),
                     axis.title.y = element_text(size=18,face="bold"),
                     axis.title.x =element_blank(),
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
aob.min.rich.pwc.plot

# adding treatment significance from the stats analyses
aob.rich.pwc.all <- aob.min.meta.df %>%
                    group_by(Type, Date, Irrigation) %>%
                    pairwise_t_test(Richness ~ Treatment, p.adjust.method = "BH") %>%
                    select(-p, -p.signif)
# plotting the treatment significances
aob.rich.pwc.all <- aob.rich.pwc.all %>% 
                    add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.rich.pwc.plot2 <- aob.min.rich.pwc.plot + 
                          stat_pvalue_manual(aob.rich.pwc.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
                          scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.rich.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich_mean_boxplot.eps",
       aob.min.rich.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# richness between irrigations
aob.rich.pwc.irri.plot <- ggplot(aob.min.meta.df, aes(x=Date, y=Richness)) +
                          geom_boxplot(aes(group = var3, fill = Irrigation))+
                          theme_bw() +
                          labs(y="AOB Richness")+
                          scale_fill_manual(values = c("#e5f5e0","#31a354"))+
                          #scale_x_discrete(labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))+
                          facet_grid(Type~ Treatment,scales="free_x")+
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
aob.rich.pwc.irri.plot
# adding irrigation significance from the stats analyses
aob.rich.pwc.irri <- aob.min.meta.df %>%
                    group_by(Type, Date, Treatment) %>%
                    pairwise_t_test(Richness ~ Irrigation, p.adjust.method = "BH") %>%
                    select(-p, -p.signif) #not significant
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_rich_irri_boxplot.eps",
       aob.rich.pwc.irri.plot, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

# 2. Shannon

# Shannon: plotting the significance across treatment

aob.min.sha.pwc.plot <- ggplot(aob.min.meta.df, aes(x=Irrigation, y=Shannon)) +
               geom_boxplot(aes(fill = Treatment))+
               #stat_summary(fun.y=mean, geom="point", aes(group=Treatment), position=position_dodge(0.75), 
               #color="red", size=1.5)+
               theme_bw() +
               labs(y="AOB Shannon")+
               #expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               facet_grid(Type~ Date,scales="free_x")+
                     theme(legend.title = element_text(size=15, face='bold'),
                     legend.text = element_text(size=15),
                     strip.text = element_text(size=18),
                     axis.text = element_text(size = 18),
                     axis.title.y = element_text(size=18,face="bold"),
                     axis.title.x =element_blank(),
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
aob.min.sha.pwc.plot
# adding asterix from the stats analyses

aob.sha.pwc.all <- aob.min.meta.df %>%
                   group_by(Type, Date, Irrigation) %>%
                   pairwise_t_test(Shannon ~ Treatment, p.adjust.method = "BH") %>%
                   select(-p, -p.signif)
# plotting the significance
aob.sha.pwc.all <- aob.sha.pwc.all %>%  
                   add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.sha.pwc.plot2 <- aob.min.sha.pwc.plot + 
                         stat_pvalue_manual(aob.sha.pwc.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
                         scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.sha.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha_boxplot.eps",
       aob.min.sha.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# Shannon: plotting the significance across irrigation
#aob.sha.df <- aob.min.meta.df %>% 
    #group_by(Date, Irrigation, Treatment, Type, var3) %>% 
    #summarise(y0 = min(Shannon), 
        #y25 = quantile(Shannon, 0.25), 
        #y50 = mean(Shannon), 
        #y75 = quantile(Shannon, 0.75), 
        #y100 = max(Shannon))
#str(aob.sha.df)
#aob.sha.df$Date  <- as.Date(aob.sha.df$Date , "%m/%d/%y")
str(aob.min.meta.df)
aob.sha.pwc.irri.plot <- ggplot(aob.min.meta.df, aes(x=Date, y=Shannon)) +
                geom_boxplot(aes(group = var3, fill = Irrigation))+
               #geom_boxplot(aes(group = var3, ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
               #stat = "identity")+
               #ggplot(aob.sha.df, aes(x=Date, fill=Irrigation)) +
               #stat_summary(fun.y=mean, geom="point", aes(group=Irrigation), position=position_dodge(0.45), 
               #color="red", size=1.5)+
               theme_bw() +
               labs(y="AOB Shannon")+
               #expand_limits(y = 0)+
               #labs(pattern="Irrigation")+
               scale_fill_manual(values = c("#e5f5e0","#31a354"))+
               scale_x_discrete(labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))+
               facet_grid(Type~ Treatment,scales="free_x")+
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
aob.sha.pwc.irri.plot
# adding irrigation significance from the stats analyses
aob.sha.pwc.irri <- aob.min.meta.df %>%
                    group_by(Type, Date, Treatment) %>%
                    pairwise_t_test(Shannon ~ Irrigation, p.adjust.method = "BH") %>%
                    select(-p, -p.signif)
# plotting the irrigation significance
aob.sha.pwc.irri <- aob.sha.pwc.irri %>% 
                    add_xy_position(x = "Date", dodge = 0.5)
aob.sha.pwc.irri.plot2 <- aob.sha.pwc.irri.plot + 
                          stat_pvalue_manual(aob.sha.pwc.irri,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0.8, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
                          scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.sha.pwc.irri.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_sha_irri_boxplot.eps",
       aob.sha.pwc.irri.plot2, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

# Simpson

# Simpson: plotting the significance across treatment
aob.min.simp.pwc.plot <- ggplot(aob.min.meta, aes(x=Irrigation, y=Simpson)) +
               geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Simpson")+
               #expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               facet_grid(Type~ Date,scales="free_x")+
                     theme(legend.title = element_text(size=15, face='bold'),
                     legend.text = element_text(size=15),
                     strip.text = element_text(size=18),
                     axis.text = element_text(size = 18),
                     axis.title.y = element_text(size=18,face="bold"),
                     axis.title.x =element_blank(),
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
aob.min.simp.pwc.plot
# adding asterix from the stats analyses
aob.simp.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(Simpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.simp.pwc.all <- aob.simp.pwc.all %>% add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.simp.pwc.plot2 <- aob.min.simp.pwc.plot + 
 stat_pvalue_manual(aob.simp.pwc.all,label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.simp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_simp_all.eps",
       aob.min.simp.pwc.plot2, device = "eps",
       width = 16, height =7, 
       units= "in", dpi = 600)


# Inverse Simpson

# Inverse Simpson: plotting the significance across treatment
aob.min.invsimp.pwc.plot <- ggplot(aob.min.meta.df, aes(x=Irrigation, y=InvSimpson)) +
               geom_boxplot(aes(fill = Treatment))+
               theme_bw() +
               labs(y="AOB Inverse Simpson")+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               facet_grid(Type~ Date,scales="free_x")+
                     theme(legend.title = element_text(size=15, face='bold'),
                     legend.text = element_text(size=15),
                     strip.text = element_text(size=18),
                     axis.text = element_text(size = 18),
                     axis.title.y = element_text(size=18,face="bold"),
                     axis.title.x =element_blank(),
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
aob.min.invsimp.pwc.plot
# adding asterix from the stats analyses
aob.invsimp.pwc.all <- aob.min.meta.df %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(InvSimpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.invsimp.pwc.all <- aob.invsimp.pwc.all %>% 
                       add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.invsimp.pwc.plot2 <- aob.min.invsimp.pwc.plot + 
                             stat_pvalue_manual(aob.invsimp.pwc.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
                             scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.invsimp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp_all.eps",
       aob.min.invsimp.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# inverse simpson between irrigations
aob.invsimp.pwc.irri.plot <- ggplot(aob.min.meta.df, aes(x=Date, y=InvSimpson)) +
                             geom_boxplot(aes(group = var3, fill = Irrigation))+
                             theme_bw() +
                             labs(y="AOB Inverse Simpson")+
                             scale_fill_manual(values = c("#e5f5e0","#31a354"))+
                             #scale_x_discrete(labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))+
                             facet_grid(Type~ Treatment,scales="free_x")+
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
aob.invsimp.pwc.irri.plot
# adding irrigation significance from the stats analyses
aob.invsimp.pwc.irri <- aob.min.meta.df %>%
                    group_by(Type, Date, Treatment) %>%
                    pairwise_t_test(InvSimpson ~ Irrigation, p.adjust.method = "BH") %>%
                    select(-p, -p.signif)
# plotting the irrigation significance
aob.invsimp.pwc.irri <- aob.invsimp.pwc.irri %>% 
                        add_xy_position(x = "Date", dodge = 0.5)
aob.invsimp.pwc.irri.plot2 <- aob.invsimp.pwc.irri.plot + 
                              stat_pvalue_manual(aob.invsimp.pwc.irri,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0.8, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
                              scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.invsimp.pwc.irri.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_invsimp_irri_boxplot.eps",
       aob.invsimp.pwc.irri.plot2, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

###################################################################################
# Beta Diversity Analyses on Rarefied Data: AOB
###################################################################################

# FOR ALL SAMPLES
# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

#aob.asv.min_PA <- 1*(aob.asv.min>0)
#aob.asv.min_PA
# Bray-Curtis using rarefied data:
aob.asv.min_dist <- vegdist(t(aob.asv.min), method = "bray")
# jaccard using rarefied data:
aob.asv.min_dist_jac <- vegdist(t(aob.asv.min), binary = TRUE, method = "jaccard") 
# Weighted UniFrac using rarefied data:
aob.rare_tree = phy_tree(aob.rare.min.physeq)
sprintf("Is tree binary: %s", is.binary(aob.rare_tree)) ## Check if binary (dichotomy) & multifurcating (polytomy) trees (please check: https://github.com/joey711/phyloseq/issues/1643)
phy_tree(aob.rare.min.physeq) = multi2di(aob.rare_tree)#If FASLE, randomly resolve polytomies and replace tree in "aob.rare.min.physeq"
sprintf("Is tree binary: %s", is.binary(phy_tree(aob.rare.min.physeq)))
aob.wUF.rare_dist <- UniFrac(aob.rare.min.physeq, weighted=TRUE, normalized = TRUE)
aob.wUF.rare_dist
# Unweighted UniFrac using rarefied data:
aob.uwUF.rare_dist <- UniFrac(aob.rare.min.physeq, weighted=FALSE, normalized = TRUE)
aob.uwUF.rare_dist

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis using rarefied data:
aob.asv.min_pcoa <- cmdscale(aob.asv.min_dist, eig=T)
# jaccard using rarefied data:
aob.asv.min_pcoa.jac <- cmdscale(aob.asv.min_dist_jac, eig=T)
# Weighted UniFrac using rarefied data:
aob.rare_pcoa.wUF <- cmdscale(aob.wUF.rare_dist, eig=T)
# Unweighted UniFrac using rarefied data:
aob.rare_pcoa.uwUF <- cmdscale(aob.uwUF.rare_dist, eig=T)

# 3. scores of PC1 and PC2

# bray-curtis:
ax1.scores <- aob.asv.min_pcoa$points[,1]
ax2.scores <- aob.asv.min_pcoa$points[,2] 
# jaccard:
ax1.scores.j <- aob.asv.min_pcoa.jac$points[,1]
ax2.scores.j <- aob.asv.min_pcoa.jac$points[,2]
# Weighted UniFrac using rarefied data:
ax1.scores.wUF.rare <- aob.rare_pcoa.wUF$points[,1]
ax2.scores.wUF.rare <- aob.rare_pcoa.wUF$points[,2]
# Unweighted UniFrac using rarefied data:
ax1.scores.uwUF.rare <- aob.rare_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rare <- aob.rare_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot

#Bray-curtis:
ax1 <- aob.asv.min_pcoa$eig[1]/sum(aob.asv.min_pcoa$eig)
ax2 <- aob.asv.min_pcoa$eig[2]/sum(aob.asv.min_pcoa$eig)
aob.map.pcoa <- cbind(aob.min.meta,ax1.scores,ax2.scores)
# jaccard
ax1.j <- aob.asv.min_pcoa.jac$eig[1]/sum(aob.asv.min_pcoa.jac$eig)
ax2.j <- aob.asv.min_pcoa.jac$eig[2]/sum(aob.asv.min_pcoa.jac$eig)
aob.map.pcoa.j <- cbind(aob.min.meta,ax1.scores.j,ax2.scores.j)
# Weighted UniFrac using rarefied:
ax1.wUF.rare <- aob.rare_pcoa.wUF$eig[1]/sum(aob.rare_pcoa.wUF$eig)
ax2.wUF.rare <- aob.rare_pcoa.wUF$eig[2]/sum(aob.rare_pcoa.wUF$eig)
aob.map.pcoa.wUF.rare <- cbind(aob.min.meta,ax1.scores.wUF.rare,ax2.scores.wUF.rare)
# unweighted UniFrac
ax1.uwUF <- aob.asv.min_pcoa.uwUF$eig[1]/sum(aob.asv.min_pcoa.uwUF$eig)
ax2.uwUF <- aob.asv.min_pcoa.uwUF$eig[2]/sum(aob.asv.min_pcoa.uwUF$eig)
aob.map.pcoa.uwUF <- cbind(aob.min.meta,ax1.scores.uwUF,ax2.scores.uwUF)

# simple plot
aob.pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
library(viridis)

# a. Bray-Curtis:

set.seed(13)
aob.pcoa_plot <- ggplot(data = aob.map.pcoa, aes(x=ax1.scores, y=ax2.scores))+
            theme_bw()+
            geom_point(data = aob.map.pcoa, aes(x = ax1.scores, y = ax2.scores, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_PCoA_rare.tiff",
       aob.pcoa_plot, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# b. Jaccard:

set.seed(13)
aob.pcoa_plot.jac <- ggplot(data = aob.map.pcoa.j, aes(x=ax1.scores.j, y=ax2.scores.j, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.j, aes(x = ax1.scores.j, y = ax2.scores.j, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.jac
 jaccard
set.seed(13)

# c. Weighted UniFrac using rarefied data:

aob.pcoa_plot.wUF.rare <- ggplot(data = aob.map.pcoa.wUF.rare, aes(x=ax1.scores.wUF.rare, y=ax2.scores.wUF.rare, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.wUF.rare, aes(x = ax1.scores.wUF.rare, y = ax2.scores.wUF.rare , shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.wUF.rare,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.wUF.rare,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Weighted UniFrac (rarefied)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.wUF.rare
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_WeightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.wUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# d. UnWeighted UniFrac using rarefied data:

aob.pcoa_plot.uwUF.rare <- ggplot(data = aob.map.pcoa.uwUF.rare, aes(x=ax1.scores.uwUF.rare, y=ax2.scores.uwUF.rare, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.uwUF.rare, aes(x = ax1.scores.uwUF.rare, y = ax2.scores.uwUF.rare, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.uwUF.rare,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.uwUF.rare,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Unweighted UniFrac (rarefied)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.uwUF.rare
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_UnweightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.uwUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

######################################################################################
# SEPARATE BETWEEN BULK SOIL AND RHIZOSPHERE
######################################################################################

# BULK SOIL

# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# Bray-Curtis - Bulk Soil :
aob.asv.min.bulk <- aob.asv.min[,1:120]
aob.asv.min.bulk1 <- aob.asv.min.bulk[rowSums(aob.asv.min.bulk)>0,]
sort(rowSums(aob.asv.min.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
aob.bulk_dist_bc <- vegdist(t(aob.asv.min.bulk1), method = "bray")
# jaccard - Bulk Soil :
aob.bulk_dist_jac <- vegdist(t(aob.asv.min.bulk1), binary = TRUE, method = "jaccard") 
# Weighted UniFrac (rarefied) - Bulk Soil:
aob.physeq_bulk <- subset_samples(aob.rare.min.physeq, Type=="BS")
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
sort(taxa_sums(aob.physeq_bulk1), decreasing =F) #checking
aob.bulk_dist_wUF <- UniFrac(aob.physeq_bulk1, weighted=TRUE, normalized = TRUE)
aob.bulk_dist_wUF
# Unweighted UniFrac (rarefied) -  Bulk Soil:
aob.bulk_dist_uwUF <- UniFrac(aob.physeq_bulk1, weighted=FALSE, normalized = TRUE)
aob.bulk_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Bulk Soil:
aob.bulk_pcoa_bc <- cmdscale(aob.bulk_dist_bc, eig=T)
# Jaccard - Bulk Soil:
aob.bulk_pcoa_jac <- cmdscale(aob.bulk_dist_jac, eig=T)
# Weighted UniFrac - Bulk Soil:
aob.bulk_pcoa_wUF <- cmdscale(aob.bulk_dist_wUF, eig=T)
# Unweighted UniFrac - Bulk Soil:
aob.bulk_pcoa.uwUF <- cmdscale(aob.bulk_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Bulk Soil:
ax1.scores.bulk <- aob.bulk_pcoa_bc$points[,1]
ax2.scores.bulk <- aob.bulk_pcoa_bc$points[,2] 
# Jaccard - Bulk Soil:
ax1.scores.j.bulk <- aob.bulk_pcoa_jac$points[,1]
ax2.scores.j.bulk <- aob.bulk_pcoa_jac$points[,2]
# Weighted UniFrac - Bulk Soil:
ax1.scores.wUF.bulk <- aob.bulk_pcoa_wUF$points[,1]
ax2.scores.wUF.bulk <- aob.bulk_pcoa_wUF$points[,2]
# Unweighted UniFrac - Bulk Soil:
ax1.scores.uwUF.bulk <- aob.bulk_pcoa.uwUF$points[,1]
ax2.scores.uwUF.bulk <- aob.bulk_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot

# Bray-curtis - Bulk Soil:
ax1.bulk <- aob.bulk_pcoa_bc$eig[1]/sum(aob.bulk_pcoa_bc$eig)
ax2.bulk <- aob.bulk_pcoa_bc$eig[2]/sum(aob.bulk_pcoa_bc$eig)
aob.map.pcoa.bulk <- cbind(aob.min.meta.bulk,ax1.scores.bulk,ax2.scores.bulk)
# Jaccard - Bulk Soil:
ax1.j.bulk <- aob.bulk_pcoa_jac$eig[1]/sum(aob.bulk_pcoa_jac$eig)
ax2.j.bulk <- aob.bulk_pcoa_jac$eig[2]/sum(aob.bulk_pcoa_jac$eig)
aob.map.pcoa.j.bulk <- cbind(aob.min.meta.bulk,ax1.scores.j.bulk,ax2.scores.j.bulk)
# Weighted UniFrac - Bulk Soil:
ax1.wUF.bulk <- aob.bulk_pcoa_wUF$eig[1]/sum(aob.bulk_pcoa_wUF$eig)
ax2.wUF.bulk <- aob.bulk_pcoa_wUF$eig[2]/sum(aob.bulk_pcoa_wUF$eig)
aob.map.pcoa.wUF.bulk <- cbind(aob.min.meta.bulk,ax1.scores.wUF.bulk,ax2.scores.wUF.bulk)
# Unweighted UniFrac - Bulk Soil:
ax1.uwUF.bulk <- aob.bulk_pcoa.uwUF$eig[1]/sum(aob.bulk_pcoa.uwUF$eig)
ax2.uwUF.bulk <- aob.bulk_pcoa.uwUF$eig[2]/sum(aob.bulk_pcoa.uwUF$eig)
aob.map.pcoa.uwUF.bulk <- cbind(aob.min.meta.bulk,ax1.scores.uwUF.bulk,ax2.scores.uwUF.bulk)

#################################################################################################

# RHIZOSPHERE

# Bray-Curtis - Rhizosphere :
aob.asv.min.rh <- aob.asv.min[,121:192]
aob.asv.min.rh1 <- aob.asv.min.rh[rowSums(aob.asv.min.rh)>0,]
sort(rowSums(aob.asv.min.rh1, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(aob.asv.min.rh1) #796
aob.rh_dist_bc <- vegdist(t(aob.asv.min.rh1), method = "bray")
# jaccard - Rhizosphere :
aob.rh_dist_jac <- vegdist(t(aob.asv.min.rh1), binary = TRUE, method = "jaccard") 
# Weighted UniFrac (rarefied) - Rhizosphere :
aob.physeq_rh <- subset_samples(aob.rare.min.physeq, Type=="RS")
aob.physeq_rh1 <- prune_taxa(taxa_sums(aob.physeq_rh)>0, aob.physeq_rh)
aob.physeq_rh1
sort(taxa_sums(aob.physeq_rh1), decreasing =F) #checking
aob.rh_dist_wUF <- UniFrac(aob.physeq_rh1, weighted=TRUE, normalized = TRUE)
aob.rh_dist_wUF
# Unweighted UniFrac (rarefied) -  Rhizosphere :
aob.rh_dist_uwUF <- UniFrac(aob.physeq_rh1, weighted=FALSE, normalized = TRUE)
aob.rh_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Rhizosphere :
aob.rh_pcoa_bc <- cmdscale(aob.rh_dist_bc, eig=T)
# Jaccard - Rhizosphere :
aob.rh_pcoa_jac <- cmdscale(aob.rh_dist_jac, eig=T)
# Weighted UniFrac - Rhizosphere :
aob.rh_pcoa_wUF <- cmdscale(aob.rh_dist_wUF, eig=T)
# Unweighted UniFrac - Rhizosphere :
aob.rh_pcoa.uwUF <- cmdscale(aob.rh_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Rhizosphere :
ax1.scores.rh <- aob.rh_pcoa_bc$points[,1]
ax2.scores.rh <- aob.rh_pcoa_bc$points[,2] 
# Jaccard - Rhizosphere :
ax1.scores.j.rh <- aob.rh_pcoa_jac$points[,1]
ax2.scores.j.rh <- aob.rh_pcoa_jac$points[,2]
# Weighted UniFrac - Rhizosphere :
ax1.scores.wUF.rh <- aob.rh_pcoa_wUF$points[,1]
ax2.scores.wUF.rh <- aob.rh_pcoa_wUF$points[,2]
# Unweighted UniFrac - Rhizosphere :
ax1.scores.uwUF.rh <- aob.rh_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rh <- aob.rh_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot

# Bray-curtis - Rhizosphere :
ax1.rh <- aob.rh_pcoa_bc$eig[1]/sum(aob.rh_pcoa_bc$eig)
ax2.rh <- aob.rh_pcoa_bc$eig[2]/sum(aob.rh_pcoa_bc$eig)
aob.map.pcoa.rh <- cbind(aob.min.meta.rh,ax1.scores.rh,ax2.scores.rh)
# Jaccard - Rhizosphere :
ax1.j.rh <- aob.rh_pcoa_jac$eig[1]/sum(aob.rh_pcoa_jac$eig)
ax2.j.rh <- aob.rh_pcoa_jac$eig[2]/sum(aob.rh_pcoa_jac$eig)
aob.map.pcoa.j.rh <- cbind(aob.min.meta.rh,ax1.scores.j.rh,ax2.scores.j.rh)
# Weighted UniFrac - Rhizosphere :
ax1.wUF.rh <- aob.rh_pcoa_wUF$eig[1]/sum(aob.rh_pcoa_wUF$eig)
ax2.wUF.rh <- aob.rh_pcoa_wUF$eig[2]/sum(aob.rh_pcoa_wUF$eig)
aob.map.pcoa.wUF.rh <- cbind(aob.min.meta.rh,ax1.scores.wUF.rh,ax2.scores.wUF.rh)
# Unweighted UniFrac - Rhizosphere :
ax1.uwUF.rh <- aob.rh_pcoa.uwUF$eig[1]/sum(aob.rh_pcoa.uwUF$eig)
ax2.uwUF.rh <- aob.rh_pcoa.uwUF$eig[2]/sum(aob.rh_pcoa.uwUF$eig)
aob.map.pcoa.uwUF.rh <- cbind(aob.min.meta.rh,ax1.scores.uwUF.rh,ax2.scores.uwUF.rh)

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
library(viridis)

# A. Bray-Curtis - Bulk Soil :
aob.pcoa_bulk.plot <- ggplot(data = aob.map.pcoa.bulk, aes(x=ax1.scores.bulk, y=ax2.scores.bulk, colour=Treatment))+
                             theme_bw()+
                      geom_point(data = aob.map.pcoa.bulk, aes(x = ax1.scores.bulk, y = ax2.scores.bulk, shape=Irrigation),size=5, alpha= 0.8)+
                      scale_color_viridis(discrete = T) +
                      scale_x_continuous(name=paste("PCoA1:\n",round(ax1.bulk,3)*100,"% var. explained", sep=""))+
                      scale_y_continuous(name=paste("PCoA2:\n",round(ax2.bulk,3)*100,"% var. explained", sep=""))+
                      labs(colour = "Treatment",  title = "A. Bulk Soil")+
                      theme(legend.position="none",
                            legend.title = element_text(size=15, face='bold'),
                            plot.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.title = element_text(size = 20, face="bold"),
                            axis.text=element_text(size=16), 
                            axis.title=element_text(size=17,face="bold"),
                            legend.text=element_text(size=15),
                            legend.spacing.x = unit(0.05, 'cm'))+
                      stat_ellipse()
aob.pcoa_bulk.plot

# B. Bray-Curtis - Rhizosphere :
aob.pcoa_rh.plot <- ggplot(data = aob.map.pcoa.rh, aes(x=ax1.scores.rh, y=ax2.scores.rh, colour=Treatment))+
                           theme_bw()+
                    geom_point(data = aob.map.pcoa.rh, aes(x = ax1.scores.rh, y = ax2.scores.rh, shape=Irrigation),size=5, alpha= 0.8)+
                    scale_color_viridis(discrete = T) +
                    scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rh,3)*100,"% var. explained", sep=""))+
                    scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rh,3)*100,"% var. explained", sep=""))+
                    labs(colour = "Treatment",  title = "B. Rhizosphere")+
                    theme(legend.position="right",
                          legend.title = element_text(size=15, face='bold'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 20, face="bold"),
                          axis.text=element_text(size=16), 
                          axis.title=element_text(size=17,face="bold"),
                          legend.text=element_text(size=15),
                          legend.spacing.x = unit(0.05, 'cm'))+
                    stat_ellipse()
aob.pcoa_rh.plot

library(work)

aob.bray.plot <- aob.pcoa_bulk.plot |  aob.pcoa_rh.plot
aob.bray.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.bray.tiff",
       aob.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Jaccard - Bulk Soil :
aob.pcoa_bulk.jac <- ggplot(data = aob.map.pcoa.j.bulk, aes(x=ax1.scores.j.bulk, y=ax2.scores.j.bulk, colour=Treatment))+
                             theme_bw()+
                      geom_point(data = aob.map.pcoa.j.bulk, aes(x = ax1.scores.j.bulk, y = ax2.scores.j.bulk, shape=Irrigation),size=5, alpha= 0.8)+
                      scale_color_viridis(discrete = T) +
                      scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j.bulk,3)*100,"% var. explained", sep=""))+
                      scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j.bulk,3)*100,"% var. explained", sep=""))+
                      labs(colour = "Treatment",  title = "A. Bulk Soil")+
                      theme(legend.position="none",
                            legend.title = element_text(size=15, face='bold'),
                            plot.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.title = element_text(size = 20, face="bold"),
                            axis.text=element_text(size=16), 
                            axis.title=element_text(size=17,face="bold"),
                            legend.text=element_text(size=15),
                            legend.spacing.x = unit(0.05, 'cm'))+
                      stat_ellipse()
aob.pcoa_bulk.jac

# B. Jaccard - Rhizosphere :
aob.pcoa_rh.jac <- ggplot(data = aob.map.pcoa.j.rh, aes(x=ax1.scores.j.rh, y=ax2.scores.j.rh, colour=Treatment))+
                           theme_bw()+
                    geom_point(data = aob.map.pcoa.j.rh, aes(x = ax1.scores.j.rh, y = ax2.scores.j.rh, shape=Irrigation),size=5, alpha= 0.8)+
                    scale_color_viridis(discrete = T) +
                    scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j.rh,3)*100,"% var. explained", sep=""))+
                    scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j.rh,3)*100,"% var. explained", sep=""))+
                    labs(colour = "Treatment",  title = "B. Rhizosphere")+
                    theme(legend.position="right",
                          legend.title = element_text(size=15, face='bold'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 20, face="bold"),
                          axis.text=element_text(size=16), 
                          axis.title=element_text(size=17,face="bold"),
                          legend.text=element_text(size=15),
                          legend.spacing.x = unit(0.05, 'cm'))+
                    stat_ellipse()
aob.pcoa_rh.jac

library(work)
aob.jac.plot <- aob.pcoa_bulk.jac |  aob.pcoa_rh.jac
aob.jac.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.jac.tiff",
       aob.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Weighted UniFrac - Bulk Soil :
aob.pcoa_bulk.wUF <- ggplot(data = aob.map.pcoa.wUF.bulk, aes(x=ax1.scores.wUF.bulk, y=ax2.scores.wUF.bulk, colour=Treatment))+
                             theme_bw()+
                      geom_point(data = aob.map.pcoa.wUF.bulk, aes(x = ax1.scores.wUF.bulk, y = ax2.scores.wUF.bulk, shape=Irrigation),size=5, alpha= 0.8)+
                      scale_color_viridis(discrete = T) +
                      scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.bulk,3)*100,"% var. explained", sep=""))+
                      scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.bulk,3)*100,"% var. explained", sep=""))+
                      labs(colour = "Treatment",  title = "A. Bulk Soil")+
                      theme(legend.position="none",
                            legend.title = element_text(size=15, face='bold'),
                            plot.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.title = element_text(size = 20, face="bold"),
                            axis.text=element_text(size=16), 
                            axis.title=element_text(size=17,face="bold"),
                            legend.text=element_text(size=15),
                            legend.spacing.x = unit(0.05, 'cm'))+
                      stat_ellipse()
aob.pcoa_bulk.wUF

# B. Weighted UniFrac - Rhizosphere :
aob.pcoa_rh.wUF <- ggplot(data = aob.map.pcoa.wUF.rh, aes(x=ax1.scores.wUF.rh, y=ax2.scores.wUF.rh, colour=Treatment))+
                           theme_bw()+
                    geom_point(data = aob.map.pcoa.wUF.rh, aes(x = ax1.scores.wUF.rh, y = ax2.scores.wUF.rh, shape=Irrigation),size=5, alpha= 0.8)+
                    scale_color_viridis(discrete = T) +
                    scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.rh,3)*100,"% var. explained", sep=""))+
                    scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.rh,3)*100,"% var. explained", sep=""))+
                    labs(colour = "Treatment",  title = "B. Rhizosphere")+
                    theme(legend.position="right",
                          legend.title = element_text(size=15, face='bold'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 20, face="bold"),
                          axis.text=element_text(size=16), 
                          axis.title=element_text(size=17,face="bold"),
                          legend.text=element_text(size=15),
                          legend.spacing.x = unit(0.05, 'cm'))+
                    stat_ellipse()
aob.pcoa_rh.wUF

library(work)
aob.wUF.plot <- aob.pcoa_bulk.wUF |  aob.pcoa_rh.wUF
aob.wUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.wUF.tiff",
       aob.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Unweighted UniFrac - Bulk Soil :
aob.pcoa_bulk.uwUF <- ggplot(data = aob.map.pcoa.uwUF.bulk, aes(x=ax1.scores.uwUF.bulk, y=ax2.scores.uwUF.bulk, colour=Treatment))+
                             theme_bw()+
                      geom_point(data = aob.map.pcoa.uwUF.bulk, aes(x = ax1.scores.uwUF.bulk, y = ax2.scores.uwUF.bulk, shape=Irrigation),size=5, alpha= 0.8)+
                      scale_color_viridis(discrete = T) +
                      scale_x_continuous(name=paste("PCoA1:\n",round(ax1.uwUF.bulk,3)*100,"% var. explained", sep=""))+
                      scale_y_continuous(name=paste("PCoA2:\n",round(ax2.uwUF.bulk,3)*100,"% var. explained", sep=""))+
                      labs(colour = "Treatment",  title = "A. Bulk Soil")+
                      theme(legend.position="none",
                            legend.title = element_text(size=15, face='bold'),
                            plot.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.title = element_text(size = 20, face="bold"),
                            axis.text=element_text(size=16), 
                            axis.title=element_text(size=17,face="bold"),
                            legend.text=element_text(size=15),
                            legend.spacing.x = unit(0.05, 'cm'))+
                      stat_ellipse()
aob.pcoa_bulk.uwUF
aob.pcoa_bulk.uwUF.id <- aob.pcoa_bulk.uwUF+geom_text_repel(aes(label = SampleID),size = 3, max.overlaps = Inf)
aob.pcoa_bulk.uwUF.id
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_UnweightedUniFrac_bulk_id.tiff",
       aob.pcoa_bulk.uwUF.id, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# B. Unweighted UniFrac - Rhizosphere :
aob.pcoa_rh.uwUF <- ggplot(data = aob.map.pcoa.uwUF.rh, aes(x=ax1.scores.uwUF.rh, y=ax2.scores.uwUF.rh, colour=Treatment))+
                           theme_bw()+
                    geom_point(data = aob.map.pcoa.uwUF.rh, aes(x = ax1.scores.uwUF.rh, y = ax2.scores.uwUF.rh, shape=Irrigation),size=5, alpha= 0.8)+
                    scale_color_viridis(discrete = T) +
                    scale_x_continuous(name=paste("PCoA1:\n",round(ax1.uwUF.rh,3)*100,"% var. explained", sep=""))+
                    scale_y_continuous(name=paste("PCoA2:\n",round(ax2.uwUF.rh,3)*100,"% var. explained", sep=""))+
                    labs(colour = "Treatment",  title = "B. Rhizosphere")+
                    theme(legend.position="right",
                          legend.title = element_text(size=15, face='bold'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 20, face="bold"),
                          axis.text=element_text(size=16), 
                          axis.title=element_text(size=17,face="bold"),
                          legend.text=element_text(size=15),
                          legend.spacing.x = unit(0.05, 'cm'))+
                    stat_ellipse()
aob.pcoa_rh.uwUF

library(work)
aob.uwUF.plot <- aob.pcoa_bulk.uwUF |  aob.pcoa_rh.uwUF
aob.uwUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.uwUF.tiff",
       aob.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

############################################################################################
# PERMANOVA FOR BULK SOIL AND RHIZOSPHERE
############################################################################################

# A. Bray-Curtis - Bulk Soil : 
set.seed(13)
aob.adonis.bulk <- adonis2(aob.bulk_dist_bc ~ Irrigation*Treatment*Date, data=aob.min.meta.bulk, 
                           permutation=999,
                           method="bray", 
                           strata = NULL) # only treatment is significant
aob.adonis.bulk
set.seed(13)
aob.adonis.bulk.irri <- adonis2(aob.bulk_dist_bc ~ Irrigation, data=aob.min.meta.bulk, 
                           permutation=999,
                           method="bray", 
                           strata = NULL) # not significant
aob.adonis.bulk.irri
set.seed(13)
#aob.adonis.bulk.irri2 <- adonis2(aob.bulk_dist_bc ~ Irrigation*Treatment*Date, data=aob.min.meta.bulk, 
                         #permutation=perm,
                         #method="bray") # not significant
#aob.adonis.bulk.irri2
#perm = how(nperm = 999,
           #within = Within(type="free"), 
           #plots = with(aob.min.meta.bulk, Plots(strata=Treatment, type="free")),
           #blocks = aob.min.meta.bulk$Date)
set.seed(13)
aob.adonis.bulk.trt <- adonis2(aob.bulk_dist_bc ~ Treatment, data=aob.min.meta.bulk, 
                           permutation=999,
                           method="bray", 
                           strata = NULL) # significant (p val = 0.001***)
aob.adonis.bulk.trt
set.seed(13)
aob.adonis.bulk.date <- adonis2(aob.bulk_dist_bc ~ Date, data=aob.min.meta.bulk, 
                           permutation=999,
                           method="bray", 
                           strata = NULL) # not significant
aob.adonis.bulk.date

# B. Bray-Curtis - Rhizosphere : 
set.seed(13)
aob.adonis.rh <- adonis2(aob.rh_dist_bc ~ Irrigation*Treatment*Date, data=aob.min.meta.rh, 
                         permutation=999,
                         method="bray", 
                         strata = NULL) # only treatment is significant
aob.adonis.rh
set.seed(13)
aob.adonis.rh.irri <- adonis2(aob.rh_dist_bc ~ Irrigation, data=aob.min.meta.rh, 
                         permutation=999,
                         method="bray", 
                         strata = NULL) # not significant
aob.adonis.rh.irri
set.seed(13)
aob.adonis.rh.irri2 <- adonis2(aob.rh_dist_bc ~ Irrigation, data=aob.min.meta.rh, 
                           permutation=999,
                           method="bray", 
                           strata = aob.min.meta.rh$Treatment) # not significant
aob.adonis.rh.irri2
set.seed(13)
aob.adonis.rh.trt <- adonis2(aob.rh_dist_bc ~ Treatment, data=aob.min.meta.rh, 
                         permutation=999,
                         method="bray", 
                         strata = NULL) # treatment is significant ( p val = 0.001***)
aob.adonis.rh.trt
set.seed(13)
aob.adonis.rh.date <- adonis2(aob.rh_dist_bc ~ Date, data=aob.min.meta.rh, 
                         permutation=999,
                         method="bray", 
                         strata = NULL) # not significant
aob.adonis.rh.date

# A. Jaccard - Bulk Soil : 
set.seed(13)
aob.adonis.jac.bulk <- adonis2(aob.bulk_dist_jac ~ Irrigation*Treatment*Date, data=aob.min.meta.bulk, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
aob.adonis.jac.bulk
# B. Jaccard - Rhizosphere : 
set.seed(13)
aob.adonis.jac.rh <- adonis2(aob.rh_dist_jac ~ Irrigation*Treatment*Date, data=aob.min.meta.rh, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
aob.adonis.jac.rh

# A. Weighted UniFrac - Bulk Soil : 
set.seed(13)
aob.adonis.wuF.bulk <- adonis2(aob.bulk_dist_wUF ~ Irrigation*Treatment*Date, data=aob.min.meta.bulk, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.wuF.bulk
# B. Weighted UniFrac - Rhizosphere : 
set.seed(13)
aob.adonis.wuF.rh <- adonis2(aob.rh_dist_wUF ~ Irrigation*Treatment*Date, data=aob.min.meta.rh, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.wuF.rh

# A. Unweighted UniFrac - Bulk Soil : 
set.seed(13)
aob.adonis.uwuF.bulk <- adonis2(aob.bulk_dist_uwUF ~ Irrigation*Treatment*Date, data=aob.min.meta.bulk, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.uwuF.bulk
# B. Unweighted UniFrac - Rhizosphere : 
set.seed(13)
aob.adonis.uwuF.rh <- adonis2(aob.rh_dist_uwUF ~ Irrigation*Treatment*Date, data=aob.min.meta.rh, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.uwuF.rh

########################################################################################
# Pairwise ccomparison analyses accross treatments and between irrigation within date
########################################################################################
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Plotting the Distance Matrices

# 1. Weighted UniFrac

# wrangle distance matrix into a longer data frame
aob.bulk_dist_wUF
aob.bulk_dist_wUF.matrix = melt(as.matrix(aob.bulk_dist_wUF))
# remove self-comparisons
aob.bulk_dist_wUF.matrix = aob.bulk_dist_wUF.matrix[aob.bulk_dist_wUF.matrix$X1 != aob.bulk_dist_wUF.matrix$X2,]
# select sample data
tmp_sam_data = tibble("sample"=aob.physeq_bulk1@sam_data$SampleID,
                      "irrigation"=aob.physeq_bulk1@sam_data$Irrigation,
                      "treatment"=aob.physeq_bulk1@sam_data$Treatment,
                      "date"=aob.physeq_bulk1@sam_data$Date,)
# combined distance matrix with sample data
colnames(tmp_sam_data) = c("X1", "irrigation1", "treatment1", "date1")
tmp_data <- left_join(aob.bulk_dist_wUF.matrix, tmp_sam_data, by = "X1")
colnames(tmp_sam_data) = c("X2", "irrigation2", "treatment2", "date2")
tmp_data <- left_join(tmp_data, tmp_sam_data, by = "X2")
# plot
tmp_data$date1 <- factor(tmp_data$date1, levels = unique(tmp_data$date1))
tmp_data$date2 <- factor(tmp_data$date2, levels = unique(tmp_data$date2))
ggplot(tmp_data, aes(x = date2, y = value)) +
    theme_bw() +
    #geom_point(position = position_dodge(width = .75)) +
    #geom_boxplot()+
    geom_boxplot(aes(fill=irrigation2), alpha=0.6) +
    #scale_fill_viridis(discrete = T)+
    scale_fill_manual(values = c("#e5f5e0","#31a354"))+
    facet_wrap(~treatment2, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))

# 2.  Bray - Curtis

aob.bulk_dist_bc
# wrangle distance matrix into a longer dataframe
aob.bulk_dist_bc
aob.bulk_dist_bc.matrix = melt(as.matrix(aob.bulk_dist_bc))
# remove self-comparisons
aob.bulk_dist_bc.matrix = aob.bulk_dist_bc.matrix[aob.bulk_dist_bc.matrix$X1 != aob.bulk_dist_bc.matrix$X2,]
# select sample data
tmp_sam_data = tibble("sample"=aob.physeq_bulk1@sam_data$SampleID,
                      "irrigation"=aob.physeq_bulk1@sam_data$Irrigation,
                      "treatment"=aob.physeq_bulk1@sam_data$Treatment,
                      "date"=aob.physeq_bulk1@sam_data$Date,)
# combined distance matrix with sample data
colnames(tmp_sam_data) = c("X1", "irrigation1", "treatment1", "date1")
tmp_data_bc <- left_join(aob.bulk_dist_bc.matrix, tmp_sam_data, by = "X1")
colnames(tmp_sam_data) = c("X2", "irrigation2", "treatment2", "date2")
tmp_data_bc <- left_join(tmp_data_bc, tmp_sam_data, by = "X2")
tmp_data_bc$date1 <- factor(tmp_data_bc$date1, levels = unique(tmp_data_bc$date1))
tmp_data_bc$date2 <- factor(tmp_data_bc$date2, levels = unique(tmp_data_bc$date2))
# Plot
ggplot(tmp_data_bc, aes(x = date2, y = value)) +
    theme_bw() +
    #geom_point(position = position_dodge(width = .75)) +
    #geom_boxplot()+
    geom_boxplot(aes(fill=irrigation2), alpha=0.6) +
    #scale_fill_viridis(discrete = T)+
    scale_fill_manual(values = c("#e5f5e0","#31a354"))+
    facet_wrap(~treatment2, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))


m04dist_bc
# wrangle distance matrix into a longer dataframe
m04dist_bc.matrix = melt(as.matrix(m04dist_bc))
# remove self-comparisons
m04dist_bc.matrix = m04dist_bc.matrix[m04dist_bc.matrix$X1 != m04dist_bc.matrix$X2,]
# select sample data
m04physeq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
m04physeq1 <- prune_taxa(taxa_sums(m04physeq)>0, m04physeq)
sort(taxa_sums(m04physeq1), decreasing =F)
m04physeq1

tmp_sam_data.m04 = tibble("sample"=m04physeq1@sam_data$SampleID,
                      "irrigation"=m04physeq1@sam_data$Irrigation,
                      "treatment"=m04physeq1@sam_data$Treatment,
                      "date"=m04physeq1@sam_data$Date,)
# combined distance matrix with sample data
colnames(tmp_sam_data.m04) = c("X1", "irrigation1", "treatment1", "date1")
tmp_data_m04 <- left_join(m04dist_bc.matrix, tmp_sam_data.m04, by = "X1")
colnames(tmp_sam_data.m04) = c("X2", "irrigation2", "treatment2", "date2")
tmp_data_m04 <- left_join(tmp_data_m04, tmp_sam_data.m04, by = "X2")

tmp_data_m04$date1 <- factor(tmp_data_m04$date1, levels = unique(tmp_data_m04$date1))
tmp_data_m04$date2 <- factor(tmp_data_m04$date2, levels = unique(tmp_data_m04$date2))
ggplot(tmp_data_m04, aes(x = irrigation2, y = value)) +
    theme_bw() +
    #geom_point(position = position_dodge(width = .75)) +
    #geom_boxplot()+
    geom_boxplot(aes(fill=irrigation2), alpha=0.6) +
    #scale_fill_viridis(discrete = T)+
    scale_fill_manual(values = c("#e5f5e0","#31a354"))+
    #facet_wrap(~treatment2, scales = "free_x") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))


# A. Pairwise Adonis Among Treatment

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
pw.bulk.trt_bc <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_bc, 
                                                  aob.min.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
pw.bulk.trt_bc # all pairwise comparisons are significant (p val =0.001**)

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.trt_wUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_wUF, 
                                                  aob.min.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
pw.bulk.trt_wUF # all pairwise comparisons are significant (p val =0.001**)

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.trt_uwUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_uwUF, 
                                                  aob.min.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
pw.bulk.trt_uwUF # all pairwise comparisons are significant (p val =0.001**)

# B. Pairwise Adonis Among Date

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
pw.bulk.dat_bc <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_bc, 
                                                  aob.min.meta.bulk$Date,
                                                  p.adjust.m = "BH")
pw.bulk.dat_bc # none are significant 

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_wUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_wUF, 
                                                  aob.min.meta.bulk$Date,
                                                  p.adjust.m = "BH")
pw.bulk.dat_wUF # none are significant 

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_uwUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_uwUF, 
                                                  aob.min.meta.bulk$Date,
                                                  p.adjust.m = "BH")
pw.bulk.dat_uwUF # none are significant 

# C. Pairwise Adonis Between Irrigation within Treatment and Date

# separate the aob.asv table by date and treatment:
# 1.
m04.asv <- aob.asv.min.bulk1 %>% select(S1, S2, S9, S10, S15, S16, S23, S24)
m04.asv1 <- m04.asv[rowSums(m04.asv)>0,]
d04.asv <- aob.asv.min.bulk1 %>% select(S3, S4, S11, S12, S13, S14, S21, S22)
d04.asv1 <- d04.asv[rowSums(d04.asv)>0,]
k04.asv <- aob.asv.min.bulk1 %>% select(S5, S6, S7, S8, S17, S18, S19, S20)
k04.asv1 <- k04.asv[rowSums(k04.asv)>0,]

# 2. 
m06.asv <- aob.asv.min.bulk1 %>% select(S25, S26, S33, S34, S39, S40, S47, S48)
m06.asv1 <- m06.asv[rowSums(m06.asv)>0,]
d06.asv <- aob.asv.min.bulk1 %>% select(S27, S28, S35, S36, S37, S38, S45, S46)
d06.asv1 <- d06.asv[rowSums(d06.asv)>0,]
k06.asv <- aob.asv.min.bulk1 %>% select(S29, S30, S31, S32, S41, S42, S43, S44)
k06.asv1 <- k06.asv[rowSums(k06.asv)>0,]

# 3. 
m0705.asv <- aob.asv.min.bulk1 %>% select(S49, S50, S57, S58, S63, S64, S71, S72)
m0705.asv1 <- m0705.asv[rowSums(m0705.asv)>0,]
d0705.asv <- aob.asv.min.bulk1 %>% select(S51, S52, S59, S60, S61, S62, S69, S70)
d0705.asv1 <- d0705.asv[rowSums(d0705.asv)>0,]
k0705.asv <- aob.asv.min.bulk1 %>% select(S53, S54, S55, S56, S65, S66, S67, S68)
k0705.asv1 <- k0705.asv[rowSums(k0705.asv)>0,]

# 4.
m0720.asv <- aob.asv.min.bulk1 %>% select(S73, S74, S81, S82, S87, S88, S95, S96)
m0720.asv1 <- m0720.asv[rowSums(m0720.asv)>0,]
d0720.asv <- aob.asv.min.bulk1 %>% select(S75, S76, S83, S84, S85, S86, S93, S94)
d0720.asv1 <- d0720.asv[rowSums(d0720.asv)>0,]
k0720.asv <- aob.asv.min.bulk1 %>% select(S77, S78, S79, S80, S89, S90, S91, S92)
k0720.asv1 <- k0720.asv[rowSums(k0720.asv)>0,]

# 5.
m09.asv <- aob.asv.min.bulk1 %>% select(S97, S98, S105, S106, S111, S112, S119, S120)
m09.asv1 <- m09.asv[rowSums(m09.asv)>0,]
d09.asv <- aob.asv.min.bulk1 %>% select(S99, S100, S107, S108, S109, S110, S117, S118)
d09.asv1 <- d09.asv[rowSums(d09.asv)>0,]
k09.asv <- aob.asv.min.bulk1 %>% select(S101, S102, S103, S104, S113, S114, S115, S116)
k09.asv1 <- k09.asv[rowSums(k09.asv)>0,]

# separate the metadata by date and treatment:
# 1. 
m04.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "04-28-22" 
                                  & aob.min.meta.bulk$Treatment == "Mineral fertilized"),]
d04.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "04-28-22" 
                                  & aob.min.meta.bulk$Treatment == "Biodynamic"),]
k04.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "04-28-22" 
                                  & aob.min.meta.bulk$Treatment == "Conventional"),]
# 2.  
m06.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "06-01-22" 
                                  & aob.min.meta.bulk$Treatment == "Mineral fertilized"),]
d06.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "06-01-22" 
                                  & aob.min.meta.bulk$Treatment == "Biodynamic"),]
k06.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "06-01-22" 
                                  & aob.min.meta.bulk$Treatment == "Conventional"),]
# 3. 
m0705.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-05-22" 
                                  & aob.min.meta.bulk$Treatment == "Mineral fertilized"),]
d0705.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-05-22" 
                                  & aob.min.meta.bulk$Treatment == "Biodynamic"),]
k0705.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-05-22" 
                                  & aob.min.meta.bulk$Treatment == "Conventional"),]
# 4. 
m0720.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-20-22" 
                                  & aob.min.meta.bulk$Treatment == "Mineral fertilized"),]
d0720.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-20-22" 
                                  & aob.min.meta.bulk$Treatment == "Biodynamic"),]
k0720.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "07-20-22" 
                                  & aob.min.meta.bulk$Treatment == "Conventional"),]
# 5. 
m09.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "09-13-22" 
                                  & aob.min.meta.bulk$Treatment == "Mineral fertilized"),]
d09.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "09-13-22" 
                                  & aob.min.meta.bulk$Treatment == "Biodynamic"),]
k09.meta <- aob.min.meta.bulk[which(aob.min.meta.bulk$Date == "09-13-22" 
                                  & aob.min.meta.bulk$Treatment == "Conventional"),]

# calculate distance matrix - Bray Curtis
# 1.
m04dist_bc <- vegdist(t(m04.asv1), method = "bray")
d04dist_bc <- vegdist(t(d04.asv1), method = "bray")
k04dist_bc <- vegdist(t(k04.asv1), method = "bray")
# 2.
m06dist_bc <- vegdist(t(m06.asv1), method = "bray")
d06dist_bc <- vegdist(t(d06.asv1), method = "bray")
k06dist_bc <- vegdist(t(k06.asv1), method = "bray")
# 3.
m0705dist_bc <- vegdist(t(m0705.asv1), method = "bray")
d0705dist_bc <- vegdist(t(d0705.asv1), method = "bray")
k0705dist_bc <- vegdist(t(k0705.asv1), method = "bray")
# 4.
m0720dist_bc <- vegdist(t(m0720.asv1), method = "bray")
d0720dist_bc <- vegdist(t(d0720.asv1), method = "bray")
k0720dist_bc <- vegdist(t(k0720.asv1), method = "bray")
# 5.
m09dist_bc <- vegdist(t(m09.asv1), method = "bray")
d09dist_bc <- vegdist(t(d09.asv1), method = "bray")
k09dist_bc <- vegdist(t(k09.asv1), method = "bray")

# pairwise adonis between irrigation within treatment and date
# 1.
set.seed(13)
m04pwc.bc <- pairwiseAdonis::pairwise.adonis(m04dist_bc, m04.meta$Irrigation, p.adjust.m = "BH")
m04pwc.bc # not significant
set.seed(13)
d04pwc.bc <- pairwiseAdonis::pairwise.adonis(d04dist_bc, d04.meta$Irrigation, p.adjust.m = "BH")
d04pwc.bc # not significant
set.seed(13)
k04pwc.bc <- pairwiseAdonis::pairwise.adonis(k04dist_bc, k04.meta$Irrigation, p.adjust.m = "BH")
k04pwc.bc # not significant

# 2.
set.seed(13)
m06pwc.bc <- pairwiseAdonis::pairwise.adonis(m06dist_bc, m06.meta$Irrigation, p.adjust.m = "BH")
m06pwc.bc # not significant
set.seed(13)
d06pwc.bc <- pairwiseAdonis::pairwise.adonis(d06dist_bc, d06.meta$Irrigation, p.adjust.m = "BH")
d06pwc.bc # not significant
set.seed(13)
k06pwc.bc <- pairwiseAdonis::pairwise.adonis(k06dist_bc, k06.meta$Irrigation, p.adjust.m = "BH")
k06pwc.bc # not significant

# 3.
set.seed(13)
m0705pwc.bc <- pairwiseAdonis::pairwise.adonis(m0705dist_bc, m0705.meta$Irrigation, p.adjust.m = "BH")
m0705pwc.bc # not significant
set.seed(13)
d0705pwc.bc <- pairwiseAdonis::pairwise.adonis(d0705dist_bc, d0705.meta$Irrigation, p.adjust.m = "BH")
d0705pwc.bc # not significant
set.seed(13)
k0705pwc.bc <- pairwiseAdonis::pairwise.adonis(k0705dist_bc, k0705.meta$Irrigation, p.adjust.m = "BH")
k0705pwc.bc # not significant

# 4.
set.seed(13)
m0720pwc.bc <- pairwiseAdonis::pairwise.adonis(m0720dist_bc, m0720.meta$Irrigation, p.adjust.m = "BH")
m0720pwc.bc # not significant
set.seed(13)
d0720pwc.bc <- pairwiseAdonis::pairwise.adonis(d0720dist_bc, d0720.meta$Irrigation, p.adjust.m = "BH")
d0720pwc.bc # not significant
set.seed(13)
k0720pwc.bc <- pairwiseAdonis::pairwise.adonis(k0720dist_bc, k0720.meta$Irrigation, p.adjust.m = "BH")
k0720pwc.bc # not significant

# 5.
set.seed(13)
m09pwc.bc <- pairwiseAdonis::pairwise.adonis(m09dist_bc, m09.meta$Irrigation, p.adjust.m = "BH")
m09pwc.bc # not significant
set.seed(13)
d09pwc.bc <- pairwiseAdonis::pairwise.adonis(d09dist_bc, d09.meta$Irrigation, p.adjust.m = "BH")
d09pwc.bc # not significant
set.seed(13)
k09pwc.bc <- pairwiseAdonis::pairwise.adonis(k09dist_bc, k09.meta$Irrigation, p.adjust.m = "BH")
k09pwc.bc # not significant
  











#########################################################################################
# AOB Community Composition
########################################################################################
# Phyloseq object of rarefied data and unrarefied data:
# 1. rarefied data
aob.rare.min.physeq
tax_table(aob.rare.min.physeq)
# merge taxa by species
aob.sp <- tax_glom(aob.rare.min.physeq, taxrank = "Species", NArm = F)
aob.sp.ra <- transform_sample_counts(aob.sp, function(x) x/sum(x))
sample_data(aob.sp.ra)

aob.sp.df <- psmelt(aob.sp.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

colours <- ColourPalleteMulti(aob.sp.df, "Genus", "Species")
colours
install.packages("Polychrome")
library(Polychrome)
# build-in color palette
mycol = glasbey.colors(21)
mycol=c("#FE8F42","#FFFFFF","#0000FF","#FF0000","#009FFF","#00FF00","#000033",
        
        "#FF00B6","#005300","#9A4D42","#00FFBE","#783FC1","#F1085C","#DD00FF",
        
        "#201A01","#1F9698","#FFD300","#FFACFD","#B1CC71", "#720055", "#766C95")

sp_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC", "#737373", "#117744", "#88CCAA",
          "#777711", "#fc8d59", "#fed976", 
          "#771122", "#AA4455", "#DD7788","#774411", "#AA7744", "#DDAA77")


str(aob.sp.df)

install.packages("ggh4x")
library(ggh4x)

aob.sp.df$Date <- as.Date(aob.sp.df$Date , "%m/%d/%Y")
aob.sp.df$Type <- factor(aob.sp.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere")
                  )
aob.sp.df$Treatment <- factor(aob.sp.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral"))
legend <- "AOB Taxa"
library(scales)
#x_cols <- rep(hue_pal()(length(unique(interaction(aob.sp.df$Date, aob.sp.df$Irrigation)))))
aob.sp.plot <- ggplot(aob.sp.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(legend, values=sp_col)+
                     #scale_fill_manual(values=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0"))+
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     facet_nested(~Type+Treatment, nest_line = element_line(linetype = 1), scales="free")+
                     theme(legend.direction = "vertical",legend.position="right") + 
                     guides(fill=guide_legend(ncol=1))+
                     labs(y= "Mean Relative Abundance")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=13, face = "bold"),
                           axis.line.x = element_blank(),
                           axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
                           #axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 12),
                           legend.title = element_text(size=13, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.background = element_blank(),
                           strip.text.x = element_text(size = 14, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           #guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           scale_y_continuous(expand = c(0,0))+
                           guides(x="axis_nested")
                           
aob.sp.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_meanRA_barplot2.eps",
       aob.sp.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)

# 2. unrarefied data
aob.physeq
# merge taxa by species
aob.sp.unrare <- tax_glom(aob.physeq, taxrank = "Species", NArm = F)
aob.sp.unrare.ra <- transform_sample_counts(aob.sp.unrare, function(x) x/sum(x))
aob.sp.unrare.ra

aob.sp.unrare.df <- psmelt(aob.sp.unrare.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

str(aob.sp.unrare.df)
sp_col_x=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC", "#737373", "#117744", "#88CCAA",
          "#777711", "#fc8d59", "#fed976", 
          "#771122", "#AA4455", "#DD7788","#774411", "#AA7744", "#DDAA77","black")
install.packages("ggh4x")
library(ggh4x)

aob.sp.unrare.df$Date <- as.Date(aob.sp.unrare.df$Date , "%m/%d/%Y")
aob.sp.unrare.df$Type <- factor(aob.sp.unrare.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere")
                  )
aob.sp.unrare.df$Treatment <- factor(aob.sp.unrare.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral"))
legend <- "AOB Taxa"
library(scales)
#x_cols <- rep(hue_pal()(length(unique(interaction(aob.sp.df$Date, aob.sp.df$Irrigation)))))
aob.sp.unrare.plot <- ggplot(aob.sp.unrare.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(legend, values=sp_col_x)+
                     #scale_fill_manual(values=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0"))+
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     facet_nested(~Type+Treatment, nest_line = element_line(linetype = 1), scales="free")+
                     theme(legend.direction = "vertical",legend.position="right") + 
                     guides(fill=guide_legend(ncol=1))+
                     labs(y= "Mean Relative Abundance")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=13, face = "bold"),
                           axis.line.x = element_blank(),
                           axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
                           #axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 12),
                           legend.title = element_text(size=13, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.background = element_blank(),
                           strip.text.x = element_text(size = 14, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           #guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           scale_y_continuous(expand = c(0,0))+
                           guides(x="axis_nested")
                           
aob.sp.unrare.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_meanRA_unrare_barplot.eps",
       aob.sp.unrare.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)


###################################################################################################
# Compare the analyses with rarefaction to 1282 reads (remove sample with reads < 1000)
###################################################################################################

# ASV Table
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F)
set.seed(13)
aob.rare.1282.seq <- rarefy_even_depth(aob.physeq, sample.size = 1282,
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.1282.seq # 1 samples removed (S11), 116 ASVs were removed
# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
aob.rare.1282 <- ggrare(aob.rare.1282.seq, step = 1, color = "Type", label = "SampleID", se = FALSE)
#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare.1282 <- aob.rare.1282 + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
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
 
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve_1282.jpg",
       plot.aob.rare.1282, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)
aob.rare.1282.seq


# rarefy to minimum sequencing depth
set.seed(13)
aob.rare.min.physeq <- rarefy_even_depth(aob.physeq, sample.size = min(sample_sums(aob.physeq)),
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min.physeq), decreasing = F) # 172 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed

# run the ggrare function attached in the file "generating_rarecurve.r"
aob.rare.min <- ggrare(aob.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare.min <- aob.rare.min + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
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
 
plot.aob.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve_min.jpg",
       plot.aob.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)


















