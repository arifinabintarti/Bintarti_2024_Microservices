#############################################################################################
# Analysis of amoA gene of AOB Illumina MiSeq Data - rarefied 1282
#############################################################################################

# Date : 13 July 2023
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

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(aob.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")

## make a phyloseq object of the asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aob.asv))
aob.asv.ord  <- aob.asv[ ,re_order]
aob.asv.physeq = phyloseq::otu_table(aob.asv.ord, taxa_are_rows = TRUE) # asv table
phyloseq::sample_names(aob.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
phyloseq::sample_names(aob.asv.physeq) <- paste0("S", phyloseq::sample_names(aob.asv.physeq))


# phyloseq object of the taxonomy table
aob.tax <- column_to_rownames(aob.tax, var = "ASVid")
aob.tax.physeq = phyloseq::tax_table(as.matrix(aob.tax)) # taxonomy table

# phyloseq object of the metadata
str(meta_micro)
meta_micro$Date <- factor(meta_micro$Date, levels = c("4/28/22", "6/1/22", "7/5/22", "7/20/22", "9/13/22"),
                          labels = c("Apr-28", "Jun-01", "Jul-05", "Jul-20", "Sep-13"))
rownames(meta_micro) <- phyloseq::sample_names(aob.asv.physeq)
aob.meta.physeq <- phyloseq::sample_data(meta_micro)# meta data
phyloseq::sample_names(aob.meta.physeq)

# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/AOB-rooted-tree/')
AOB_rooted_tree <- ape::read.tree("tree.nwk")
AOB_rooted_tree
# tree visualization
library("treeio")
library("ggtree")
p <- ggtree(AOB_rooted_tree, layout = "circular") + 
  geom_tiplab(size=3, color="purple")
p 


# make phyloseq object
library(phyloseq)
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq,AOB_rooted_tree)
aob.physeq
aob.asv.ord <- as.data.frame(otu_table(aob.physeq))
aob.asv.ord
#write.csv(aob.asv.ord, file = "aob.asv.ord.csv")
sample_data(aob.physeq)$SampleID <- paste0("S", sample_data(aob.physeq)$SampleID)
sample_data(aob.physeq)

# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
aob.rare <- ggrare(aob.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
#install.packages("colorBlindness")
library(colorBlindness)
displayAvailablePalette(color="white")
#Palette <- c("#1F968BFF","#FDE725FF")
PairedColor12Steps
Brown2Blue10Steps
Palette <- c("#FF7F00", "#662F00")
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

setwd('D:/Fina/INRAE_Project/microservices_fig/AOB/')
ggsave("AOB_rarecurve.tiff",
       plot.aob.rare, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

###################################################################################################
# rarefaction to 1282 reads (remove sample with reads < 1000)
###################################################################################################

# ASV Table
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F)
set.seed(333)
aob.rare.1282.seq <- rarefy_even_depth(aob.physeq, sample.size = 1282,
                                       rngseed = 333, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.1282.seq # 1 samples removed (S11), 116 ASVs were removed
 sort(rowSums(otu_table(aob.rare.1282.seq), na.rm = FALSE, dims = 1), decreasing = F)
######################################################################################################

# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 
#### AOB alpha diversity using rarefied data to the minimum sequencing depth (reads=858) ###
colSums(otu_table(aob.rare.1282.seq))
aob.asv.rare1k <- as.data.frame(otu_table(aob.rare.1282.seq))
dim(aob.asv.rare1k) # 1222 ASVs
aob.asv.rare1k_pa <- 1*(aob.asv.rare1k>0)
aob.s <- specnumber(aob.asv.rare1k, MARGIN = 2) # richness
aob.richness <- as.data.frame(aob.s) 
aob.h <- diversity(t(aob.asv.rare1k), index = 'shannon') # Shannon index
aob.shannon <- as.data.frame(aob.h)
aob.d <- diversity(t(aob.asv.rare1k), index = 'simpson') # Simpson index
aob.simpson <- as.data.frame(aob.d)
aob.inv.d <- diversity(t(aob.asv.rare1k), index = 'invsimpson')

# 1. Richness

# Line plot of AOB richness 
aob.meta.df <- data.frame(meta_micro)
aob.meta.df.sub <- aob.meta.df %>% filter(SampleID != 11)# filter out S11 from the metadata
aob.meta.df.sub$Richness <- aob.s
aob.meta.df.sub$Shannon <- aob.h
aob.meta.df.sub$Simpson <- aob.d
aob.meta.df.sub$InvSimpson <- aob.inv.d
#aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
str(aob.meta.df.sub)
aob.meta.df.sub$Type <- factor(aob.meta.df.sub$Type, levels = c("BS", "RS"),
                               labels = c("Bulk Soil", "Rhizosphere"))
aob.meta.df.sub$Treatment <- factor(aob.meta.df.sub$Treatment, levels = c("D", "K", "M"),
                                    labels = c("BIODYN", "CONFYM", "CONMIN"))
aob.meta.df.sub$SampleID<-factor(aob.meta.df.sub$SampleID)
aob.meta.df.sub$PlotID<-factor(aob.meta.df.sub$PlotID)
aob.meta.df.sub$Irrigation<-factor(aob.meta.df.sub$Irrigation)
aob.meta.df.sub$Block<-factor(aob.meta.df.sub$Block)
aob.meta.df.sub$x<-factor(aob.meta.df.sub$x)
aob.meta.df.sub$DxD<-factor(aob.meta.df.sub$DxD)
aob.meta.df.sub$period<-factor(aob.meta.df.sub$period)
aob.meta.df.sub$rep<-factor(aob.meta.df.sub$rep)
aob.meta.df.sub$rep2<-factor(aob.meta.df.sub$rep2)
aob.meta.df.sub$var2<-factor(aob.meta.df.sub$var2)
aob.meta.df.sub$var3<-factor(aob.meta.df.sub$var3)
aob.meta.df.sub[sapply(aob.meta.df.sub, is.character)] <- 
 lapply(aob.meta.df.sub[sapply(aob.meta.df.sub, is.character)], as.numeric)
aob.meta.df.sub[sapply(aob.meta.df.sub, is.integer)] <- 
 lapply(aob.meta.df.sub[sapply(aob.meta.df.sub, is.integer)], as.numeric)

# tidy up the data frame
aob.meta.df.tidy <- aob.meta.df.sub %>%
  group_by(Irrigation, Treatment, Date,  Type, var2,var3) %>%
  summarize(Mean.Rich=mean(Richness),
            Mean.Sha=mean(Shannon),
            Mean.Simp=mean(Simpson),
            Mean.invsimp=mean(InvSimpson))
str(aob.meta.df.tidy)

# plotting line chart
#remotes::install_github("Nowosad/rcartocolor")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
carto_pal(n = NULL, 'Vivid')
display_carto_pal(7, "Vivid")
display_carto_pal(12, "Safe")
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("viridis")
library(viridis)
aob.meta.df.sub$Date <- factor(aob.meta.df.sub$Date, levels = unique(aob.meta.df.sub$Date))
aob.rich.plot <- ggplot(aob.meta.df.tidy, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
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
aob.rich.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich3.eps",
       aob.min.rich.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_rich.tiff",
       aob.rich.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)


dev.off()

# 2. Shannon

#Line plot of AOB Shannon
library(scales)
aob.sha.plot <- ggplot(aob.meta.df.tidy, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
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

aob.sha.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha2.eps",
       aob.min.sha.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_sha.tiff",
       aob.sha.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#. 3. Inverse Simpson

#Line plot of AOB Inverse Simpson
aob.invsimp.plot <- ggplot(aob.meta.df.tidy, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
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

aob.invsimp.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp.eps",
       aob.min.invsimp.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_invsimp.tiff",
       aob.invsimp.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#############################################################################

# 1. Richness

# Box plot of AOB Richness

#install.packages('remotes')
#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("rstatix")
#install.packages("sf")
library(remotes)
library(rstatix)
library(sf)
library(ggpattern)

aob.meta.df.sub$x
aob.meta.df.sub.ed <- aob.meta.df.sub %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

# Richness: plotting the pairwise comparisons across treatment 
aob.rich.pwc.plot <- ggplot(aob.meta.df.sub, aes(x=Irrigation, y=Richness)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw() +
  labs(y="AOB Richness")+
  expand_limits(y = 0)+
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
aob.rich.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
xy.rich.bulk <- emm.rich.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
xy.rich.rh <- emm.rich.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
df.xy.rich.bulk <- as.data.frame(xy.rich.bulk)
df.xy.rich.rh <- as.data.frame(xy.rich.rh)
df.xy.rich.all <- rbind(df.xy.rich.bulk, df.xy.rich.rh) 
df.xy.rich.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aob.rich.pwc.plot2 <- aob.rich.pwc.plot + 
  stat_pvalue_manual(df.xy.rich.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.rich.pwc.plot2

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich_mean_boxplot.eps",
       aob.min.rich.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_rich_mean_boxplot.tiff",
       aob.rich.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# richness between irrigations
aob.meta.df.sub.ed2 <- aob.meta.df.sub.ed
aob.meta.df.sub.ed2$Date <- factor(aob.meta.df.sub.ed2$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
aob.meta.df.sub.ed2$Date <- as.Date(aob.meta.df.sub.ed2$Date)

ann_text<-data.frame(Date = as.Date("2022-07-30"), Richness = 0, Type="Bulk Soil",label="Rewetting")
stat_text.BS <- data.frame(Date = as.Date("2022-08-25"), Richness = 140, Type="Bulk Soil",Treatment="CONMIN", label="C = **")
stat_text.RS <- data.frame(Date = as.Date("2022-08-25"), Richness = 140, Type="Rhizosphere",Treatment="CONMIN", label="C = **")
# plot BS and RS together
aob.rich.pwc.irri.plot <- ggplot(aob.meta.df.sub.ed2, aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  labs(y="AOB Richness")+
  facet_grid(Type~ Treatment,scales="free_x")+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_vline(data=filter(aob.meta.df.sub.ed2, Type=="Bulk Soil"),aes(xintercept = as.Date("2022-07-14")), linetype="dashed", colour="darkgrey") +
 geom_text(data = ann_text,label=ann_text$label, colour="grey25")+
 geom_text(data = stat_text.BS,label=stat_text$label, colour="black", size=4, fontface="bold")+
 geom_text(data = stat_text.RS,label=stat_text$label, colour="black", size=4, fontface="bold")
#annotate(geom = "text", x = as.Date("2022-07-10"), y = -5, hjust = 0, size = 4, label = "Rewetting", color = "black")
aob.rich.pwc.irri.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_rich_irri_boxplot.tiff",
       aob.rich.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
#setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
# 2. Shannon

# Shannon: plotting the significance across treatment

aob.sha.pwc.plot <- ggplot(aob.meta.df.sub, aes(x=Irrigation, y=Shannon)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="AOB Shannon")+
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
aob.sha.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
xy.sha.bulk <- emm.sha.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
xy.sha.rh <- emm.sha.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
df.xy.sha.bulk <- as.data.frame(xy.sha.bulk)
df.xy.sha.rh <- as.data.frame(xy.sha.rh)
df.xy.sha.all <- rbind(df.xy.sha.bulk, df.xy.sha.rh) 
df.xy.sha.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aob.sha.pwc.plot2 <- aob.sha.pwc.plot + 
  stat_pvalue_manual(df.xy.sha.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.sha.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha_boxplot.eps",
       aob.min.sha.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_sha_boxplot.tiff",
       aob.sha.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# shannon between irrigations

aob.sha.pwc.irri.plot <- ggplot(aob.meta.df.sub.ed, aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  labs(y="AOB Shannon")+
  facet_grid(Type~ Treatment,scales="free_x")+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
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
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_sha_irri_boxplot.tiff",
       aob.sha.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
#setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')

#############################################################################################################################
# Separate between Bulk soil and Rhizosphere Alpha Diversity

#_Bulk Soil Richness_________________________________________________________________________________________________________________
# plot just BS
aob.meta.bulk <- aob.meta.df.sub[1:119,]
str(aob.meta.bulk)
aob.meta.bulk$SampleID<-factor(aob.meta.bulk$SampleID)
aob.meta.bulk$PlotID<-factor(aob.meta.bulk$PlotID)
aob.meta.bulk$Irrigation<-factor(aob.meta.bulk$Irrigation)
aob.meta.bulk$Block<-factor(aob.meta.bulk$Block)
aob.meta.bulk$x<-factor(aob.meta.bulk$x)
aob.meta.bulk$rep<-factor(aob.meta.bulk$rep)
aob.meta.bulk$rep2<-factor(aob.meta.bulk$rep2)
aob.meta.bulk$var2<-factor(aob.meta.bulk$var2)
aob.meta.bulk$DxD<-factor(aob.meta.bulk$DxD)
aob.meta.bulk$period<-factor(aob.meta.bulk$period)
aob.meta.bulk$period2<-factor(aob.meta.bulk$period2)

aob.meta.bulk.edit <- aob.meta.bulk
str(aob.meta.bulk.edit)
#aob.meta.bulk.edit$Date <- factor(aob.meta.bulk.edit$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          #labels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
#aob.meta.bulk.edit$Date <- factor(aob.meta.bulk.edit$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          #labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
#aob.meta.bulk.edit$Date <- as.Date(aob.meta.bulk.edit$Date)

#aob.meta.bulk.edit$Date <- factor(aob.meta.bulk.edit$Date)
aob.meta.bulk.edit <- aob.meta.bulk.edit %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

#stat_text.BS2 <- data.frame(Date = as.Date("2022-05-01"), Richness = 100,Treatment="BIODYN", label="C **")

stat_text.BS2 <- data.frame(Date = 0.5, Richness = 9,Treatment="BIODYN", label="C **")
AOB.BS.Richness.plot <- ggplot(aob.meta.bulk.edit , aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('BIODYN.control', 'BIODYN.drought', 'CONFYM.control', 
                             'CONFYM.drought', 'CONMIN.control', 'CONMIN.drought'))+
  labs(y="AOB Richness", title="Bulk Soil\nA", subtitle = "C**")+# tag="Bulk Soil")+
  facet_wrap(~ Treatment)+
  #scale_y_continuous(limits = c(0, NA))+
  scale_y_continuous(limits = c(0, 120))+
  theme(legend.position = "none",
        legend.title = element_text(size=26, face='bold'),
        legend.text = element_text(size=26),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank(),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                        size = 25,
                        lineheight = 1,
                        padding = margin(5.5, 5.5, 5.5, 5.5),
                        margin = margin(0, 0, 5.5, 0),
                        linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") 
 #geom_vline(data=filter(aob.meta.df.sub.ed2, Type=="Bulk Soil"),aes(xintercept = as.Date("2022-07-14")), linetype="dashed", colour="darkgrey") +
 #geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", colour="darkgrey") 
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 #geom_label(data = stat_text.BS2,label=stat_text.BS2$label, hjust=0, colour="black", size=6, fontface="bold")
 
AOB.BS.Richness.plot


#_Rizhosphere Richness________________________________________________________________________________________________________________________
# plot just RS
aob.meta.rh <- aob.meta.df.sub[120:191,]
str(aob.meta.rh)
aob.meta.rh$SampleID<-factor(aob.meta.rh$SampleID)
aob.meta.rh$PlotID<-factor(aob.meta.rh$PlotID)
aob.meta.rh$Irrigation<-factor(aob.meta.rh$Irrigation)
aob.meta.rh$Block<-factor(aob.meta.rh$Block)
aob.meta.rh$x<-factor(aob.meta.rh$x)
aob.meta.rh$rep<-factor(aob.meta.rh$rep)
aob.meta.rh$rep2<-factor(aob.meta.rh$rep2)
aob.meta.rh$var2<-factor(aob.meta.rh$var2)
aob.meta.rh$var3<-factor(aob.meta.rh$var3)
aob.meta.rh$DxD<-factor(aob.meta.rh$DxD)
aob.meta.rh$period<-factor(aob.meta.rh$period)
#aob.meta.rh$period2<-factor(aob.meta.rh$period2)

aob.meta.rh.edit <- aob.meta.rh
aob.meta.rh.edit$Date <- factor(aob.meta.rh.edit$Date, levels = c("Apr-28", "Jun-01", "Jul-05"))
aob.meta.rh.edit <- aob.meta.rh.edit %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

aob.stat_text.RS.rich <- data.frame(Date = 0.5, Richness = 9,Treatment="BIODYN", label="C **")

AOB.RS.Richness.plot <- ggplot(aob.meta.rh.edit , aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('BIODYN.control', 'BIODYN.drought', 'CONFYM.control', 
                             'CONFYM.drought', 'CONMIN.control', 'CONMIN.drought'))+
  labs(y="AOB Richness", title="Rhizosphere\nG", subtitle = "C**")+
  facet_wrap(~ Treatment)+
  #ylim(0,139)+
  scale_y_continuous(limits = c(0, 140))+
  theme(legend.position = "none",
        legend.title = element_text(size=22, face='bold'),
        legend.text = element_text(size=26),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank(),
        strip.text = element_text(size=25),
        axis.text.y = element_text(size = 23),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                        size = 25,
                        lineheight = 1,
                        padding = margin(5.5, 5.5, 5.5, 5.5),
                        margin = margin(0, 0, 5.5, 0),
                        linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 #geom_label(data = aob.stat_text.RS.rich,label=aob.stat_text.RS.rich$label, hjust=0,colour="black", size=6, fontface="bold")

 AOB.RS.Richness.plot


#_Bulk Soil Shannon Index________________________________________________________________________________________________________________________________

#aob.stat_text.BS.shan <- data.frame(Date = as.Date("2022-04-28"), Shannon = 2.3,Treatment="BIODYN", label="C **")
aob.stat_text.BS.shan <- data.frame(Date = 0.5, Shannon = 0.35,Treatment="BIODYN", label="C **")
AOB.BS.Shannon.plot <- ggplot(aob.meta.bulk.edit , aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('BIODYN.control', 'BIODYN.drought', 'CONFYM.control', 
                             'CONFYM.drought', 'CONMIN.control', 'CONMIN.drought'))+
  labs(y="AOB Shannon", title="D", subtitle = "C**")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_text(size=22, face='bold'),
        legend.text = element_text(size=26),
        #axis.text.x  = element_blank(),
        #axis.ticks.x  = element_blank(),
        #strip.text = element_text(size=25),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                        size = 25,
                        lineheight = 1,
                        padding = margin(5.5, 5.5, 5.5, 5.5),
                        margin = margin(0, 0, 5.5, 0),
                        linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))+
 #geom_vline(data=filter(aob.meta.df.sub.ed2, Type=="Bulk Soil"),aes(xintercept = as.Date("2022-07-14")), linetype="dashed", colour="darkgrey") +
 #geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", colour="darkgrey") +
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 annotate(geom = "text", x = 3.6, y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 #geom_label(data = aob.stat_text.BS.shan,label=aob.stat_text.BS.shan$label, hjust=0,colour="black", size=6, fontface="bold")
 
AOB.BS.Shannon.plot

#_Rhizosphere Shannon Index_____________________________________________________________________________________________________________________________________

aob.stat_text.RS.shan <- data.frame(Date = 0.5, Shannon = 0.35,Treatment="BIODYN", label="C *")

AOB.RS.Shannon.plot <- ggplot(aob.meta.rh.edit , aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('BIODYN.control', 'BIODYN.drought', 'CONFYM.control', 
                             'CONFYM.drought', 'CONMIN.control', 'CONMIN.drought'))+
  labs(y="AOB Shannon", title="J", subtitle = "C*")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_text(size=22, face='bold'),
        legend.text = element_text(size=26),
        #axis.text.x  = element_blank(),
        #axis.ticks.x  = element_blank(),
        #strip.text = element_text(size=25),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23,angle = 45, hjust = 1),
        axis.title.y = element_text(size=24),
        axis.title.x =element_blank(),
        plot.title = element_text(size = 27, face = "bold"),
        plot.subtitle = element_textbox_simple(face = "italic",
                        size = 25,
                        lineheight = 1,
                        padding = margin(5.5, 5.5, 5.5, 5.5),
                        margin = margin(0, 0, 5.5, 0),
                        linetype = 1),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0,'lines'))
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 #geom_label(data = aob.stat_text.RS.shan,label=aob.stat_text.RS.shan$label, hjust=0,colour="black", size=6, fontface="bold")

 AOB.RS.Shannon.plot

# SAVE ALL GROUPS SHANNON PLOTS 
library(patchwork)

all.alpha <- ((AOB.BS.Richness.plot/AOB.BS.Shannon.plot/AOB.RS.Richness.plot/AOB.RS.Shannon.plot)|(AOA.BS.Richness.plot/AOA.BS.Shannon.plot/AOA.RS.Richness.plot/AOA.RS.Shannon.plot)|(COM.BS.Richness.plot/COM.BS.Shannon.plot/COM.RS.Richness.plot/COM.RS.Shannon.plot))+
 plot_layout(guides = "collect") & theme(legend.position = 'bottom',legend.title = element_blank())
all.alpha             

setwd('/Users/arifinabintarti/Documents/France/Figures')               
ggsave("Fig.1dpi600.tiff",
       t.all, device = "tiff",
       width = 23.6, height =24,
       #width = 19.8, height =13.5, 
       units= "in", dpi = 600, compression="lzw")




t1=((AOB.BS.Richness.plot|AOA.BS.Richness.plot|COM.BS.Richness.plot)/
    (AOB.BS.Shannon.plot|AOA.BS.Shannon.plot|COM.BS.Shannon.plot))
t2=((AOB.RS.Richness.plot|AOA.RS.Richness.plot|COM.RS.Richness.plot)/
   (AOB.RS.Shannon.plot|AOA.RS.Shannon.plot|COM.RS.Shannon.plot))
t.all=t1/t2

t1=(AOB.BS.Richness.plot/AOB.BS.Shannon.plot/AOB.RS.Richness.plot/AOB.RS.Shannon.plot)
t2=(AOA.BS.Richness.plot/AOA.BS.Shannon.plot/AOA.RS.Richness.plot/AOA.RS.Shannon.plot)
t3=(COM.BS.Richness.plot/COM.BS.Shannon.plot/COM.RS.Richness.plot/COM.RS.Shannon.plot)
t=(t1|t2|t3)+plot_layout(guides = "collect") & theme(legend.position = 'bottom',legend.title = element_blank())

ggsave("t.tiff",
       t, device = "tiff",
       #width = 28.5, height =24,
       width = 20, height =21, 
       units= "in", dpi = 300, compression="lzw")



#((AOB.BS.Richness.plot/AOB.BS.Shannon.plot/AOB.RS.Richness.plot/AOB.RS.Shannon.plot)|(AOA.BS.Richness.plot/AOA.BS.Shannon.plot/AOA.RS.Richness.plot/AOA.RS.Shannon.plot))

#_____________________________________________________________________________________________________________________________________
# Inverse Simpson

# Inverse Simpson: plotting the significance across treatment
aob.invsimp.pwc.plot <- ggplot(aob.meta.df.sub, aes(x=Irrigation, y=InvSimpson)) +
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
aob.invsimp.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
xy.invsimp.bulk <- emm.invsimp.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
xy.invsimp.rh <- emm.invsimp.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
df.xy.invsimp.bulk <- as.data.frame(xy.invsimp.bulk)
df.xy.invsimp.rh <- as.data.frame(xy.invsimp.rh)
df.xy.invsimp.all <- rbind(df.xy.invsimp.bulk, df.xy.invsimp.rh) 
df.xy.invsimp.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aob.invsimp.pwc.plot2 <- aob.invsimp.pwc.plot + 
  stat_pvalue_manual(df.xy.invsimp.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.invsimp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp_all.eps",
       aob.min.invsimp.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_invsimp_all.tiff",
       aob.invsimp.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

###################################################################################
# Beta Diversity Analyses on Rarefied Data: AOB
###################################################################################

# FOR ALL SAMPLES
# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
aob.asv.rare1k <- as.data.frame(otu_table(aob.rare.1282.seq))
#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/AOB.ASV-analysis')
#write.csv(aob.asv.rare1k, file = "aob.rare1k.otutab.csv")
# Bray-Curtis using rarefied data:
aob_dist_bc <- vegdist(t(aob.asv.rare1k), method = "bray")
# Jaccard using rarefied data:
aob.dist_jac <- vegdist(t(aob.asv.rare1k), binary = TRUE, method = "jaccard")
# Weighted UniFrac using rarefied data:
aob.wUF_dist <- UniFrac(aob.rare.1282.seq, weighted=TRUE, normalized = TRUE)
aob.wUF_dist
# Unweighted UniFrac using rarefied data:
aob.uwUF_dist <- UniFrac(aob.rare.1282.seq, weighted=FALSE, normalized = TRUE)
aob.uwUF_dist

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis using rarefied data:
aob_pcoa.bc <- cmdscale(aob_dist_bc, eig=T)
# jaccard using rarefied data:
aob_pcoa.jac <- cmdscale(aob.dist_jac, eig=T)
# Weighted UniFrac using rarefied data:
aob_pcoa.wUF <- cmdscale(aob.wUF_dist, eig=T)
# Unweighted UniFrac using rarefied data:
aob_pcoa.uwUF <- cmdscale(aob.uwUF_dist, eig=T)

# 3. scores of PC1 and PC2

# bray-curtis:
ax1.scores <- aob_pcoa.bc$points[,1]
ax2.scores <- aob_pcoa.bc$points[,2] 
# jaccard:
ax1.scores.j <- aob_pcoa.jac$points[,1]
ax2.scores.j <- aob_pcoa.jac$points[,2]
# Weighted UniFrac using rarefied data:
ax1.scores.wUF <- aob_pcoa.wUF$points[,1]
ax2.scores.wUF <- aob_pcoa.wUF$points[,2]
# Unweighted UniFrac using rarefied data:
ax1.scores.uwUF <- aob_pcoa.uwUF$points[,1]
ax2.scores.uwUF <- aob_pcoa.uwUF$points[,2]

# 4. calculate percent variance explained, then add to plot

#Bray-curtis:
ax1 <- aob_pcoa.bc$eig[1]/sum(aob_pcoa.bc$eig)
ax2 <- aob_pcoa.bc$eig[2]/sum(aob_pcoa.bc$eig)
aob.map.pcoa <- cbind(aob.meta.df.sub,ax1.scores,ax2.scores)
# jaccard
ax1.j <- aob_pcoa.jac$eig[1]/sum(aob_pcoa.jac$eig)
ax2.j <- aob_pcoa.jac$eig[2]/sum(aob_pcoa.jac$eig)
aob.map.pcoa.j <- cbind(aob.meta.df.sub,ax1.scores.j,ax2.scores.j)
# Weighted UniFrac using rarefied:
ax1.wUF <- aob_pcoa.wUF$eig[1]/sum(aob_pcoa.wUF$eig)
ax2.wUF <- aob_pcoa.wUF$eig[2]/sum(aob_pcoa.wUF$eig)
aob.map.pcoa.wUF <- cbind(aob.meta.df.sub,ax1.scores.wUF,ax2.scores.wUF)
# unweighted UniFrac
ax1.uwUF <- aob_pcoa.uwUF$eig[1]/sum(aob_pcoa.uwUF$eig)
ax2.uwUF <- aob_pcoa.uwUF$eig[2]/sum(aob_pcoa.uwUF$eig)
aob.map.pcoa.uwUF <- cbind(aob.meta.df.sub,ax1.scores.uwUF,ax2.scores.uwUF)

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
library(viridis)

set.seed(13)
aob.pcoa_plot <- ggplot(data = aob.map.pcoa, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
  geom_point(data = aob.map.pcoa, aes(x = ax1.scores, y = ax2.scores, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
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
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.jac
jaccard

# c. Weighted UniFrac using rarefied data:

aob.pcoa_plot.wUF <- ggplot(data = aob.map.pcoa.wUF, aes(x=ax1.scores.wUF, y=ax2.scores.wUF, color=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.wUF, aes(x = ax1.scores.wUF, y = ax2.scores.wUF , shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PC1:\n",round(ax1.wUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PC2:\n",round(ax2.wUF,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.wUF

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_WeightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.wUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# d. UnWeighted UniFrac using rarefied data:

aob.pcoa_plot.uwUF <- ggplot(data = aob.map.pcoa.uwUF, aes(x=ax1.scores.uwUF, y=ax2.scores.uwUF, color=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.uwUF, aes(x = ax1.scores.uwUF, y = ax2.scores.uwUF, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PC1:\n",round(ax1.uwUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PC2:\n",round(ax2.uwUF,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB Unweighted UniFrac (rarefied)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.uwUF
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
aob.asv.rare1k <- as.data.frame(otu_table(aob.rare.1282.seq))
aob.asv.bulk <- aob.asv.rare1k[,1:119]
aob.asv.bulk1 <- aob.asv.bulk[rowSums(aob.asv.bulk)>0,]
sort(rowSums(aob.asv.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out/')
#write.csv(aob.asv.bulk1, file = "aob.bulk1.otutab.csv")
aob.bulk_dist_bc <- vegdist(t(aob.asv.bulk1), method = "bray")
# jaccard - Bulk Soil :
aob.bulk_dist_jac <- vegdist(t(aob.asv.bulk1), binary = TRUE, method = "jaccard")
# Weighted UniFrac (rarefied) - Bulk Soil:
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS")
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
ax1.scores.aob.BS <- aob.bulk_pcoa_bc$points[,1]
ax2.scores.aob.BS <- aob.bulk_pcoa_bc$points[,2] 
# Jaccard - Bulk Soil:
ax1.scores.j.bulk <- aob.bulk_pcoa_jac$points[,1]
ax2.scores.j.bulk <- aob.bulk_pcoa_jac$points[,2]
# Weighted UniFrac - Bulk Soil:
ax1.scores.wUF.bulk <- aob.bulk_pcoa_wUF$points[,1]
ax2.scores.wUF.bulk <- aob.bulk_pcoa_wUF$points[,2]
# Unweighted UniFrac - Bulk Soil:
ax1.scores.uwUF.bulk <- aob.bulk_pcoa.uwUF$points[,1]
ax2.scores.uwUF.bulk <- aob.bulk_pcoa.uwUF$points[,2]

# 4. Envfit
env.aob.bulk <- aob.meta.df.sub[1:119,c(14:20, 23,27:29)]
str(env.aob.bulk)
env.aob.bulk <- env.aob.bulk %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg'), as.numeric)
# bray-curtis
set.seed(149)
env_fit.aob.bc.bulk <- envfit(aob.bulk_pcoa_bc, env.aob.bulk, na.rm=TRUE)
env_fit.aob.bc.bulk
# Jaccard 
env_fit.aob.jac <- envfit(aob.bulk_pcoa_jac, env.aob.bulk, na.rm=TRUE)
# Weighted UniFrac
env_fit.aob.wuF <- envfit(aob.bulk_pcoa_wUF, env.aob.bulk, na.rm=TRUE)
# UnWeighted UniFrac
env_fit.aob.uwuF <- envfit(aob.bulk_pcoa.uwUF, env.aob.bulk, na.rm=TRUE)


# 5. calculate percent variance explained, then add to plot
aob.meta.bulk <- aob.meta.df.sub[1:119,]
# Bray-curtis - Bulk Soil:
ax1.aob.BS <- aob.bulk_pcoa_bc$eig[1]/sum(aob.bulk_pcoa_bc$eig)
ax2.aob.BS <- aob.bulk_pcoa_bc$eig[2]/sum(aob.bulk_pcoa_bc$eig)
aob.map.pcoa.bulk <- cbind(aob.meta.bulk,ax1.scores.aob.BS,ax2.scores.aob.BS)
# Jaccard - Bulk Soil:
ax1.j.bulk <- aob.bulk_pcoa_jac$eig[1]/sum(aob.bulk_pcoa_jac$eig)
ax2.j.bulk <- aob.bulk_pcoa_jac$eig[2]/sum(aob.bulk_pcoa_jac$eig)
aob.map.pcoa.j.bulk <- cbind(aob.meta.bulk,ax1.scores.j.bulk,ax2.scores.j.bulk)
# Weighted UniFrac - Bulk Soil:
ax1.wUF.bulk <- aob.bulk_pcoa_wUF$eig[1]/sum(aob.bulk_pcoa_wUF$eig)
ax2.wUF.bulk <- aob.bulk_pcoa_wUF$eig[2]/sum(aob.bulk_pcoa_wUF$eig)
aob.map.pcoa.wUF.bulk <- cbind(aob.meta.bulk,ax1.scores.wUF.bulk,ax2.scores.wUF.bulk)
# Unweighted UniFrac - Bulk Soil:
ax1.uwUF.bulk <- aob.bulk_pcoa.uwUF$eig[1]/sum(aob.bulk_pcoa.uwUF$eig)
ax2.uwUF.bulk <- aob.bulk_pcoa.uwUF$eig[2]/sum(aob.bulk_pcoa.uwUF$eig)
aob.map.pcoa.uwUF.bulk <- cbind(aob.meta.bulk,ax1.scores.uwUF.bulk,ax2.scores.uwUF.bulk)

#################################################################################################

# RHIZOSPHERE

# Bray-Curtis - Rhizosphere :
aob.asv.rh <- aob.asv.rare1k[,120:191]
aob.asv.rh1 <- aob.asv.rh[rowSums(aob.asv.rh)>0,]
sort(rowSums(aob.asv.rh1, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(aob.asv.rh1) #831
aob.rh_dist_bc <- vegdist(t(aob.asv.rh1), method = "bray")
# jaccard - Rhizosphere :
aob.rh_dist_jac <- vegdist(t(aob.asv.rh1), binary = TRUE, method = "jaccard") 
# Weighted UniFrac (rarefied) - Rhizosphere :
aob.physeq_rh <- subset_samples(aob.rare.1282.seq, Type=="RS")
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
ax1.scores.aob.RS <- aob.rh_pcoa_bc$points[,1]
ax2.scores.aob.RS <- aob.rh_pcoa_bc$points[,2] 
# Jaccard - Rhizosphere :
ax1.scores.j.rh <- aob.rh_pcoa_jac$points[,1]
ax2.scores.j.rh <- aob.rh_pcoa_jac$points[,2]
# Weighted UniFrac - Rhizosphere :
ax1.scores.wUF.rh <- aob.rh_pcoa_wUF$points[,1]
ax2.scores.wUF.rh <- aob.rh_pcoa_wUF$points[,2]
# Unweighted UniFrac - Rhizosphere :
ax1.scores.uwUF.rh <- aob.rh_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rh <- aob.rh_pcoa.uwUF$points[,2]

## 4. Envfit
env.aob.rh <- aob.meta.df.sub[120:191,c(31:40)]
str(env.aob.rh)
env.aob.rh <- env.aob.rh %>% mutate_at(colnames(env.aob.rh), as.numeric)
# bray-curtis
env_fit.aob.bc.rh <- envfit(aob.rh_pcoa_bc, env.aob.rh, na.rm=TRUE)
env_fit.aob.bc.rh
# Jaccard 
env_fit.aob.jac.rh <- envfit(aob.rh_pcoa_jac, env.aob.rh, na.rm=TRUE)
# Weighted UniFrac
env_fit.aob.wuF.rh <- envfit(aob.rh_pcoa_wUF, env.aob.rh, na.rm=TRUE)
# UnWeighted UniFrac
env_fit.aob.uwuF.rh <- envfit(aob.rh_pcoa.uwUF, env.aob.rh, na.rm=TRUE)

# 5. calculate percent variance explained, then add to plot
aob.meta.rh <- aob.meta.df.sub[120:191,]
# Bray-curtis - Rhizosphere :
ax1.aob.RS <- aob.rh_pcoa_bc$eig[1]/sum(aob.rh_pcoa_bc$eig)
ax2.aob.RS <- aob.rh_pcoa_bc$eig[2]/sum(aob.rh_pcoa_bc$eig)
aob.map.pcoa.rh <- cbind(aob.meta.rh,ax1.scores.aob.RS,ax2.scores.aob.RS)
# Jaccard - Rhizosphere :
ax1.j.rh <- aob.rh_pcoa_jac$eig[1]/sum(aob.rh_pcoa_jac$eig)
ax2.j.rh <- aob.rh_pcoa_jac$eig[2]/sum(aob.rh_pcoa_jac$eig)
aob.map.pcoa.j.rh <- cbind(aob.meta.rh,ax1.scores.j.rh,ax2.scores.j.rh)
# Weighted UniFrac - Rhizosphere :
ax1.wUF.rh <- aob.rh_pcoa_wUF$eig[1]/sum(aob.rh_pcoa_wUF$eig)
ax2.wUF.rh <- aob.rh_pcoa_wUF$eig[2]/sum(aob.rh_pcoa_wUF$eig)
aob.map.pcoa.wUF.rh <- cbind(aob.meta.rh,ax1.scores.wUF.rh,ax2.scores.wUF.rh)
# Unweighted UniFrac - Rhizosphere :
ax1.uwUF.rh <- aob.rh_pcoa.uwUF$eig[1]/sum(aob.rh_pcoa.uwUF$eig)
ax2.uwUF.rh <- aob.rh_pcoa.uwUF$eig[2]/sum(aob.rh_pcoa.uwUF$eig)
aob.map.pcoa.uwUF.rh <- cbind(aob.meta.rh,ax1.scores.uwUF.rh,ax2.scores.uwUF.rh)

###############################################################################
# 6. PCoA Plot 

#require("ggrepel")
library(ggrepel)
#install.packages("viridis")
library(viridis)

# A. Bray-Curtis - Bulk Soil :
# a. Bray-Curtis:
set.seed(33)
A.bc <- as.list(env_fit.aob.bc.bulk$vectors) #shortcutting ef$vectors
pvals.bc<-as.data.frame(A.bc$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.bc <- as.data.frame(scores(env_fit.aob.bc.bulk, display="vectors"))
env.scores2.bc <- cbind(env.scores1.bc, pvals.bc)
env.scores3.bc <- cbind(env.scores2.bc,Variable=rownames(env.scores2.bc))
env.scores4.bc <- subset(env.scores3.bc,pvals.bc<0.05)
set.seed(33)
mult <-.53
paletteMartin

AOB.PCoA.BS <- ggplot(data = aob.map.pcoa.bulk, aes(x=ax1.scores.aob.BS, y=ax2.scores.aob.BS, colour=Treatment))+
  theme_bw()+
  geom_point(aes(color = aob.map.pcoa.bulk$Treatment, shape = aob.map.pcoa.bulk$Irrigation), size = 4,stroke=2) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#E69F00","#E69F00","#009E73","#009E73","#FF618C","#FF618C","#FF618C","#FF618C",
                            "#E69F00","#E69F00","#009E73","#009E73","#009E73","#009E73","#E69F00","#E69F00",
                            "#FF618C","#FF618C","#FF618C","#FF618C","#009E73","#009E73","#E69F00","#E69F00")) +
 #scale_fill_manual(values = c("#E69F00","#009E73","#FF618C","#FF618C",
                            # "#E69F00","#009E73","#009E73","#E69F00",
                            # "#FF618C","#FF618C","#009E73","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.aob.BS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.aob.BS,3)*100,"% var. explained", sep=""))+
  labs(title = "Bulk Soil", subtitle = "A. AOB")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  geom_mark_ellipse(aes(fill = aob.map.pcoa.bulk$PlotID, label = aob.map.pcoa.bulk$Block),label.fontsize = 20, 
                    expand = 0, linewidth = NA, show.legend = FALSE)
AOB.PCoA.BS
#setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB.PCoA.BS.tiff",
       AOB.PCoA.BS, device = "tiff",
       width = 15, height = 15, 
       units= "in", dpi = 600)

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("aob.bray.plotid_block.tiff",
     aob.pcoa_bulk.plot, device = "tiff",
       width = 6, height =5, 
       units= "in", dpi = 600)

AOB.PCoA.BS.fert <- ggplot(data = aob.map.pcoa.bulk, aes(x=ax1.scores.aob.BS, y=ax2.scores.aob.BS, colour=Treatment))+
  theme_classic()+
  geom_point(aes(color = aob.map.pcoa.bulk$Treatment, shape = aob.map.pcoa.bulk$Irrigation), size = 4,stroke=1) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.aob.BS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.aob.BS,3)*100,"% var. explained", sep=""))+
  labs(title="Bulk Soil", subtitle = "A. AOB")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face="bold"),
        plot.subtitle = element_text(size = 20, face="bold"),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=19))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  stat_ellipse(aes(colour = aob.map.pcoa.bulk$Treatment))
AOB.PCoA.BS.fert

# B. Bray-Curtis - Rhizosphere :
AOB.PCoA.RS <- ggplot(data = aob.map.pcoa.rh, aes(x=ax1.scores.aob.RS, y=ax2.scores.aob.RS, colour=Treatment))+
  theme_bw()+
  geom_point(aes(color = aob.map.pcoa.rh$Treatment, shape = aob.map.pcoa.rh$Irrigation), size = 4, stroke=2) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#E69F00","#E69F00","#009E73","#009E73","#FF618C","#FF618C","#FF618C","#FF618C",
                            "#E69F00","#E69F00","#009E73","#009E73","#009E73","#009E73","#E69F00","#E69F00",
                            "#FF618C","#FF618C","#FF618C","#FF618C","#009E73","#009E73","#E69F00","#E69F00")) +
 #scale_fill_manual(values = c("#E69F00","#009E73","#FF618C","#FF618C",
                            # "#E69F00","#009E73","#009E73","#E69F00",
                            # "#FF618C","#FF618C","#009E73","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.aob.RS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.aob.RS,3)*100,"% var. explained", sep=""))+
  labs(title = "Rhizosphere", subtitle="B. AOB")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  geom_mark_ellipse(aes(fill = aob.map.pcoa.rh$PlotID, label = aob.map.pcoa.rh$Block),label.fontsize = 20, 
                    expand = 0, linewidth = NA, show.legend = FALSE)
AOB.PCoA.RS

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("aob.bray.plotid.rh_block.tiff",
     aob.pcoa_rh.plot, device = "tiff",
       width = 6, height =5, 
       units= "in", dpi = 600)

install.packages("patchwork")
library(patchwork)

aob.bray.plot.envfit <- aob.pcoa_bulk.plot|aob.pcoa_rh.plot
  
aob.bray.plot.envfit
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.bray.tiff",
       aob.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)


AOB.PCoA.RS.fert <- ggplot(data = aob.map.pcoa.rh, aes(x=ax1.scores.aob.RS, y=ax2.scores.aob.RS, colour=Treatment))+
  theme_classic()+
  geom_point(aes(color = aob.map.pcoa.rh$Treatment, shape = aob.map.pcoa.rh$Irrigation), size = 4,stroke=1) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.aob.RS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.aob.RS,3)*100,"% var. explained", sep=""))+
  labs(title="Rhizosphere", subtitle = "B. AOB")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face="bold"),
        plot.subtitle = element_text(size = 20, face="bold"),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=19))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  stat_ellipse(aes(colour = aob.map.pcoa.rh$Treatment))
                    #expand = 0, linewidth = NA, show.legend = FALSE)
AOB.PCoA.RS.fert

#setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
##############################################################################################################################

# Combines all PCoA figures
PCoA.All <- ((AOB.PCoA.BS / AOA.PCoA.BS / COM.PCoA.BS) | (AOB.PCoA.RS / AOA.PCoA.RS / COM.PCoA.RS)) +
 plot_layout(guides = "collect") & theme(legend.position = 'bottom',legend.title = element_blank(),
                                         legend.text = element_text(size=25))
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Supp.Fig.3x.tiff",
       PCoA.All, device = "tiff",
       width = 27, height = 30, 
       units= "in", dpi = 600,compression = "lzw")
ggsave("Supp.Fig.3dpi300.tiff",
       PCoA.All, device = "tiff",
       width = 27, height = 30, 
       units= "in", dpi =300,compression = "lzw")

# Combines all PCoA figures
PCoA.Fert.All <- ((AOB.PCoA.BS.fert / AOA.PCoA.BS.fert / COM.PCoA.BS.fert) | (AOB.PCoA.RS.fert / AOA.PCoA.RS.fert / COM.PCoA.RS.fert)) +
 plot_layout(guides = "collect") & theme(legend.position = 'bottom',legend.title = element_blank(),
                                         legend.text = element_text(size=17))
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Supp.Fig.6dpi300.tiff",
       PCoA.Fert.All, device = "tiff",
       width = 12, height = 14, 
       units= "in", dpi =300,compression = "lzw")


###############################################################################################################################
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
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("aob.jac.tiff",
       aob.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Weighted UniFrac - Bulk Soil :
set.seed(33)
A.wu <- as.list(env_fit.aob.wuF$vectors) #shortcutting ef$vectors
pvals.wu<-as.data.frame(A.wu$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.wu <- as.data.frame(scores(env_fit.aob.wuF, display="vectors"))
env.scores2.wu <- cbind(env.scores1.wu, pvals)
env.scores3.wu <- cbind(env.scores2.wu,Variable=rownames(env.scores2.wu))
env.scores4.wu <- subset(env.scores3.wu,pvals<0.05)
set.seed(33)
mult <-.25
aob.pcoa_bulk.wUF <- ggplot(data = aob.map.pcoa.wUF.bulk, aes(x=ax1.scores.wUF.bulk, y=ax2.scores.wUF.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.wUF.bulk, aes(x = ax1.scores.wUF.bulk, y = ax2.scores.wUF.bulk, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.bulk,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  geom_segment(data=env.scores4.wu,
               aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), 
               arrow = arrow(length = unit(0.3, "cm")),
               colour = "grey",inherit.aes = FALSE)+
  geom_text_repel(data = env.scores4.wu,
                  aes(x = mult*Dim1, y = mult*Dim2, label = Variable),
                  size = 5,fontface="bold",
                  position=position_jitter(width=0.03,height=0.001), inherit.aes = FALSE)+
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
  geom_point(data = aob.map.pcoa.wUF.rh, aes(x = ax1.scores.wUF.rh, y = ax2.scores.wUF.rh, shape=Irrigation),size=5, alpha= 0.6)+
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

aob.wUF.plot.envfit <- aob.pcoa_bulk.wUF |  aob.pcoa_rh.wUF
aob.wUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.wUF.tiff",
       aob.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("aob.wUF.tiff",
       aob.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
ggsave("aob.wUF.envfit.tiff",
       aob.wUF.plot.envfit, device = "tiff",
       width = 16, height = 6, 
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
#aob.pcoa_bulk.uwUF.id <- aob.pcoa_bulk.uwUF+geom_text_repel(aes(label = SampleID),size = 3, max.overlaps = Inf)
#aob.pcoa_bulk.uwUF.id
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

aob.uwUF.plot <- aob.pcoa_bulk.uwUF |  aob.pcoa_rh.uwUF
aob.uwUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("aob.uwUF.tiff",
       aob.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("aob.uwUF.tiff",
       aob.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
############################################################################################
# PERMANOVA FOR BULK SOIL AND RHIZOSPHERE
############################################################################################
# A. Bray-Curtis - Bulk Soil : 

block.aob=as.factor(aob.meta.bulk$Block)
plot.aob=as.factor(aob.meta.bulk$PlotID)
TxI.aob=as.factor(aob.meta.bulk$x)
trt.aob=as.factor(aob.meta.bulk$Treatment)
irri.aob=as.factor(aob.meta.bulk$Irrigation)

## Betadisper for treatment
aob.trt.mod <- betadisper(aob.bulk_dist_bc, trt.aob)
aob.trt.mod
boxplot(aob.trt.mod)
# Null hypothesis of no difference in dispersion between groups
set.seed(13)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.aob.bs <- permutest(aob.trt.mod, permutations = 999, pairwise = T)
permod.aob.bs # there is significant differences in dispersion between groups
# the variances among groups are not homogeneous,
hsd.aob.bs <- TukeyHSD(aob.trt.mod) #which groups differ in relation to their variances
hsd.aob.bs

## Betadisper for irrigation
aob.irri.mod <- betadisper(aob.bulk_dist_bc, irri.aob)
aob.irri.mod
boxplot(aob.irri.mod)
# Null hypothesis of no difference in dispersion between groups
set.seed(13)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.aob.bs.irri <- permutest(aob.irri.mod, permutations = 999, pairwise = T)
permod.aob.bs.irri # there are no significant differences in dispersion between groups
# the variances among groups are homogeneous,
hsd.aob.bs.irri <- TukeyHSD(aob.irri.mod) #which groups differ in relation to their variances
hsd.aob.bs.irri

# 1. Using adonis2 package with defined perm to restrict the permutation 
set.seed(13)
aob.adonis.bulk.bc <- adonis2(aob.bulk_dist_bc ~ Irrigation*Treatment+Date, strata=block.aob,
                              data=aob.meta.bulk, 
                              permutations = 999) # significant
aob.adonis.bulk.bc

# similar with below:
perm1.aob = how(within = Within(type="free"), 
            plots = Plots(type = "none"),
            blocks = block.aob,
            nperm = 999,
            observed = TRUE)
set.seed(13)
aob.adonis.bulk.bc.perm1 <- adonis2(aob.bulk_dist_bc ~ Irrigation, data=aob.meta.bulk, 
                                    permutations = perm1.aob)
aob.adonis.bulk.bc.perm1

# another way to use how()

# Since our intent is to focus on the variation among treatments, 
# we need to restrict the permutations so that plots are permuted within each block, but plots are not permuted across blocks.
# these two ways are equivalent:
CTRL.t1.aob <- how(within = Within(type = "free"),
               plots = Plots(type = "none"),
               blocks = block.aob,
               nperm = 999,
               observed = TRUE)
# and
CTRL.t2.aob <- how(within = Within(type = "free"),
               plots = Plots(strata = block.aob, type = "none"),
               nperm = 999,
               observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(13)
aob.adonis.bulk.bc.CTRL.t2 <- adonis2(aob.bulk_dist_bc ~ Irrigation, data=aob.meta.bulk, 
                                      permutations = CTRL.t2.aob)
aob.adonis.bulk.bc.CTRL.t2
set.seed(13)
aob.adonis.bulk.wUF.CTRL.t2 <- adonis2(aob.bulk_dist_wUF ~ Irrigation, data=aob.meta.bulk, 
                                      permutations = CTRL.t2.aob)
aob.adonis.bulk.wUF.CTRL.t2
# 2. Using ANOSIM package and define the strata
set.seed(13)
aob.bc.anosim <- anosim(aob.bulk_dist_bc,
                        grouping = aob.meta.bulk$x, permutations = 999, strata = block.aob)
summary(aob.bc.anosim) # NOT SIGNIFICANT


# test the permanova for farming system
set.seed(13)
aob.adonis.bulk <- adonis2(aob.bulk_dist_bc ~ Treatment*Irrigation*Date, data=aob.meta.bulk, 
                           permutation=999) # only treatment is significant
aob.adonis.bulk
####################################################################################################

# B. Bray-Curtis - Rhizosphere : 
block.aob.rh=as.factor(aob.meta.rh$Block)
plot.aob.rh=as.factor(aob.meta.rh$PlotID)
TxI.aob.rh=as.factor(aob.meta.rh$x)
trt.aob.rh=as.factor(aob.meta.rh$Treatment)
irri.aob.rh=as.factor(aob.meta.rh$Irrigation)

# and
CTRL.t2.rh.aob <- how(within = Within(type = "free"),
                  plots = Plots(strata = block.aob.rh, type = "none"),
                  nperm = 999,
                  observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(13)
aob.adonis.rh.bc.CTRL.t2 <- adonis2(aob.rh_dist_bc ~ block.aob.rh+Irrigation, data=aob.meta.rh, 
                                    permutations = CTRL.t2.rh.aob)
aob.adonis.rh.bc.CTRL.t2
set.seed(13)
aob.adonis.rh.wUF.CTRL.t2 <- adonis2(aob.rh_dist_wUF ~ Irrigation, data=aob.meta.rh, 
                                    permutations = CTRL.t2.rh.aob)
aob.adonis.rh.wUF.CTRL.t2
# 2. Using ANOSIM package and define the strata
set.seed(13)
aob.bc.anosim.rh <- anosim(aob.rh_dist_wUF,
                           grouping = irri.aob.rh, permutations = 999, strata = block.aob.rh)
summary(aob.bc.anosim.rh) # SIGNIFICANT

set.seed(13)
aob.adonis.rh <- adonis2(aob.rh_dist_bc ~ Treatment*Irrigation*Date, data=aob.meta.rh,
                         permutation=999) # only treatment is significant
aob.adonis.rh






00set.seed(13)
aob.adonis.rh.irri <- adonis2(aob.rh_dist_bc ~ Irrigation, data=aob.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
aob.adonis.rh.irri

set.seed(13)
aob.adonis.rh.irri2 <- adonis2(aob.rh_dist_bc ~ Irrigation, data=aob.meta.rh, 
                               permutation=999,
                               method="bray", 
                               strata = aob.meta.rh$Treatment) # not significant
aob.adonis.rh.irri2

set.seed(13)
aob.adonis.rh.trt <- adonis2(aob.rh_dist_bc ~ Treatment, data=aob.meta.rh, 
                             permutation=999,
                             method="bray", 
                             strata = NULL) # treatment is significant ( p val = 0.001***)
aob.adonis.rh.trt

set.seed(13)
aob.adonis.rh.date <- adonis2(aob.rh_dist_bc ~ Date, data=aob.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
aob.adonis.rh.date

# A. Jaccard - Bulk Soil : 
set.seed(13)
aob.adonis.jac.bulk <- adonis2(aob.bulk_dist_jac ~ Irrigation*Treatment*Date, data=aob.meta.bulk, 
                               permutation=999,
                               method="jaccard", 
                               strata = NULL)
aob.adonis.jac.bulk
# B. Jaccard - Rhizosphere : 
set.seed(13)
aob.adonis.jac.rh <- adonis2(aob.rh_dist_jac ~ Irrigation*Treatment*Date, data=aob.meta.rh, 
                             permutation=999,
                             method="jaccard", 
                             strata = NULL)
aob.adonis.jac.rh

# A. Weighted UniFrac - Bulk Soil : 
set.seed(13)
aob.adonis.wuF.bulk <- adonis2(aob.bulk_dist_wUF ~ Irrigation*Treatment*Date, data=aob.meta.bulk, 
                               permutation=999, 
                               strata = NULL)
aob.adonis.wuF.bulk
# B. Weighted UniFrac - Rhizosphere : 
set.seed(13)
aob.adonis.wuF.rh <- adonis2(aob.rh_dist_wUF ~ Irrigation*Treatment*Date, data=aob.meta.rh, 
                             permutation=999, 
                             strata = NULL)
aob.adonis.wuF.rh

# A. Unweighted UniFrac - Bulk Soil : 
set.seed(13)
aob.adonis.uwuF.bulk <- adonis2(aob.bulk_dist_uwUF ~ Irrigation*Treatment*Date, data=aob.meta.bulk, 
                                permutation=999, 
                                strata = NULL)
aob.adonis.uwuF.bulk
# B. Unweighted UniFrac - Rhizosphere : 
set.seed(13)
aob.adonis.uwuF.rh <- adonis2(aob.rh_dist_uwUF ~ Irrigation*Treatment*Date, data=aob.meta.rh, 
                              permutation=999, 
                              strata = NULL)
aob.adonis.uwuF.rh

########################################################################################
# Pairwise comparison analyses across treatments and between irrigation within date
########################################################################################
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
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


# A. Pairwise Adonis Among Treatment (overall)

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
pw.bulk.trt_bc <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_bc, 
                                                  aob.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
pw.bulk.trt_bc # all pairwise comparisons are significant (p val =0.001**)

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.trt_wUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_wUF, 
                                                   aob.meta.bulk$Treatment,
                                                   p.adjust.m = "BH")
pw.bulk.trt_wUF # all pairwise comparisons are significant (p val =0.001**)

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.trt_uwUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_uwUF, 
                                                    aob.meta.bulk$Treatment,
                                                    p.adjust.m = "BH")
pw.bulk.trt_uwUF # all pairwise comparisons are significant (p val =0.001**)

# B. Pairwise Adonis Among Date

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
pw.bulk.dat_bc <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_bc, 
                                                  aob.meta.bulk$Date,
                                                  p.adjust.m = "BH")
pw.bulk.dat_bc # none are significant 

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_wUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_wUF, 
                                                   aob.meta.bulk$Date,
                                                   p.adjust.m = "BH")
pw.bulk.dat_wUF # none are significant 

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_uwUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_uwUF, 
                                                    aob.meta.bulk$Date,
                                                    p.adjust.m = "BH")
pw.bulk.dat_uwUF # none are significant 

# C. Pairwise Adonis Between Irrigation within Treatment and Date

# separate the aob.asv table by date and treatment:

# 1.04-28-2022
# mineral:
m04.asv <- aob.asv.bulk1 %>% select(S1, S2, S9, S10, S15, S16, S23, S24)
m04.asv1 <- m04.asv[rowSums(m04.asv)>0,]
# biodynamic:
d04.asv <- aob.asv.bulk1 %>% select(S3, S4, S12, S13, S14, S21, S22)
d04.asv1 <- d04.asv[rowSums(d04.asv)>0,]
# conventional:
k04.asv <- aob.asv.bulk1 %>% select(S5, S6, S7, S8, S17, S18, S19, S20)
k04.asv1 <- k04.asv[rowSums(k04.asv)>0,]

# 2. 06-01-2022
# mineral:
m06.asv <- aob.asv.bulk1 %>% select(S25, S26, S33, S34, S39, S40, S47, S48)
m06.asv1 <- m06.asv[rowSums(m06.asv)>0,]
# biodynamic:
d06.asv <- aob.asv.bulk1 %>% select(S27, S28, S35, S36, S37, S38, S45, S46)
d06.asv1 <- d06.asv[rowSums(d06.asv)>0,]
# conventional:
k06.asv <- aob.asv.bulk1 %>% select(S29, S30, S31, S32, S41, S42, S43, S44)
k06.asv1 <- k06.asv[rowSums(k06.asv)>0,]

# 3. 07-05-2022
# mineral:
m0705.asv <- aob.asv.bulk1 %>% select(S49, S50, S57, S58, S63, S64, S71, S72)
m0705.asv1 <- m0705.asv[rowSums(m0705.asv)>0,]
# biodynamic:
d0705.asv <- aob.asv.bulk1 %>% select(S51, S52, S59, S60, S61, S62, S69, S70)
d0705.asv1 <- d0705.asv[rowSums(d0705.asv)>0,]
# conventional:
k0705.asv <- aob.asv.bulk1 %>% select(S53, S54, S55, S56, S65, S66, S67, S68)
k0705.asv1 <- k0705.asv[rowSums(k0705.asv)>0,]

# 4. 07-20-2022
# mineral:
m0720.asv <- aob.asv.bulk1 %>% select(S73, S74, S81, S82, S87, S88, S95, S96)
m0720.asv1 <- m0720.asv[rowSums(m0720.asv)>0,]
# biodynamic:
d0720.asv <- aob.asv.bulk1 %>% select(S75, S76, S83, S84, S85, S86, S93, S94)
d0720.asv1 <- d0720.asv[rowSums(d0720.asv)>0,]
# conventional:
k0720.asv <- aob.asv.bulk1 %>% select(S77, S78, S79, S80, S89, S90, S91, S92)
k0720.asv1 <- k0720.asv[rowSums(k0720.asv)>0,]

# 5. 09-13-2022
# mineral:
m09.asv <- aob.asv.bulk1 %>% select(S97, S98, S105, S106, S111, S112, S119, S120)
m09.asv1 <- m09.asv[rowSums(m09.asv)>0,]
# biodynamic:
d09.asv <- aob.asv.bulk1 %>% select(S99, S100, S107, S108, S109, S110, S117, S118)
d09.asv1 <- d09.asv[rowSums(d09.asv)>0,]
# conventional:
k09.asv <- aob.asv.bulk1 %>% select(S101, S102, S103, S104, S113, S114, S115, S116)
k09.asv1 <- k09.asv[rowSums(k09.asv)>0,]

# separate the metadata by date and treatment:
# 1. 
m04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                    & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                    & aob.meta.bulk$Treatment == "Biodynamic"),]
k04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                    & aob.meta.bulk$Treatment == "Conventional"),]
# 2.  
m06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                    & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                    & aob.meta.bulk$Treatment == "Biodynamic"),]
k06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                    & aob.meta.bulk$Treatment == "Conventional"),]
# 3. 
m0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                      & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                      & aob.meta.bulk$Treatment == "Biodynamic"),]
k0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                      & aob.meta.bulk$Treatment == "Conventional"),]
# 4. 
m0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                      & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                      & aob.meta.bulk$Treatment == "Biodynamic"),]
k0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                      & aob.meta.bulk$Treatment == "Conventional"),]
# 5. 
m09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                    & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                    & aob.meta.bulk$Treatment == "Biodynamic"),]
k09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                    & aob.meta.bulk$Treatment == "Conventional"),]

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

# calculate distance matrix - Weighted UniFrac

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

#########################################################################################
# AOB Community Composition
########################################################################################
# Phyloseq object of rarefied data and unrarefied data:
# 1. rarefied data
aob.rare.1282.seq
phyloseq::tax_table(aob.rare.1282.seq)
# merge taxa by species
aob.sp <- tax_glom(aob.rare.1282.seq, taxrank = "Species", NArm = F)
aob.sp.ra <- transform_sample_counts(aob.sp, function(x) x/sum(x))
sample_data(aob.sp.ra)

aob.sp.df <- psmelt(aob.sp.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

test.df <- psmelt(aob.sp.ra) %>%
  group_by(Type, Date, OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

aob.abund.trt.subcla <- psmelt(aob.gen.ra) %>%
  group_by(Type,Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)


colours <- ColourPalleteMulti(aob.sp.df, "Genus", "Species")
colours
#install.packages("Polychrome")
#install.packages("colorBlindness")
library(colorBlindness)
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
color <- c("#260F99","#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF",
           "#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF",
           "#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2",
           "#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2","#99540F", "#B26F2C")
displayAvailablePalette(color="white")
displayAvailablePalette()
SteppedSequential5Steps
str(aob.sp.df)
#install.packages("ggh4x")
library(ggh4x)

aob.sp.df$Type <- factor(aob.sp.df$Type, levels = c("BS", "RS"),
                         labels = c("Bulk Soil", "Rhizosphere"))
aob.sp.df$Treatment <- factor(aob.sp.df$Treatment, levels = c("D", "K", "M"),
                              labels = c("Biodynamic", "Conventional", "Mineral"))
legend <- "AOB Taxa"
library(scales)
#x_cols <- rep(hue_pal()(length(unique(interaction(aob.sp.df$Date, aob.sp.df$Irrigation)))))
#aob.sp.df$Date <- factor(aob.sp.df$Date, levels = unique(aob.sp.df$Date))
aob.sp.plot <- ggplot(aob.sp.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(legend, values=SteppedSequential5Steps)+
  facet_nested(~Type+Treatment, 
               nest_line = element_line(linetype = 1, linewidth = 1), 
               scales="free",
               resect = unit(0.4, "cm"))+
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
  scale_y_continuous(expand = c(0,0))+ guides(x="axis_nested")
  #geom_vline(xintercept = 5, linetype="dotted", 
             #color = "blue", linewidth=1.5)
aob.sp.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB_meanRA_barplot2.eps",
       aob.sp.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)

#________________________________________________________________________________________________________________________________-
# genus level

aob.gen <- tax_glom(aob.rare.1282.seq, taxrank = "Genus", NArm = F)
aob.gen.ra <- transform_sample_counts(aob.gen, function(x) x/sum(x))
aob.genus.df <- psmelt(aob.gen.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

str(aob.genus.df)
view(aob.genus.df)
aob.genus.df$Type <- factor(aob.genus.df$Type, levels = c("BS", "RS"),
                         labels = c("Bulk Soil", "Rhizosphere"))
aob.genus.df$Treatment <- factor(aob.genus.df$Treatment, levels = c("D", "K", "M"),
                              labels = c("BIODYN", "CONFYM", "CONMIN"))
aob.genus.df$Irrigation <- factor(aob.genus.df$Irrigation, levels = c("Rainout", "Control"),
                         labels = c("Drought", "Control"))
aob.genus.df$Date <- factor(aob.genus.df$Date, levels = c("Apr-28", "Jun-01", "Jul-05", "Jul-20", "Sep-13"))
                          #labels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))

library(RColorBrewer)

brewer.pal(12, "Paired")
display.brewer.pal(n, name)
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
colorblindFriendly=T)

set.seed(13)
aob.genus.plot <- ggplot(aob.genus.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  theme_bw()+
  #scale_fill_manual(legend, values=c("#1B9E77","#D95F02","#7570B3","#E7298A"))+
  scale_fill_manual(values=c("#33A02C","#FF7F00","#6A3D9A","#1F78B4"))+
  facet_nested(~Type+Treatment,
               #nest_line = element_line(linetype = 1, linewidth = 0.5), 
               scales="free")+
               #resect = unit(5, "cm"))+
  guides(fill=guide_legend(ncol=1))+
  labs(y= "Mean Relative Abundance", title = "A. AOB")+
  theme(
        legend.direction = "vertical",
        legend.position="right",
        plot.title = element_text(size = 25, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=20),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y =element_text(size=22),
        legend.text=element_text(size = 20),
        legend.title = element_text(size=22, face="bold"),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        #strip.background = element_blank(),
        strip.background=element_rect(color="grey30", fill="white"),
        strip.text.x = element_text(size = 20),
        #panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2),
        panel.border=element_rect(color="grey30",fill = NA))+
  scale_y_continuous(expand = c(0,0))+ guides(x="axis_nested")
  #geom_vline(xintercept = 5, linetype="dotted", 
             #color = "blue", linewidth=1.5)
aob.genus.plot

# Combine All Relative Abundance Figures
Composition.All <- (aob.genus.plot / aoa.clade.plot / com.clade.plot) & 
 theme(legend.justification = "left")
Composition.All
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Supp.Fig.4dpi3002.tiff",
       Composition.All, device = "tiff",
       width = 18, height = 17, 
       units= "in", dpi = 600, compression="lzw")

#________________________________________________________________________________________________________________________________
#setwd('D:/Fina/INRAE_Project/microservices_fig/AOB'

# 2. unrarefied data
head(aob.tax)
aob.tax.physeq = phyloseq::tax_table(as.matrix(aob.tax)) # taxonomy table
dim(aob.tax.physeq)
# phyloseq object of the metadata
str(meta_micro)
meta_micro.row <- column_to_rownames(meta_micro, var = "SampleID")

aob.meta.physeq <- sample_data(meta_micro.row)# meta data
sample_names(aob.meta.physeq)

# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB-rooted-tree/')
AOB_rooted_tree <- ape::read.tree("tree.nwk")
AOB_rooted_tree
# make phyloseq object
aob.physeq.unrare <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq,AOB_rooted_tree)
aob.physeq.unrare
sample_data(aob.physeq.unrare)$SampleID <- paste0("S", sample_data(aob.physeq.unrare)$SampleID)
sample_data(aob.physeq.unrare)
# merge taxa by species
aob.sp.unrare <- tax_glom(aob.physeq.unrare, taxrank = "Species", NArm = F)
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

#aob.sp.unrare.df$Date <- as.Date(aob.sp.unrare.df$Date , "%m/%d/%Y")
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
  scale_fill_manual(legend, values=SteppedSequential5Steps)+
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
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_meanRA_unrare_barplot.tiff",
       aob.sp.unrare.plot, device = "tiff",
       width = 15, height =6, 
       units= "in", dpi = 600)














