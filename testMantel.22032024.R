
### Mantel correlation analysis with Pearson's correlation method

# Project: Microservices, BIODIVERSA
# Author: Ari Fina Bintarti
# Date: 22/03/2024

# Data set: AOB, AOA, Comammox, Soil properties, NO3, NH4, N2O (Remove TS data)
# Note: This analysis was performed using a data set only until July 20th (without the last sampling point), to match the N2O data set

library(dplyr)
library(vegan)
library(ggplot2)
library(devtools)
#install.packages("microeco")
library(microeco)
library(magrittr)
###devtools::install_github("hannet91/ggcor")
#devtools::install_github("houyunhuang/ggcor")
#.rs.restartR()
#install.packages("igraph")
library(ggcor)


# 1. CONTROL

# Preparing the micro table class
# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
# select only data set until Jul 20th
aob.BS.sub  <- aob.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
aob.BS.sub1 <- prune_taxa(taxa_sums(aob.BS.sub)>0, aob.BS.sub)
aob.BS.sub1
# subset AOB CONTROL
aob.BS.cont_seq <- subset_samples(aob.BS.sub1, Irrigation=="Control")
aob.BS.cont_seq1 <- prune_taxa(taxa_sums(aob.BS.cont_seq)>0, aob.BS.cont_seq)
sort(rowSums(otu_table(aob.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aob.BS.cont_df <- as.data.frame(otu_table(aob.BS.cont_seq1))
dim(aob.BS.cont_df) # 616 otus, 47 samples
view(aob.BS.cont_df)
# taxonomy table AOB CONTROL
aob.BS.cont_tax <- as.data.frame(phyloseq::tax_table(aob.BS.cont_seq1))
dim(aob.BS.cont_tax) # 616 2
aob.BS.cont_tax$Group <- 'AOB'#adding a new column with AOB as value
view(aob.BS.cont_tax)
# environmental data
aob.env.Jul <- aob.meta.df.sub[1:95,] # no S11 already
colnames(aob.env.Jul)
dim(aob.env.Jul) # 95 48
aob.env.Jul <- aob.env.Jul %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOB_richness" = "Richness",
                 "AOB_Shannon" = "Shannon")
# adding N2O data set
setwd('/Users/arifinabintarti/Documents/France/microservices/')
soilprop_n2o <- read.csv("soilprop.n2o.csv")
colnames(soilprop_n2o)
soilprop_n2o.sub <- soilprop_n2o %>% filter(SampleID != "S11")# filter out S11 from the metadata
aob.env.Jul$mean.N2Oflux <- soilprop_n2o.sub$mean.N2Oflux
aob.env.Jul$N2Oflux.cbrt <- soilprop_n2o.sub$N2Oflux.cbrt
aob.env.Jul$log10NH4 <- soilprop_n2o.sub$log10NH4
aob.env.Jul$logNO3 <- soilprop_n2o.sub$logNO3
view(aob.env.Jul)
# relocate N2O flux
#aob.env.Jul <- aob.env.Jul %>% relocate("mean.N2Oflux", .after="AOA_Shannon")
#aob.env.Jul <- aob.cont_env.ed4 %>% relocate("AOB_Shannon", .after="AOB_richness")
# environmental data CONTROL
view(aob.env.Jul)
aob.BS.cont_env <- aob.env.Jul %>% filter(Irrigation == "Control")
colnames(aob.BS.cont_env)
view(aob.BS.cont_env)
dim(aob.BS.cont_env) # 47 48
#extract alpha diversity CONTROL
aob.BS.cont_alpha <- aob.BS.cont_env[,45:46]
head(aob.BS.cont_alpha)
aob.BS.cont_alpha <- rownames_to_column(aob.BS.cont_alpha, var = "SampleID")
dim(aob.BS.cont_alpha) # 47 3
# adding qPCR data
setwd('/Users/arifinabintarti/Documents/France/microservices/')
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
# adding qPCR data CONTROL
qpcr.BS.control <- qpcr.BS %>% filter(irrigation == "control")
view(qpcr.BS.control)
colnames(qpcr.BS.control)
# filter only data set until Jul20
qpcr.BS.cont.jul <- qpcr.BS.control %>% filter(sampling.date != "Sept 13th")
view(qpcr.BS.cont.jul)
# select only log gram/dws AOB, AOA, Coma A, and Coma B
qpcr.BS.cont_logdws <- qpcr.BS.cont.jul[,29:32]
dim(qpcr.BS.cont_logdws) # 48 4
view(qpcr.BS.cont_logdws)
qpcr.BS.cont_logdws.x <- rownames_to_column(qpcr.BS.cont_logdws, var = "SampleID")
aob.BS.cont_logdws <- qpcr.BS.cont_logdws.x %>% filter(SampleID != "S11")# filter out S11 from the metadata
aob.BS.cont_logdws <- column_to_rownames(aob.BS.cont_logdws, var = "SampleID")
dim(aob.BS.cont_logdws) # 47 4
view(aob.BS.cont_logdws)
# combine ONLY CONTROL
aob.cont_env.ed <- cbind(aob.BS.cont_env,aob.BS.cont_logdws)
head(aob.cont_env.ed)
dim(aob.cont_env.ed) # 47 56
str(aob.cont_env.ed)
# rename CONTROL environment
aob.cont_env.ed <- aob.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names
aob.cont_env.ed$UniqueID <- factor(aob.cont_env.ed$UniqueID)
aob.cont_env.ed$PlotID <- factor(aob.cont_env.ed$PlotID)
aob.cont_env.ed$Block <- factor(aob.cont_env.ed$Block)
aob.cont_env.ed$Irrigation <- factor(aob.cont_env.ed$Irrigation)
aob.cont_env.ed$Treatment <- factor(aob.cont_env.ed$Treatment)
aob.cont_env.ed$Type <- factor(aob.cont_env.ed$Type)
aob.cont_env.ed$rep <- factor(aob.cont_env.ed$rep)
aob.cont_env.ed$rep2 <- factor(aob.cont_env.ed$rep2)
aob.cont_env.ed$Date <- factor(aob.cont_env.ed$Date)
aob.cont_env.ed$x <- factor(aob.cont_env.ed$x)
aob.cont_env.ed$var3 <- factor(aob.cont_env.ed$var3)
aob.cont_env.ed$var2 <- factor(aob.cont_env.ed$var2)

# adding alpha data of AOA CONTROL
dim(aoa.BS.cont_alpha) # 48 2 still contain S11
view(aoa.BS.cont_alpha)
#aoa.BS.cont_alpha <- rownames_to_column(aoa.BS.cont_alpha, var = "SampleID")
aoa.BS.cont_alpha.ed <- aoa.BS.cont_alpha %>% filter(SampleID != "S11") # Remove S11 to match the OTU table
dim(aoa.BS.cont_alpha.ed) #47 3
head(aoa.BS.cont_alpha.ed)
aoa.BS.cont_alpha.ed.x <- column_to_rownames(aoa.BS.cont_alpha.ed, var = "SampleID")
aob.cont_env.ed2 <- cbind(aob.cont_env.ed,aoa.BS.cont_alpha.ed.x)
dim(aob.cont_env.ed2) # 47 58
head(aob.cont_env.ed2)

# adding alpha data of COMA CONTROL
aob.cont_env.ed2 <- rownames_to_column(aob.cont_env.ed2, var = "SampleID")
head(com.BS.cont_alpha)
dim(com.BS.cont_alpha) # 48 3 still contain S11
com.BS.cont_alpha.ed <- com.BS.cont_alpha %>% filter(SampleID != "S11") # Remove S11 to match the OTU table
dim(com.BS.cont_alpha.ed) # 47 3
head(com.BS.cont_alpha.ed)
#com.BS.cont_alpha.x <- rownames_to_column(com.BS.cont_alpha, var = "SampleID")
aob.cont_env.ed3 <- aob.cont_env.ed2 %>%
 left_join(com.BS.cont_alpha.ed, by="SampleID") 
dim(aob.cont_env.ed3) #47 61
colnames(aob.cont_env.ed3)
view(aob.cont_env.ed3)
# relocate alpha data of AOB in CONTROL data
aob.cont_env.ed4 <- aob.cont_env.ed3 %>% relocate("AOB_richness", .after="AOA_Shannon") %>%
                                         relocate("AOB_Shannon", .after="AOB_richness") %>%
                                         relocate("N2O", .after="GWC")
                                         
colnames(aob.cont_env.ed4)
view(aob.cont_env.ed4) # there are nno NAs in the dataset
dim(aob.cont_env.ed4) #47 61
aob.cont_env.ed4 <- column_to_rownames(aob.cont_env.ed4, var = "SampleID")
# create  a microtable CONTROL
aob.cont_microdata <- microtable$new(sample_table = aob.cont_env.ed4, otu_table = aob.BS.cont_df, tax_table = aob.BS.cont_tax)
# calculate beta diversity in CONTROL
aob.cont_microdata$tidy_dataset()
set.seed(13)
aob.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test in CONTROL
aob.t1.cont <- trans_env$new(dataset = aob.cont_microdata, env_cols = c(17,23:24,27,31:33,48:60))
set.seed(13)
aob.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
aob.t1.cont$res_mantel
# extract a part of the results in CONTROL
aob.x1.cont <- data.frame(aob.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aob.x1.cont[aob.x1.cont=="All"] <- "AOB"
# rename columns in CONTROL
colnames(aob.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(aob.x1.cont)

################################################################################
# rarefied AOA ASV table of bulk soil CONTROL

aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# select only data set until Jul 20th
aoa.BS.sub  <- aoa.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
aoa.BS.sub1 <- prune_taxa(taxa_sums(aoa.BS.sub)>0, aoa.BS.sub)
aoa.BS.sub1 # 445 taxa, 96 samples
# subset AOA CONTROL
aoa.BS.cont_seq <- subset_samples(aoa.BS.sub1, Irrigation=="Control")
aoa.BS.cont_seq1 <- prune_taxa(taxa_sums(aoa.BS.cont_seq)>0, aoa.BS.cont_seq)
sort(rowSums(otu_table(aoa.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aoa.BS.cont_df <- as.data.frame(otu_table(aoa.BS.cont_seq1))
dim(aoa.BS.cont_df) # 347 otus, 48 samples nothing is missing
view(aoa.BS.cont_df)
# taxonomy table AOA CONTROL
aoa.BS.cont_tax <- as.data.frame(phyloseq::tax_table(aoa.BS.cont_seq1))
dim(aoa.BS.cont_tax) # 347 6
aoa.BS.cont_tax$Group <- 'AOA'#adding a new column with AOB as value
view(aoa.BS.cont_tax)
# environmental data
aoa.env.Jul <- aoa.meta.df[1:96,] 
view(aoa.env.Jul)
dim(aoa.env.Jul) # 96 48
aoa.env.Jul <- aoa.env.Jul %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOA_richness" = "Richness",
                 "AOA_Shannon" = "Shannon")
# adding N2O data set
setwd('/Users/arifinabintarti/Documents/France/microservices/')
soilprop_n2o <- read.csv("soilprop.n2o.csv")
colnames(soilprop_n2o)
aoa.env.Jul$mean.N2Oflux <- soilprop_n2o$mean.N2Oflux
aoa.env.Jul$N2Oflux.cbrt <- soilprop_n2o$N2Oflux.cbrt
aoa.env.Jul$log10NH4 <- soilprop_n2o$log10NH4
aoa.env.Jul$logNO3 <- soilprop_n2o$logNO3
view(aoa.env.Jul)
# environmental data CONTROL
aoa.BS.cont_env <- aoa.env.Jul %>% filter(Irrigation == "Control")
colnames(aoa.BS.cont_env)
view(aoa.BS.cont_env)
dim(aoa.BS.cont_env) # 48, 52
#extract alpha diversity CONTROL
aoa.BS.cont_alpha <- aoa.BS.cont_env[,45:46]
head(aoa.BS.cont_alpha)
aoa.BS.cont_alpha <- rownames_to_column(aoa.BS.cont_alpha, var = "SampleID")
dim(aoa.BS.cont_alpha) # 48 3
# adding qPCR data CONTROL
view(qpcr.BS.cont_logdws)
aoa.BS.cont_logdws <- qpcr.BS.cont_logdws
dim(aoa.BS.cont_logdws) # 48 6
view(aoa.BS.cont_logdws)
# combine ONLY CONTROL
aoa.cont_env.ed <- cbind(aoa.BS.cont_env,aoa.BS.cont_logdws)
head(aoa.cont_env.ed)
dim(aoa.cont_env.ed) # 48 56
str(aoa.cont_env.ed)
# rename CONTROL environment
aoa.cont_env.ed <- aoa.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names
aoa.cont_env.ed$UniqueID <- factor(aoa.cont_env.ed$UniqueID)
aoa.cont_env.ed$PlotID <- factor(aoa.cont_env.ed$PlotID)
aoa.cont_env.ed$Block <- factor(aoa.cont_env.ed$Block)
aoa.cont_env.ed$Irrigation <- factor(aoa.cont_env.ed$Irrigation)
aoa.cont_env.ed$Treatment <- factor(aoa.cont_env.ed$Treatment)
aoa.cont_env.ed$Type <- factor(aoa.cont_env.ed$Type)
aoa.cont_env.ed$rep <- factor(aoa.cont_env.ed$rep)
aoa.cont_env.ed$rep2 <- factor(aoa.cont_env.ed$rep2)
aoa.cont_env.ed$Date <- factor(aoa.cont_env.ed$Date)
aoa.cont_env.ed$x <- factor(aoa.cont_env.ed$x)
aoa.cont_env.ed$var3 <- factor(aoa.cont_env.ed$var3)
aoa.cont_env.ed$var2 <- factor(aoa.cont_env.ed$var2)

# adding alpha data of AOB of CONTROL
dim(aob.BS.cont_alpha) # 47 3, it does not contain S11
head(aob.BS.cont_alpha)
head(aoa.cont_env.ed)
aoa.cont_env.ed <- rownames_to_column(aoa.cont_env.ed, var = "SampleID")
aoa.cont_env.ed2 <- aoa.cont_env.ed %>%
 left_join(aob.BS.cont_alpha, by="SampleID") # USE left join for different length of two dataframes
colnames(aoa.cont_env.ed2)
dim(aoa.cont_env.ed2) # 48 55, BUT THERE IS A MISSING DATA ON AOB ALPHA (S11)
view(aoa.cont_env.ed2) 
# adding alpha data of COMA of CONTROL
dim(com.BS.cont_alpha) # 48 3
aoa.cont_env.ed3 <- aoa.cont_env.ed2 %>%
 left_join(com.BS.cont_alpha, by="SampleID") 
dim(aoa.cont_env.ed3) # 48 61
view(aoa.cont_env.ed3)
# relocate alpha data of AOA in CONTROL
colnames(aoa.cont_env.ed3)
aoa.cont_env.ed4 <- aoa.cont_env.ed3 %>% relocate("AOA_richness", .after="Coma_B_abundance") %>%
                                         relocate("AOA_Shannon", .after="AOA_richness") %>%
                                         relocate("N2O", .after="GWC")
head(aoa.cont_env.ed4)
view(aoa.cont_env.ed4)
dim(aoa.cont_env.ed4) # 48 61, with S11 missing AOB alpha
aoa.cont_env.ed4 <- column_to_rownames(aoa.cont_env.ed4, var = "SampleID")
# create  a microtable for CONTROL
aoa.cont_microdata <- microtable$new(sample_table = aoa.cont_env.ed4, otu_table = aoa.BS.cont_df, tax_table = aoa.BS.cont_tax)
# calculate beta diversity for CCONTROL
aoa.cont_microdata$tidy_dataset()
set.seed(13)
aoa.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test for CONTROL
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
aoa.t1.cont <- trans_env$new(dataset = aoa.cont_microdata, env_cols = c(17,23:24,27,31:33,48:60),complete_na = TRUE) 
set.seed(13)
aoa.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
aoa.t1.cont$res_mantel
# extract a part of the results CONTROL
aoa.x1.cont <- data.frame(aoa.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aoa.x1.cont[aoa.x1.cont=="All"] <- "AOA"
# rename columns CONTROL
colnames(aoa.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(aoa.x1.cont)

################################################################################
# rarefied COMAMMOX ASV table of bulk soil CONTROL

com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
dim(com.asv.tab) # 497, 118 samples (missing 2 values: SS26 & S52)
# select only data set until Jul 20th
com.BS.sub  <- com.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
com.BS.sub1 <- prune_taxa(taxa_sums(com.BS.sub)>0, com.BS.sub)
com.BS.sub1 # 458 taxa, 94 samples
# subset AOA CONTROL
com.BS.cont_seq <- subset_samples(com.BS.sub1, Irrigation=="Control")
com.BS.cont_seq1 <- prune_taxa(taxa_sums(com.BS.cont_seq)>0, com.BS.cont_seq)
sort(rowSums(otu_table(com.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
com.BS.cont_df <- as.data.frame(otu_table(com.BS.cont_seq1))
dim(com.BS.cont_df) # 335 otus, 48 samples (COMPLETE because SS26 & S52 are not in the CONTROL)
view(com.BS.cont_df)
# taxonomy table COM CONTROL
com.BS.cont_tax <- as.data.frame(phyloseq::tax_table(com.BS.cont_seq1))
dim(com.BS.cont_tax) # 335 4
com.BS.cont_tax$Group <- 'COM'#adding a new column with AOB as value
view(com.BS.cont_tax)
# environmental data
com.env.Jul <- com.meta.df[1:94,] # 
view(com.env.Jul)
dim(com.env.Jul) # 95 48
com.env.Jul <- com.env.Jul %>% 
  dplyr::rename("UniqueID"="SampleID",
                "COM_richness" = "Richness",
                 "COM_Shannon" = "Shannon")
# adding N2O data set
setwd('/Users/arifinabintarti/Documents/France/microservices/')
soilprop_n2o <- read.csv("soilprop.n2o.csv")
colnames(soilprop_n2o)
soilprop_n2o.sub2 <- soilprop_n2o %>% filter(SampleID != "S26",SampleID != "S52")# filter out S26 & S52 from the metadata
com.env.Jul$mean.N2Oflux <- soilprop_n2o.sub2$mean.N2Oflux
com.env.Jul$N2Oflux.cbrt <- soilprop_n2o.sub2$N2Oflux.cbrt
com.env.Jul$log10NH4 <- soilprop_n2o.sub2$log10NH4
com.env.Jul$logNO3 <- soilprop_n2o.sub2$logNO3
dim(com.env.Jul) # 94 52
# environmental data CONTROL
com.BS.cont_env <- com.env.Jul %>% filter(Irrigation == "Control")
colnames(com.BS.cont_env)
view(com.BS.cont_env)
dim(com.BS.cont_env) # 48 52 
#extract alpha diversity CONTROL
com.BS.cont_alpha <- com.BS.cont_env[,45:46]
head(com.BS.cont_alpha)
com.BS.cont_alpha <- rownames_to_column(com.BS.cont_alpha, var = "SampleID")
dim(com.BS.cont_alpha) # 48 3 (complete)
# adding qPCR data CONTROL
dim(qpcr.BS.cont_logdws)
com.BS.cont_logdws <- qpcr.BS.cont_logdws
# combine ONLY CONTROL
com.cont_env.ed <- cbind(com.BS.cont_env,com.BS.cont_logdws)
head(com.cont_env.ed)
dim(com.cont_env.ed) # 48 56
str(com.cont_env.ed)
# rename CONTROL environment
com.cont_env.ed <- com.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names

com.cont_env.ed$UniqueID <- factor(com.cont_env.ed$UniqueID)
com.cont_env.ed$PlotID <- factor(com.cont_env.ed$PlotID)
com.cont_env.ed$Block <- factor(com.cont_env.ed$Block)
com.cont_env.ed$Irrigation <- factor(com.cont_env.ed$Irrigation)
com.cont_env.ed$Treatment <- factor(com.cont_env.ed$Treatment)
com.cont_env.ed$Type <- factor(com.cont_env.ed$Type)
com.cont_env.ed$rep <- factor(com.cont_env.ed$rep)
com.cont_env.ed$rep2 <- factor(com.cont_env.ed$rep2)
com.cont_env.ed$Date <- factor(com.cont_env.ed$Date)
com.cont_env.ed$x <- factor(com.cont_env.ed$x)
com.cont_env.ed$var3 <- factor(com.cont_env.ed$var3)
com.cont_env.ed$var2 <- factor(com.cont_env.ed$var2)

# adding alpha data of AOA CONTROL
head(aoa.BS.cont_alpha)
dim(aoa.BS.cont_alpha) # 48 3 complete
aoa.BS.cont_alpha <- column_to_rownames(aoa.BS.cont_alpha, var = "SampleID")
com.cont_env.ed2 <- cbind(com.cont_env.ed,aoa.BS.cont_alpha)
dim(com.cont_env.ed2) # 48 58
# adding alpha data of AOB CONTROL
com.cont_env.ed2 <- rownames_to_column(com.cont_env.ed2, var = "SampleID")
dim(aob.BS.cont_alpha) # 47 3 because it does not contain S11,thus use left_join
com.cont_env.ed3 <- com.cont_env.ed2 %>%
 left_join(aob.BS.cont_alpha, by="SampleID") 
dim(com.cont_env.ed3) # 48 61, but the AOB alpha of S11 are missing
view(com.cont_env.ed3)
# relocate alpha data of COMA CONTROL
com.cont_env.ed4 <- com.cont_env.ed3 %>% relocate("COM_richness", .after="AOB_Shannon") %>% 
                                         relocate("COM_Shannon", .after="COM_richness") %>%
                                         relocate("N2O", .after="GWC")
dim(com.cont_env.ed4)
str(com.cont_env.ed4) # 48 61
com.cont_env.ed4 <- column_to_rownames(com.cont_env.ed4, var = "SampleID")

# create  a microtable CONTROL
com.cont_microdata <- microtable$new(sample_table = com.cont_env.ed4, otu_table = com.BS.cont_df, tax_table = com.BS.cont_tax)
# calculate beta diversity CONTROL
com.cont_microdata$tidy_dataset()
set.seed(13)
com.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test CONTROL
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
com.t1.cont <- trans_env$new(dataset = com.cont_microdata, env_cols = c(17,23:24,27,31:33,48:60), complete_na = TRUE)
set.seed(13)
com.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
com.t1.cont$res_mantel
# extract a part of the results CONTROL
com.x1.cont <- data.frame(com.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
com.x1.cont[com.x1.cont=="All"] <- "comammox"
# rename columns CONTROL
colnames(com.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(com.x1.cont)

###############################################################################
# generate interval data CONTROL
aob.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
aoa.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
com.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# combine three tables CONTROL
plot_table.CONT <- rbind(aob.x1.cont, aoa.x1.cont,com.x1.cont)
plot_table.CONT <- plot_table.CONT %>% 
  mutate(spec= replace(spec, spec == "AOB", "AOB β diversity")) %>%
  mutate(spec= replace(spec, spec =="AOA","AOA β diversity")) %>%
  mutate(spec= replace(spec, spec=="comammox", "Comammox β diversity")) %>%# change column names
  mutate(env= replace(env, env =="Total_C","Total C")) %>%
  mutate(env= replace(env, env =="Total_N","Total N")) %>%
  mutate(env= replace(env, env =="AOA_abundance","AOA abundance")) %>%
  mutate(env= replace(env, env =="AOB_abundance","AOB abundance")) %>%
  mutate(env= replace(env, env =="Coma_A_abundance","COM-A abundance")) %>%
  mutate(env= replace(env, env =="Coma_B_abundance","COM-B abundance")) %>%
  mutate(env= replace(env, env =="AOA_richness","AOA richness")) %>%
  mutate(env= replace(env, env =="AOA_Shannon","AOA Shannon")) %>%
  mutate(env= replace(env, env =="AOB_richness","AOB richness")) %>%
  mutate(env= replace(env, env =="AOB_Shannon","AOB Shannon")) %>%
  mutate(env= replace(env, env =="COM_richness","COM richness")) %>%
  mutate(env= replace(env, env =="COM_Shannon","COM Shannon")) 
view(plot_table.CONT)
# plotting
#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
aoa.env.cont <- aoa.t1.cont$data_env # the complete one
dim(aoa.env.cont) # 48 18
com.env.cont <- com.t1.cont$data_env # this one is also complete (include all samples and no missing data as a result of interpolation by MICE)
view(com.env.cont)
#write.csv(aob.env, file = "aob.env.mantel.csv")
all.env.cont <- aoa.env.cont
#all.env.no.alpha <- all.env[,c(1:4,6:11)]
colnames(all.env.cont)
str(all.env.cont)
#qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
#colnames(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
#qpcr.BS.logdws <- qpcr.BS[,29:32]
#all.env.qpcr <- cbind(all.env.no.alpha,qpcr.BS.logdws)
#all.env.qpcr$Mg_mgkg <- as.numeric(all.env.qpcr$Mg_mgkg)
all.env.cont.ed <- all.env.cont %>% 
  dplyr::rename("AOA abundance" = "AOA_abundance",
         "AOB abundance" = "AOB_abundance",
         "COM-A abundance" = "Coma_A_abundance",
         "COM-B abundance" = "Coma_B_abundance",
         "AOA richness" = "AOA_richness",
         "AOA Shannon" = "AOA_Shannon",
         "AOB richness" = "AOB_richness",
         "AOB Shannon" = "AOB_Shannon",
         "COM richness" = "COM_richness",
         "COM Shannon" = "COM_Shannon",
         "Total C"="Total_C",
         "Total N"="Total_N") 
str(all.env.cont.ed)
head(all.env.cont.ed)
# combine all alpha diversity
#ao.alpha <- aoa.alpha %>%
 #left_join(aob.alpha, by="SampleID") %>%
 #left_join(com.alpha, by="SampleID")
#ao.alpha <- column_to_rownames(ao.alpha, var="SampleID")
#head(ao.alpha)
# combine
#all.env.qpcr.ed2 <- cbind(all.env.qpcr.ed,ao.alpha)

mantelplot.control <- quickcor(all.env.cont.ed, type = "upper", cor.test = TRUE,method="pearson",exact=FALSE) + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.CONT, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank())
mantelplot.control
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("mantel.control_withN2O.tiff",
       mantelplot.control, device = "tiff",
       width = 10, height = 10, 
       units= "in", dpi = 600)

# try with NA of AOB alpha
all.env.cont.ed2 <- rownames_to_column(all.env.cont.ed, var = "SampleID")
all.env.cont.ed3 <- all.env.cont.ed2 %>% 
  mutate(`AOB richness`=ifelse(SampleID=="S11",NA,`AOB richness`),
         `AOB Shannon`=ifelse(SampleID=="S11",NA,`AOB Shannon`))
view(all.env.cont.ed3)
all.env.cont.ed3 <- column_to_rownames(all.env.cont.ed3, var = "SampleID")

mantelplot.control.noAOBalph <- quickcor(all.env.cont.ed3, type = "upper", cor.test = TRUE,method="pearson",exact=FALSE,use="pairwise.complete.obs") + #cor.test = TRUE #use="complete.obs"
  geom_square(position = "identity") + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.CONT, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 3) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank(),
        axis.text= element_text(size=15),
        plot.margin = margin(0, 0, 0, 0, "cm"))
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.control.noAOBalph

ggsave("mantelplot.control.noAOBalph_Pearson.tiff",
       mantelplot.control.noAOBalph, device = "tiff",
       width = 13, height = 10, 
       units= "in", dpi = 600)

###################################################################################################################

# 2. DROUGHT

# Preparing the microtable class
# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
aob.asv.tab <- as.data.frame(otu_table(aob.physeq_bulk1))
dim(aob.asv.tab)
# select only data set until Jul 20th
aob.BS.sub  <- aob.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
aob.BS.sub1 <- prune_taxa(taxa_sums(aob.BS.sub)>0, aob.BS.sub)
aob.BS.sub1
# subset AOB DROUGHT
aob.BS.dro_seq <- subset_samples(aob.BS.sub1, Irrigation=="Rainout")
aob.BS.dro_seq1 <- prune_taxa(taxa_sums(aob.BS.dro_seq)>0, aob.BS.dro_seq)
sort(rowSums(otu_table(aob.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aob.BS.dro_df <- as.data.frame(otu_table(aob.BS.dro_seq1))
dim(aob.BS.dro_df) # 600 otus, 48 samples
view(aob.BS.dro_df)
# taxonomy table AOB DROUGHT
aob.BS.dro_tax <- as.data.frame(phyloseq::tax_table(aob.BS.dro_seq1))
dim(aob.BS.dro_tax) # 600 2
aob.BS.dro_tax$Group <- 'AOB'#adding a new column with AOB as value
view(aob.BS.dro_tax)
# environmental data
aob.env.Jul 
# environmental data DROUGHT
aob.BS.dro_env <- aob.env.Jul %>% filter(Irrigation == "Rainout")
colnames(aob.BS.dro_env)
view(aob.BS.dro_env)
dim(aob.BS.dro_env) # 48 52
#extract alpha diversity DROUGHT
aob.BS.dro_alpha <- aob.BS.dro_env[,45:46]
head(aob.BS.dro_alpha)
aob.BS.dro_alpha <- rownames_to_column(aob.BS.dro_alpha, var = "SampleID")
dim(aob.BS.dro_alpha) # 48 3
# adding qPCR data DROUGHT
qpcr.BS.drought <- qpcr.BS %>% filter(irrigation == "rainout")
view(qpcr.BS.drought)
colnames(qpcr.BS.drought)
# filter only data set until Jul20
qpcr.BS.dro.jul <- qpcr.BS.drought %>% filter(sampling.date != "Sept 13th")
colnames(qpcr.BS.dro.jul)
qpcr.BS.dro_logdws <- qpcr.BS.dro.jul[,29:32]
dim(qpcr.BS.dro_logdws) # 48 4
aob.BS.dro_logdws <- qpcr.BS.dro_logdws
head(aob.BS.dro_logdws)
# combine ONLY CONTROL
aob.dro_env.ed <- cbind(aob.BS.dro_env,aob.BS.dro_logdws)
head(aob.dro_env.ed)
dim(aob.dro_env.ed) # 48 56
str(aob.dro_env.ed)
# rename CONTROL environment
aob.dro_env.ed <- aob.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names
aob.dro_env.ed$UniqueID <- factor(aob.dro_env.ed$UniqueID)
aob.dro_env.ed$PlotID <- factor(aob.dro_env.ed$PlotID)
aob.dro_env.ed$Block <- factor(aob.dro_env.ed$Block)
aob.dro_env.ed$Irrigation <- factor(aob.dro_env.ed$Irrigation)
aob.dro_env.ed$Treatment <- factor(aob.dro_env.ed$Treatment)
aob.dro_env.ed$Type <- factor(aob.dro_env.ed$Type)
aob.dro_env.ed$rep <- factor(aob.dro_env.ed$rep)
aob.dro_env.ed$rep2 <- factor(aob.dro_env.ed$rep2)
aob.dro_env.ed$Date <- factor(aob.dro_env.ed$Date)
aob.dro_env.ed$x <- factor(aob.dro_env.ed$x)
aob.dro_env.ed$var3 <- factor(aob.dro_env.ed$var3)
aob.dro_env.ed$var2 <- factor(aob.dro_env.ed$var2)

# adding alpha data of AOA DROUGHT
dim(aoa.BS.dro_alpha) # 48 3 complete
view(aoa.BS.dro_alpha)
aoa.BS.dro_alpha <- column_to_rownames(aoa.BS.dro_alpha, var = "SampleID")
aob.dro_env.ed2 <- cbind(aob.dro_env.ed,aoa.BS.dro_alpha)
dim(aob.dro_env.ed2) # 48 58
head(aob.dro_env.ed2)
# adding alpha data of COMA DROUGHT
aob.dro_env.ed2 <- rownames_to_column(aob.dro_env.ed2, var = "SampleID")
head(com.BS.dro_alpha)
dim(com.BS.dro_alpha) # 46 3 it does not contain S26 and S52 (missing)
aob.dro_env.ed3 <- aob.dro_env.ed2 %>%
 left_join(com.BS.dro_alpha, by="SampleID") 
dim(aob.dro_env.ed3) # 48 61 but it has NA in the S26 and S52 for Coma alpha
head(aob.dro_env.ed3)
# relocate alpha data of AOB DROUGHT data
aob.dro_env.ed4 <- aob.dro_env.ed3 %>% relocate("AOB_richness", .after="AOA_Shannon")%>%
 relocate("AOB_Shannon", .after="AOB_richness") %>% 
 relocate("N2O", .after="GWC")
view(aob.dro_env.ed4) # there are no NAs in the dataset
dim(aob.dro_env.ed4) #48 61 contains missing values (NA)
aob.dro_env.ed4 <- column_to_rownames(aob.dro_env.ed4, var = "SampleID")
# create  a microtable DROUGHT
aob.dro_microdata <- microtable$new(sample_table = aob.dro_env.ed4, otu_table = aob.BS.dro_df, tax_table = aob.BS.dro_tax)
# calculate beta diversity in DROUGHT
aob.dro_microdata$tidy_dataset()
set.seed(13)
aob.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test in DROUGHT
set.seed(13)
aob.t1.dro <- trans_env$new(dataset = aob.dro_microdata, env_cols = c(17,23:24,27,31:33,48:60),complete_na = TRUE) # interpolation for the missing values
set.seed(13)
aob.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = T,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
aob.t1.dro$res_mantel
# extract a part of the results in DROUGHT
aob.x1.dro <- data.frame(aob.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aob.x1.dro[aob.x1.dro=="All"] <- "AOB"
# rename columns in DROUGHT
colnames(aob.x1.dro)  <- c("spec", "env", "r", "p.value","p.adj")
view(aob.x1.dro)

################################################################################
# rarefied AOA ASV table of bulk soil DROUGHT

aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# select only data set until Jul 20th
aoa.BS.sub  <- aoa.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
aoa.BS.sub1 <- prune_taxa(taxa_sums(aoa.BS.sub)>0, aoa.BS.sub)
aoa.BS.sub1 # 445 taxa, 96 samples
# subset AOA DROUGHT
aoa.BS.dro_seq <- subset_samples(aoa.BS.sub1, Irrigation=="Rainout")
aoa.BS.dro_seq1 <- prune_taxa(taxa_sums(aoa.BS.dro_seq)>0, aoa.BS.dro_seq)
sort(rowSums(otu_table(aoa.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aoa.BS.dro_df <- as.data.frame(otu_table(aoa.BS.dro_seq1))
dim(aoa.BS.dro_df) # 349 otus, 48 samples nothing is missing
view(aoa.BS.dro_df)
# taxonomy table AOA DROUGHT
aoa.BS.dro_tax <- as.data.frame(phyloseq::tax_table(aoa.BS.dro_seq1))
dim(aoa.BS.dro_tax) # 349 5
aoa.BS.dro_tax$Group <- 'AOA'#adding a new column with AOB as value
view(aoa.BS.dro_tax)
# environmental data
view(aoa.env.Jul)
# environmental data CONTROL
aoa.BS.dro_env <- aoa.env.Jul %>% filter(Irrigation == "Rainout")
colnames(aoa.BS.dro_env)
view(aoa.BS.dro_env)
dim(aoa.BS.dro_env) # 48, 52
#extract alpha diversity DROUGHT
aoa.BS.dro_alpha <- aoa.BS.dro_env[,45:46]
head(aoa.BS.dro_alpha)
aoa.BS.dro_alpha <- rownames_to_column(aoa.BS.dro_alpha, var = "SampleID")
dim(aoa.BS.dro_alpha) # 48 3
# adding qPCR data DROUGHT
view(qpcr.BS.dro_logdws)
aoa.BS.dro_logdws <- qpcr.BS.dro_logdws
dim(aoa.BS.dro_logdws) # 48 4
view(aoa.BS.dro_logdws)
# combine ONLY DROUGHT
aoa.dro_env.ed <- cbind(aoa.BS.dro_env,aoa.BS.dro_logdws)
head(aoa.dro_env.ed)
dim(aoa.dro_env.ed) # 48 56
str(aoa.dro_env.ed)
# rename DROUGHT environment
aoa.dro_env.ed <- aoa.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names
aoa.dro_env.ed$UniqueID <- factor(aoa.dro_env.ed$UniqueID)
aoa.dro_env.ed$PlotID <- factor(aoa.dro_env.ed$PlotID)
aoa.dro_env.ed$Block <- factor(aoa.dro_env.ed$Block)
aoa.dro_env.ed$Irrigation <- factor(aoa.dro_env.ed$Irrigation)
aoa.dro_env.ed$Treatment <- factor(aoa.dro_env.ed$Treatment)
aoa.dro_env.ed$Type <- factor(aoa.dro_env.ed$Type)
aoa.dro_env.ed$rep <- factor(aoa.dro_env.ed$rep)
aoa.dro_env.ed$rep2 <- factor(aoa.dro_env.ed$rep2)
aoa.dro_env.ed$Date <- factor(aoa.dro_env.ed$Date)
aoa.dro_env.ed$x <- factor(aoa.dro_env.ed$x)
aoa.dro_env.ed$var3 <- factor(aoa.dro_env.ed$var3)
aoa.dro_env.ed$var2 <- factor(aoa.dro_env.ed$var2)
# adding alpha data of AOB of DROUGHT
dim(aob.BS.dro_alpha) # 48 3 complete
head(aob.BS.dro_alpha)
head(aoa.dro_env.ed)
aoa.dro_env.ed <- rownames_to_column(aoa.dro_env.ed, var = "SampleID")
aoa.dro_env.ed2 <- aoa.dro_env.ed %>%
 left_join(aob.BS.dro_alpha, by="SampleID") 
colnames(aoa.dro_env.ed2)
dim(aoa.dro_env.ed2) # 48 59
view(aoa.dro_env.ed2) 
# adding alpha data of COMA of DROUGHT
dim(com.BS.dro_alpha) # 46 3
aoa.dro_env.ed3 <- aoa.dro_env.ed2 %>%
 left_join(com.BS.dro_alpha, by="SampleID") 
dim(aoa.dro_env.ed3) # 48 61 (it has missing values of comammox alpha)
view(aoa.dro_env.ed3)
# relocate alpha data of AOA in DROUGHT
colnames(aoa.dro_env.ed3)
aoa.dro_env.ed4 <- aoa.dro_env.ed3 %>% relocate("AOA_richness", .after="Coma_B_abundance")%>% 
                                       relocate("AOA_Shannon", .after="AOA_richness")%>%
                                       relocate("N2O", .after="GWC")

head(aoa.dro_env.ed4)
view(aoa.dro_env.ed4)
dim(aoa.dro_env.ed4) # 48 61, with S26 & S52 missing Comammox alpha
aoa.dro_env.ed4 <- column_to_rownames(aoa.dro_env.ed4, var = "SampleID")
# create  a microtable for DROUGHT
aoa.dro_microdata <- microtable$new(sample_table = aoa.dro_env.ed4, otu_table = aoa.BS.dro_df, tax_table = aoa.BS.dro_tax)
# calculate beta diversity for DROUGHT
aoa.dro_microdata$tidy_dataset()
set.seed(13)
aoa.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test for DROUGHT
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
aoa.t1.dro <- trans_env$new(dataset = aoa.dro_microdata, env_cols = c(17,23:24,27,31:33,48:60),complete_na = TRUE) 
set.seed(13)
aoa.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = T,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
aoa.t1.dro$res_mantel
# extract a part of the results DROUGHT
aoa.x1.dro <- data.frame(aoa.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aoa.x1.dro[aoa.x1.dro=="All"] <- "AOA"
# rename columns DROUGHT
colnames(aoa.x1.dro) <- c("spec", "env", "r", "p.value","p.adj")
view(aoa.x1.dro)

################################################################################
# rarefied COMAMMOX ASV table of bulk soil DROUGHT

com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
dim(com.asv.tab) # 497, 118 samples (missing 2 values: SS26 & S52)
# select only data set until Jul 20th
com.BS.sub  <- com.physeq_bulk1 %>%
  subset_samples(Date %in% c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th"))
com.BS.sub1 <- prune_taxa(taxa_sums(com.BS.sub)>0, com.BS.sub)
com.BS.sub1 # 458 taxa, 94 samples
# subset COM DROUGHT
com.BS.dro_seq <- subset_samples(com.BS.sub1, Irrigation=="Rainout")
com.BS.dro_seq1 <- prune_taxa(taxa_sums(com.BS.dro_seq)>0, com.BS.dro_seq)
sort(rowSums(otu_table(com.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
com.BS.dro_df <- as.data.frame(otu_table(com.BS.dro_seq1))
dim(com.BS.dro_df) # 357 otus, 46 samples (missing because S26 & S52 are removed)
view(com.BS.dro_df)
# taxonomy table COM DROUGHT
com.BS.dro_tax <- as.data.frame(phyloseq::tax_table(com.BS.dro_seq1))
dim(com.BS.dro_tax) # 357 4
com.BS.dro_tax$Group <- 'COM'#adding a new column with AOB as value
view(com.BS.dro_tax)
# environmental data
dim(com.env.Jul) # 94 52
# environmental data DROUGHT
com.BS.dro_env <- com.env.Jul %>% filter(Irrigation == "Rainout")
colnames(com.BS.dro_env)
view(com.BS.dro_env)
dim(com.BS.dro_env) # 46 52 (S26 and S52 were removed from the metadata to match the otu table)
#extract alpha diversity DROUGHT
com.BS.dro_alpha <- com.BS.dro_env[,45:46]
head(com.BS.dro_alpha)
com.BS.dro_alpha <- rownames_to_column(com.BS.dro_alpha, var = "SampleID")
dim(com.BS.dro_alpha) # 46 3
# adding qPCR data DROUGHT
view(qpcr.BS.dro_logdws)
dim(qpcr.BS.dro_logdws) # 48 6 still contain S26 and S52
qpcr.BS.dro_logdws.x <- rownames_to_column(qpcr.BS.dro_logdws, var = "SampleID")
com.BS.dro_logdws <- qpcr.BS.dro_logdws.x %>% filter(SampleID != "S26",
                                                  SampleID !="S52")# filter out S11 from the qPCR data
com.BS.dro_logdws <- column_to_rownames(com.BS.dro_logdws, var = "SampleID")
dim(com.BS.dro_logdws) # 46 4
view(com.BS.dro_logdws)
# combine ONLY DROUGHT
com.dro_env.ed <- cbind(com.BS.dro_env,com.BS.dro_logdws)
head(com.dro_env.ed)
dim(com.dro_env.ed) # 46 56
str(com.dro_env.ed)
# rename DROUGHT environment
com.dro_env.ed <- com.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g",
         "N2O"="mean.N2Oflux") # change column names
com.dro_env.ed$UniqueID <- factor(com.dro_env.ed$UniqueID)
com.dro_env.ed$PlotID <- factor(com.dro_env.ed$PlotID)
com.dro_env.ed$Block <- factor(com.dro_env.ed$Block)
com.dro_env.ed$Irrigation <- factor(com.dro_env.ed$Irrigation)
com.dro_env.ed$Treatment <- factor(com.dro_env.ed$Treatment)
com.dro_env.ed$Type <- factor(com.dro_env.ed$Type)
com.dro_env.ed$rep <- factor(com.dro_env.ed$rep)
com.dro_env.ed$rep2 <- factor(com.dro_env.ed$rep2)
com.dro_env.ed$Date <- factor(com.dro_env.ed$Date)
com.dro_env.ed$x <- factor(com.dro_env.ed$x)
com.dro_env.ed$var3 <- factor(com.dro_env.ed$var3)
com.dro_env.ed$var2 <- factor(com.dro_env.ed$var2)
str(com.dro_env.ed)
# adding alpha data of AOA DROUGHT
head(aoa.BS.dro_alpha)
dim(aoa.BS.dro_alpha) # 48 2 complete
dim(com.dro_env.ed) # 46 56
aoa.BS.dro_alpha <- rownames_to_column(aoa.BS.dro_alpha, var = "SampleID")
com.dro_env.ed <- rownames_to_column(com.dro_env.ed, var = "SampleID")
com.dro_env.ed2 <- com.dro_env.ed %>%
 left_join(aoa.BS.dro_alpha, by="SampleID") 
dim(com.dro_env.ed2) #46 59 #the NA are automatically removed to match the comammox meta data
view(com.dro_env.ed2)
# adding alpha data of AOB DROUGHT
dim(aob.BS.dro_alpha) # 48 3 
com.dro_env.ed3 <- com.dro_env.ed2 %>%
 left_join(aob.BS.dro_alpha, by="SampleID") 
dim(com.dro_env.ed3) # 46 61, #the NA are automatically removed to match the comammox meta data
colnames(com.dro_env.ed3)
# relocate alpha data of COMA DROUGHT
com.dro_env.ed4 <- com.dro_env.ed3 %>% relocate("COM_richness", .after="AOB_Shannon")%>%
                                       relocate("COM_Shannon", .after="COM_richness")%>%
                                       relocate("N2O", .after="GWC")
head(com.dro_env.ed4)
dim(com.dro_env.ed4) # 58 56
com.dro_env.ed4 <- column_to_rownames(com.dro_env.ed4, var = "SampleID")
# create  a micro table DROUGHT
com.dro_microdata <- microtable$new(sample_table = com.dro_env.ed4, otu_table = com.BS.dro_df, tax_table = com.BS.dro_tax)
# calculate beta diversity DROUGHT
com.dro_microdata$tidy_dataset()
set.seed(13)
com.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test DROUGHT
set.seed(13)
com.t1.dro <- trans_env$new(dataset = com.dro_microdata, env_cols = c(17,23:24,27,31:33,48:60))
set.seed(13)
com.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = T,
                  permutations=999,
                  add_matrix = NULL,
                  method = "pearson",
                  p_adjust_method = "fdr")
com.t1.dro$res_mantel
# extract a part of the results DROUGHT
com.x1.dro <- data.frame(com.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
com.x1.dro[com.x1.dro=="All"] <- "comammox"
# rename columns DROUGHT
colnames(com.x1.dro)  <- c("spec", "env", "r", "p.value","p.adj")
view(com.x1.dro)

###############################################################################
# generate interval data DROUGHT
aob.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
aoa.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
com.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# combine three tables DROUGHT
plot_table.DRO <- rbind(aob.x1.dro, aoa.x1.dro,com.x1.dro)
plot_table.DRO <- plot_table.DRO %>% 
  mutate(spec= replace(spec, spec == "AOB", "AOB β diversity")) %>%
  mutate(spec= replace(spec, spec =="AOA","AOA β diversity")) %>%
  mutate(spec= replace(spec, spec=="comammox", "Comammox β diversity")) %>%# change column names
  mutate(env= replace(env, env =="Total_C","Total C")) %>%
  mutate(env= replace(env, env =="Total_N","Total N")) %>%
  mutate(env= replace(env, env =="AOA_abundance","AOA abundance")) %>%
  mutate(env= replace(env, env =="AOB_abundance","AOB abundance")) %>%
  mutate(env= replace(env, env =="Coma_A_abundance","COM-A abundance")) %>%
  mutate(env= replace(env, env =="Coma_B_abundance","COM-B abundance")) %>%
  mutate(env= replace(env, env =="AOA_richness","AOA richness")) %>%
  mutate(env= replace(env, env =="AOA_Shannon","AOA Shannon")) %>%
  mutate(env= replace(env, env =="AOB_richness","AOB richness")) %>%
  mutate(env= replace(env, env =="AOB_Shannon","AOB Shannon")) %>%
  mutate(env= replace(env, env =="COM_richness","COM richness")) %>%
  mutate(env= replace(env, env =="COM_Shannon","COM Shannon")) 
view(plot_table.DRO)
# plotting
#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
# extract env data from one of the group
aoa.env.dro <- aoa.t1.dro$data_env # the complete one
dim(aoa.env.dro) # 48 20
aob.env.dro <- aob.t1.dro$data_env 
view(aoa.env.dro)
view(aob.env.dro)
all.env.dro <- aoa.env.dro
colnames(all.env.dro)
str(all.env.dro)
all.env.dro.ed <- all.env.dro %>% 
  dplyr::rename("AOA abundance" = "AOA_abundance",
         "AOB abundance" = "AOB_abundance",
         "COM-A abundance" = "Coma_A_abundance",
         "COM-B abundance" = "Coma_B_abundance",
         "Total C" = "Total_C",
         "Total N" = "Total_N",
         "AOA richness" = "AOA_richness",
         "AOA Shannon" = "AOA_Shannon",
         "AOB richness" = "AOB_richness",
         "AOB Shannon" = "AOB_Shannon",
         "COM richness" = "COM_richness",
         "COM Shannon" = "COM_Shannon") 
str(all.env.dro.ed)
head(all.env.dro.ed)

# plot DROUGHT
mantelplot.drought <- quickcor(all.env.dro.ed, type = "upper", cor.test = TRUE,method="pearson",exact=FALSE) + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.DRO, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 3) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank())
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.drought

setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("mantel.drought_withN2O.tiff",
       mantelplot.drought, device = "tiff",
       width = 10, height = 10, 
       units= "in", dpi = 600)

# try with NA of COMAMMOX alpha
all.env.dro.ed2 <- rownames_to_column(all.env.dro.ed, var = "SampleID")
all.env.dro.ed3 <- all.env.dro.ed2 %>% 
  mutate(`COM richness`=ifelse(SampleID=="S26",NA,`COM richness`),
         `COM Shannon`=ifelse(SampleID=="S26",NA,`COM Shannon`),
         `COM richness`=ifelse(SampleID=="S52",NA,`COM richness`),
         `COM Shannon`=ifelse(SampleID=="S52",NA,`COM Shannon`))
view(all.env.dro.ed3)
all.env.dro.ed3 <- column_to_rownames(all.env.dro.ed3, var = "SampleID")
# plot
#install.packages("latex2exp")
library(latex2exp)
mantelplot.drought.noCOMalph <- quickcor(all.env.dro.ed3, type = "upper", 
                                         cor.test = TRUE,
                                         #cluster = TRUE,
                                         method="pearson",
                                         exact=F,
                                         show.diag = F,
                                         use="pairwise.complete.obs") + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, 
            color = "black", 
            size=2.7)+
            #digits = 2,
            #nudge_x = 0, 
            #nudge_y = 0,) +
  add_link(plot_table.DRO, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 10) +
  guides(fill = guide_colorbar(title = "Spearman's r",order = 3),
         size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1))+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=15),
        plot.margin = margin(0, 0, 0, 0, "cm"))
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.drought.noCOMalph

ggsave("mantelplot.drought.noCOMalph2_Pearson.tiff",
       mantelplot.drought.noCOMalph, device = "tiff",
       width = 12, height = 8.1, 
       units= "in", dpi = 600)
