##################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST COMAMMOX
#################################################################################################

# GROUP & SEPARATE PHYLOSEQ OBJECT BY TYPE, DATE AND TREATMENT:

##### 1. BULK SOIL - RAREFIED #####
# rarefied ASV table # 632 taxa
com.rare.min.physeq
# separate the bulk soil
com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS")
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1 # 497 taxa
################################################################################
# Subset the data set per irrigation-treatment-date

# Date: 04-28-2022
# 1. MINERAL
comM04seq<- subset_samples(com.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
comM04seq1 <- prune_taxa(taxa_sums(comM04seq)>0, comM04seq)
sort(rowSums(otu_table(comM04seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
comD04seq<- subset_samples(com.physeq_bulk1, Date=="04-28-22" & Treatment=="D")
comD04seq1 <- prune_taxa(taxa_sums(comD04seq)>0, comD04seq)
# 3. CONVENTIONAL
comK04seq<- subset_samples(com.physeq_bulk1, Date=="04-28-22" & Treatment=="K")
comK04seq1 <- prune_taxa(taxa_sums(comK04seq)>0, comK04seq)

# Date: 06-01-2022
# 1. MINERAL
comM06seq<- subset_samples(com.physeq_bulk1, Date=="06-01-22" & Treatment=="M")
comM06seq1 <- prune_taxa(taxa_sums(comM06seq)>0, comM06seq)
# 2. BIODYNAMIC
comD06seq<- subset_samples(com.physeq_bulk1, Date=="06-01-22" & Treatment=="D")
comD06seq1 <- prune_taxa(taxa_sums(comD06seq)>0, comD06seq)
# 3. CONVENTIONAL
comK06seq<- subset_samples(com.physeq_bulk1, Date=="06-01-22" & Treatment=="K")
comK06seq1 <- prune_taxa(taxa_sums(comK06seq)>0, comK06seq)

# Date: 07-05-2022
# 1. MINERAL
comM0705seq<- subset_samples(com.physeq_bulk1, Date=="07-05-22" & Treatment=="M")
comM0705seq1 <- prune_taxa(taxa_sums(comM0705seq)>0, comM0705seq)
# 2. BIODYNAMIC
comD0705seq<- subset_samples(com.physeq_bulk1, Date=="07-05-22" & Treatment=="D")
comD0705seq1 <- prune_taxa(taxa_sums(comD0705seq)>0, comD0705seq)
# 3. CONVENTIONAL
comK0705seq<- subset_samples(com.physeq_bulk1, Date=="07-05-22" & Treatment=="K")
comK0705seq1 <- prune_taxa(taxa_sums(comK0705seq)>0, comK0705seq)

# Date: 07-20-2022
# 1. MINERAL
comM0720seq<- subset_samples(com.physeq_bulk1, Date=="07-20-22" & Treatment=="M")
comM0720seq1 <- prune_taxa(taxa_sums(comM0720seq)>0, comM0720seq)
# 2. BIODYNAMIC
comD0720seq<- subset_samples(com.physeq_bulk1, Date=="07-20-22" & Treatment=="D")
comD0720seq1 <- prune_taxa(taxa_sums(comD0720seq)>0, comD0720seq)
# 3. CONVENTIONAL
comK0720seq<- subset_samples(com.physeq_bulk1, Date=="07-20-22" & Treatment=="K")
comK0720seq1 <- prune_taxa(taxa_sums(comK0720seq)>0, comK0720seq)

# Date: 09-13-2022
# 1. MINERAL
comM09seq<- subset_samples(com.physeq_bulk1, Date=="09-13-22" & Treatment=="M")
comM09seq1 <- prune_taxa(taxa_sums(comM09seq)>0, comM09seq)
# 2. BIODYNAMIC
comD09seq<- subset_samples(com.physeq_bulk1, Date=="09-13-22" & Treatment=="D")
comD09seq1 <- prune_taxa(taxa_sums(comD09seq)>0, comD09seq)
# 3. CONVENTIONAL
comK09seq<- subset_samples(com.physeq_bulk1, Date=="09-13-22" & Treatment=="K")
comK09seq1 <- prune_taxa(taxa_sums(comK09seq)>0, comK09seq)

################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- comK0720seq1
physeq.subset #
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Irrigation)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Irrigation)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Irrigation)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Irrigation == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Irrigation == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s)) # using unfiltered data
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson', offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_BulkSoil_rare_prev80/')
#write.csv(glmT3s.pairwise.global, file = "COM_K09_130923.csv")
#com.K09.fil <- as.data.frame(otu_table(physeq.subset.75))
#write.csv(com.K09.fil, file = "COM_K09.tab_130923.csv")

################################################################################
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_BulkSoil_rare_prev80/')
COM_M04 <- read.csv("COM_M04_130923.csv")[,-1]
COM_M04$contrast <- paste("M_042822", COM_M04$contrast, sep="_")
COM_D04 <- read.csv("COM_D04_130923.csv")[,-1]
COM_D04$contrast <- paste("D_042822", COM_D04$contrast, sep="_")
COM_K04 <- read.csv("COM_K04_130923.csv")[,-1]
COM_K04$contrast <- paste("K_042822", COM_K04$contrast, sep="_")

COM_M06 <- read.csv("COM_M06_130923.csv")[,-1]
COM_M06$contrast <- paste("M_060122", COM_M06$contrast, sep="_")
COM_D06 <- read.csv("COM_D06_130923.csv")[,-1]
COM_D06$contrast <- paste("D_060122", COM_D06$contrast, sep="_")
COM_K06 <- read.csv("COM_K06_130923.csv")[,-1]
COM_K06$contrast <- paste("K_060122", COM_K06$contrast, sep="_")

COM_M0705 <- read.csv("COM_M0705_130923.csv")[,-1]
COM_M0705$contrast <- paste("M_070522", COM_M0705$contrast, sep="_")
COM_D0705 <- read.csv("COM_D0705_130923.csv")[,-1]
COM_D0705$contrast <- paste("D_070522", COM_D0705$contrast, sep="_")
COM_K0705 <- read.csv("COM_K0705_130923.csv")[,-1]
COM_K0705$contrast <- paste("K_070522", COM_K0705$contrast, sep="_")

COM_M0720 <- read.csv("COM_M0720_130923.csv")[,-1]
COM_M0720$contrast <- paste("M_072022", COM_M0720$contrast, sep="_")
COM_D0720 <- read.csv("COM_D0720_130923.csv")[,-1]
COM_D0720$contrast <- paste("D_072022", COM_D0720$contrast, sep="_")
COM_K0720 <- read.csv("COM_K0720_130923.csv")[,-1]
COM_K0720$contrast <- paste("K_072022", COM_K0720$contrast, sep="_")

COM_M09 <- read.csv("COM_M09_130923.csv")[,-1]
COM_M09$contrast <- paste("M_091322", COM_M09$contrast, sep="_")
COM_D09 <- read.csv("COM_D09_130923.csv")[,-1]
COM_D09$contrast <- paste("D_091322", COM_D09$contrast, sep="_")
COM_K09 <- read.csv("COM_K09_130923.csv")[,-1]
COM_K09$contrast <- paste("K_091322", COM_K09$contrast, sep="_")

glmT3s.pairwise.global.ALL <- rbind(COM_M04, COM_D04, COM_K04, COM_M06, COM_D06, COM_K06,
                                    COM_M0705, COM_D0705, COM_K0705, COM_M0720, COM_D0720, COM_K0720,
                                    COM_M09, COM_D09, COM_K09)

## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.ALL$p.value <= 0.05)
table(glmT3s.pairwise.global.ALL$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global.ALL$OTU[glmT3s.pairwise.global.ALL$p.adjust <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global.ALL[glmT3s.pairwise.global.ALL$p.adjust <=0.05,]

length(tmp_otu3s)
tmp_otu3s

# cast pvalues
contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global.ALL[,c(10,1,2)]
# numeric variable needs to be named "value" 
colnames(contrasts.glm.CBFP.T3s) <- c("value", "OTU_names", "contrast")
#contrasts.glm.CBFP.T3s <- subset(contrasts.glm.CBFP.T3s, (contrasts.glm.CBFP.T3s$OTU_names %in% BFPOTUs.T3snet.sig))
head(contrasts.glm.CBFP.T3s)
str(contrasts.glm.CBFP.T3s)
contrasts.glm.CBFP.T3s <- data.frame(cast(contrasts.glm.CBFP.T3s, contrast ~ OTU_names, value="value"))
str(contrasts.glm.CBFP.T3s)
rownames(contrasts.glm.CBFP.T3s) <- contrasts.glm.CBFP.T3s$contrast
contrasts.glm.CBFP.T3s$contrast <- NULL


# keep OTUs with at least one contrast <0.05 
contrasts.glm.CBFP.T3s.sub <- contrasts.glm.CBFP.T3s[,colSums(contrasts.glm.CBFP.T3s<0.05, na.rm=TRUE) >= 1]
dim(contrasts.glm.CBFP.T3s.sub)
head(contrasts.glm.CBFP.T3s.sub)
str(contrasts.glm.CBFP.T3s.sub)

ctrst.glm.CBFP.T3s.sub <- data.frame(t(contrasts.glm.CBFP.T3s.sub))

# replace pvalues to 0 if non significant, or 1 if significant
#ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub ==NA] <- 0
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >0.05] <- 2
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub <0.05] <- 1
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >1] <- 0
ctrst.glm.CBFP.T3s.sub[is.na(ctrst.glm.CBFP.T3s.sub)] <- 0
head(ctrst.glm.CBFP.T3s.sub)

# Calculate the OTU avg per treatment
# CHECK THE OBJECT
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
meanotus<-phyloseq_average(com.physeq_bulk1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
meanotus<-as.data.frame(otu_table(meanotus));meanotus

# same order for both meanotus and tmp_otu3s
meanotus<-meanotus[tmp_otu3s,]
#meanotus<-meanotus[c("ASV_10", "ASV_12", "ASV_28", "ASV_33", "ASV_35", "ASV_44", "ASV_52", "ASV_54", "ASV_59", "ASV_64", "ASV_71"),]
#meanotus<-meanotus[,c("ASV_10", "ASV_12", "ASV_28", "ASV_33", "ASV_35", "ASV_44", "ASV_52", "ASV_54", "ASV_59", "ASV_64", "ASV_71"),]

# Calculate log2fold ratios for all OTUs in the filtered table

meanotus$RR_M_042822 <- log2(meanotus$RMBS1 / meanotus$CMBS1)
meanotus$RR_M_060122 <- log2(meanotus$RMBS2 / meanotus$CMBS2)
meanotus$RR_M_070522 <- log2(meanotus$RMBS3 / meanotus$CMBS3)
meanotus$RR_M_072022 <- log2(meanotus$RMBS4 / meanotus$CMBS4)
meanotus$RR_M_091322 <- log2(meanotus$RMBS5 / meanotus$CMBS5)

meanotus$RR_D_042822 <- log2(meanotus$RDBS1 / meanotus$CDBS1)
meanotus$RR_D_060122 <- log2(meanotus$RDBS2 / meanotus$CDBS2)
meanotus$RR_D_070522 <- log2(meanotus$RDBS3 / meanotus$CDBS3)
meanotus$RR_D_072022 <- log2(meanotus$RDBS4 / meanotus$CDBS4)
meanotus$RR_D_091322 <- log2(meanotus$RDBS5 / meanotus$CDBS5)

meanotus$RR_K_042822 <- log2(meanotus$RKBS1 / meanotus$CKBS1)
meanotus$RR_K_060122 <- log2(meanotus$RKBS2 / meanotus$CKBS2)
meanotus$RR_K_070522 <- log2(meanotus$RKBS3 / meanotus$CKBS3)
meanotus$RR_K_072022 <- log2(meanotus$RKBS4 / meanotus$CKBS4)
meanotus$RR_K_091322 <- log2(meanotus$RKBS5 / meanotus$CKBS5)


head(meanotus)
# keep only columns containing log2fold ratios (RRs)
meanotus<-meanotus[,c(31:45)]
head(meanotus)
meanotus[meanotus == "-Inf"] <- 0
meanotus[meanotus == "Inf"] <- 0
meanotus[meanotus == "NaN"] <- 0


ctrst.glm.CBFP.T3s.sub.ed <- 
  head(ctrst.glm.CBFP.T3s.sub)

# put the same column order in ctrst.glm.CBFP.T3s.sub and in meanotus
ctrst.glm.CBFP.T3s.sub<-ctrst.glm.CBFP.T3s.sub[,c(11,12,13,14,15,
                                                  1,2,3,4,5,
                                                  6,7,8,9,10)]
ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[rownames(meanotus), ]

# Multiply the matrices to get the RR when it is significant and 0 when it is not significant
rr<-meanotus*ctrst.glm.CBFP.T3s.sub

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_BulkSoil_rare_prev80/')
write.csv(rr, file = "COMA_RR_130923.csv")


# HeatMap

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_BulkSoil_rare_prev80/')
rr.com <- read.csv("COMA_RR_130923.csv", row.names = 1)
names(rr.com)=str_sub(names(rr.com),4)
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_Rhizo_rare_prev80/')
rr.com.rhizo <- read.csv("COMA_RR_Rhizo_140923.csv", row.names = 1)
names(rr.com.rhizo)=str_sub(names(rr.com.rhizo),4)

#install.packages("colorRamp2")
library(colorRamp2)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
com.bs.hm <- Heatmap(as.matrix(rr.com), cluster_columns = F, col= col_fun)
com.bs.hm
com.rz.hm <- Heatmap(as.matrix(rr.com.rhizo), cluster_columns = F, col= col_fun)
com.rz.hm


#################################################################################

##### 2. RHIZOSPHERE - RAREFIED #####

com.physeq_rh <- subset_samples(com.rare.min.physeq, Type=="RS")
com.physeq_rh # 632 taxa 72 samples
com.physeq_rh1 <- prune_taxa(taxa_sums(com.physeq_rh)>0, com.physeq_rh)
com.physeq_rh1 # 448 taxa

# Date: 04-28-2022
# 1. MINERAL
M04.rh.seq<- subset_samples(com.physeq_rh1, Date=="04-28-22" & Treatment=="M")
M04.rh.seq1 <- prune_taxa(taxa_sums(M04.rh.seq)>0, M04.rh.seq)
sort(rowSums(otu_table(M04.rh.seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04.rh.seq<- subset_samples(com.physeq_rh1, Date=="04-28-22" & Treatment=="D")
D04.rh.seq1 <- prune_taxa(taxa_sums(D04.rh.seq)>0, D04.rh.seq)
# 3. CONVENTIONAL
K04.rh.seq<- subset_samples(com.physeq_rh1, Date=="04-28-22" & Treatment=="K")
K04.rh.seq1 <- prune_taxa(taxa_sums(K04.rh.seq)>0, K04.rh.seq)

# Date: 06-01-2022
# 1. MINERAL
M06.rh.seq<- subset_samples(com.physeq_rh1, Date=="06-01-22" & Treatment=="M")
M06.rh.seq1 <- prune_taxa(taxa_sums(M06.rh.seq)>0, M06.rh.seq)
# 2. BIODYNAMIC
D06.rh.seq<- subset_samples(com.physeq_rh1, Date=="06-01-22" & Treatment=="D")
D06.rh.seq1 <- prune_taxa(taxa_sums(D06.rh.seq)>0, D06.rh.seq)
# 3. CONVENTIONAL
K06.rh.seq<- subset_samples(com.physeq_rh1, Date=="06-01-22" & Treatment=="K")
K06.rh.seq1 <- prune_taxa(taxa_sums(K06.rh.seq)>0, K06.rh.seq)

# Date: 07-05-2022
# 1. MINERAL
M0705.rh.seq<- subset_samples(com.physeq_rh1, Date=="07-05-22" & Treatment=="M")
M0705.rh.seq1 <- prune_taxa(taxa_sums(M0705.rh.seq)>0, M0705.rh.seq)
# 2. BIODYNAMIC
D0705.rh.seq<- subset_samples(com.physeq_rh1, Date=="07-05-22" & Treatment=="D")
D0705.rh.seq1 <- prune_taxa(taxa_sums(D0705.rh.seq)>0, D0705.rh.seq)
# 3. CONVENTIONAL
K0705.rh.seq<- subset_samples(com.physeq_rh1, Date=="07-05-22" & Treatment=="K")
K0705.rh.seq1 <- prune_taxa(taxa_sums(K0705.rh.seq)>0, K0705.rh.seq)
################################################################################
###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- K0705.rh.seq1
physeq.subset 
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Irrigation)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Irrigation)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Irrigation)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Irrigation == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Irrigation == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 80,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75  # 32 taxa

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s)) # using unfiltered data
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson', offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_Rhizo_rare_prev80/')
write.csv(glmT3s.pairwise.global, file = "COM_K0705.rh_140923.csv")
com.K0705.rh.fil <- as.data.frame(otu_table(physeq.subset.75))
write.csv(com.K0705.rh.fil, file = "COM_K0705.rh.tab_140923.csv")

##############################################################################################################################
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_Rhizo_rare_prev80/')
COM_M04 <- read.csv("COM_M04.rh_140923.csv")[,-1]
COM_M04$contrast <- paste("M_042822", COM_M04$contrast, sep="_")
COM_D04 <- read.csv("COM_D04.rh_140923.csv")[,-1]
COM_D04$contrast <- paste("D_042822", COM_D04$contrast, sep="_")
COM_K04 <- read.csv("COM_K04.rh_140923.csv")[,-1]
COM_K04$contrast <- paste("K_042822", COM_K04$contrast, sep="_")

COM_M06 <- read.csv("COM_M06.rh_140923.csv")[,-1]
COM_M06$contrast <- paste("M_060122", COM_M06$contrast, sep="_")
COM_D06 <- read.csv("COM_D06.rh_140923.csv")[,-1]
COM_D06$contrast <- paste("D_060122", COM_D06$contrast, sep="_")
COM_K06 <- read.csv("COM_K06.rh_140923.csv")[,-1]
COM_K06$contrast <- paste("K_060122", COM_K06$contrast, sep="_")

COM_M0705 <- read.csv("COM_M0705.rh_140923.csv")[,-1]
COM_M0705$contrast <- paste("M_070522", COM_M0705$contrast, sep="_")
COM_D0705 <- read.csv("COM_D0705.rh_140923.csv")[,-1]
COM_D0705$contrast <- paste("D_070522", COM_D0705$contrast, sep="_")
COM_K0705 <- read.csv("COM_K0705.rh_140923.csv")[,-1]
COM_K0705$contrast <- paste("K_070522", COM_K0705$contrast, sep="_")

glmT3s.pairwise.global.ALL <- rbind(COM_M04, COM_D04, COM_K04, COM_M06, COM_D06, COM_K06,
                                    COM_M0705, COM_D0705, COM_K0705)

## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.ALL$p.value <= 0.05)
table(glmT3s.pairwise.global.ALL$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global.ALL$OTU[glmT3s.pairwise.global.ALL$p.adjust <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global.ALL[glmT3s.pairwise.global.ALL$p.adjust <=0.05,]

length(tmp_otu3s)
tmp_otu3s

# cast pvalues
contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global.ALL[,c(10,1,2)]
# numeric variable needs to be named "value" 
colnames(contrasts.glm.CBFP.T3s) <- c("value", "OTU_names", "contrast")
#contrasts.glm.CBFP.T3s <- subset(contrasts.glm.CBFP.T3s, (contrasts.glm.CBFP.T3s$OTU_names %in% BFPOTUs.T3snet.sig))
head(contrasts.glm.CBFP.T3s)
str(contrasts.glm.CBFP.T3s)
contrasts.glm.CBFP.T3s <- data.frame(cast(contrasts.glm.CBFP.T3s, contrast ~ OTU_names, value="value"))
str(contrasts.glm.CBFP.T3s)
rownames(contrasts.glm.CBFP.T3s) <- contrasts.glm.CBFP.T3s$contrast
contrasts.glm.CBFP.T3s$contrast <- NULL


# keep OTUs with at least one contrast <0.05 
contrasts.glm.CBFP.T3s.sub <- contrasts.glm.CBFP.T3s[,colSums(contrasts.glm.CBFP.T3s<0.05, na.rm=TRUE) >= 1]
dim(contrasts.glm.CBFP.T3s.sub)
head(contrasts.glm.CBFP.T3s.sub)
str(contrasts.glm.CBFP.T3s.sub)

ctrst.glm.CBFP.T3s.sub <- data.frame(t(contrasts.glm.CBFP.T3s.sub))

# replace pvalues to 0 if non significant, or 1 if significant
#ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub ==NA] <- 0
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >0.05] <- 2
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub <0.05] <- 1
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >1] <- 0
ctrst.glm.CBFP.T3s.sub[is.na(ctrst.glm.CBFP.T3s.sub)] <- 0
head(ctrst.glm.CBFP.T3s.sub)

# Calculate the OTU avg per treatment
# CHECK THE OBJECT
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
meanotus<-phyloseq_average(com.physeq_rh1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
meanotus<-as.data.frame(otu_table(meanotus));meanotus

# same order for both meanotus and tmp_otu3s
meanotus<-meanotus[tmp_otu3s,]
#meanotus<-meanotus[c("ASV_10", "ASV_12", "ASV_28", "ASV_33", "ASV_35", "ASV_44", "ASV_52", "ASV_54", "ASV_59", "ASV_64", "ASV_71"),]
#meanotus<-meanotus[,c("ASV_10", "ASV_12", "ASV_28", "ASV_33", "ASV_35", "ASV_44", "ASV_52", "ASV_54", "ASV_59", "ASV_64", "ASV_71"),]

# Calculate log2fold ratios for all OTUs in the filtered table

meanotus$RR_M_042822 <- log2(meanotus$RMRS1 / meanotus$CMRS1)
meanotus$RR_M_060122 <- log2(meanotus$RMRS2 / meanotus$CMRS2)
meanotus$RR_M_070522 <- log2(meanotus$RMRS3 / meanotus$CMRS3)

meanotus$RR_D_042822 <- log2(meanotus$RDRS1 / meanotus$CDRS1)
meanotus$RR_D_060122 <- log2(meanotus$RDRS2 / meanotus$CDRS2)
meanotus$RR_D_070522 <- log2(meanotus$RDRS3 / meanotus$CDRS3)

meanotus$RR_K_042822 <- log2(meanotus$RKRS1 / meanotus$CKRS1)
meanotus$RR_K_060122 <- log2(meanotus$RKRS2 / meanotus$CKRS2)
meanotus$RR_K_070522 <- log2(meanotus$RKRS3 / meanotus$CKRS3)

head(meanotus)
# keep only columns containing log2fold ratios (RRs)
meanotus<-meanotus[,c(19:27)]
head(meanotus)
meanotus[meanotus == "-Inf"] <- 0
meanotus[meanotus == "Inf"] <- 0
meanotus[meanotus == "NaN"] <- 0

# put the same column order in ctrst.glm.CBFP.T3s.sub and in meanotus
ctrst.glm.CBFP.T3s.sub<-ctrst.glm.CBFP.T3s.sub[,c(7,8,9,
                                                  1,2,3,
                                                  4,5,6)]
ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[rownames(meanotus), ]

# Multiply the matrices to get the RR when it is significant and 0 when it is not significant
rr<-meanotus*ctrst.glm.CBFP.T3s.sub

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/COM_Rhizo_rare_prev80/')
write.csv(rr, file = "COMA_RR_Rhizo_140923.csv")

################################################################################
##### 1. BULK SOIL - NOT RAREFIED #####

com.physeq # still contain rhizosphere!
com.bulk.rawseq <- subset_samples(com.physeq, Type=="BS")
com.bulk.rawseq1 <- prune_taxa(taxa_sums(com.bulk.rawseq)>0,com.bulk.rawseq)
com.bulk.rawseq1 # 531 taxa

# separate data by treatment and date

# Date: 04-28-2022
# 1. MINERAL
comM04rawseq<- subset_samples(com.bulk.rawseq1, Date=="04-28-22" & Treatment=="M")
comM04rawseq1 <- prune_taxa(taxa_sums(comM04rawseq)>0, comM04rawseq)
comM04_table <- as.data.frame(otu_table(comM04rawseq1))
# 2. BIODYNAMIC
comD04rawseq<- subset_samples(com.bulk.rawseq1, Date=="04-28-22" & Treatment=="D")
comD04rawseq1 <- prune_taxa(taxa_sums(comD04rawseq)>0, comD04rawseq)
comD04_table <- as.data.frame(otu_table(comD04rawseq1))
# 3. CONVENTIONAL
comK04rawseq<- subset_samples(com.bulk.rawseq1, Date=="04-28-22" & Treatment=="K")
comK04rawseq1 <- prune_taxa(taxa_sums(comK04rawseq)>0, comK04rawseq)
comK04_table <- as.data.frame(otu_table(comK04rawseq1))

# Date: 06-01-2022
# 1. MINERAL
comM06rawseq<- subset_samples(com.bulk.rawseq1, Date=="06-01-22" & Treatment=="M")
comM06rawseq1 <- prune_taxa(taxa_sums(comM06rawseq)>0, comM06rawseq)
comM06_table <- as.data.frame(otu_table(comM06rawseq1))
# 2. BIODYNAMIC
comD06rawseq<- subset_samples(com.bulk.rawseq1, Date=="06-01-22" & Treatment=="D")
comD06rawseq1 <- prune_taxa(taxa_sums(comD06rawseq)>0, comD06rawseq)
comD06_table <- as.data.frame(otu_table(comD06rawseq1))
# 3. CONVENTIONAL
comK06rawseq<- subset_samples(com.bulk.rawseq1, Date=="06-01-22" & Treatment=="K")
comK06rawseq1 <- prune_taxa(taxa_sums(comK06rawseq)>0, comK06rawseq)
comK06_table <- as.data.frame(otu_table(comK06rawseq1))

# Date: 07-05-2022
# 1. MINERAL
comM0705rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-05-22" & Treatment=="M")
comM0705rawseq1 <- prune_taxa(taxa_sums(comM0705rawseq)>0, comM0705rawseq)
# 2. BIODYNAMIC
comD0705rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-05-22" & Treatment=="D")
comD0705rawseq1 <- prune_taxa(taxa_sums(comD0705rawseq)>0, comD0705rawseq)
# 3. CONVENTIONAL
comK0705rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-05-22" & Treatment=="K")
comK0705rawseq1 <- prune_taxa(taxa_sums(comK0705rawseq)>0, comK0705rawseq)

# Date: 07-20-2022
# 1. MINERAL
comM0720rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-20-22" & Treatment=="M")
comM0720rawseq1 <- prune_taxa(taxa_sums(comM0720rawseq)>0, comM0720rawseq)
# 2. BIODYNAMIC
comD0720rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-20-22" & Treatment=="D")
comD0720rawseq1 <- prune_taxa(taxa_sums(comD0720rawseq)>0, comD0720rawseq)
# 3. CONVENTIONAL
comK0720rawseq<- subset_samples(com.bulk.rawseq1, Date=="07-20-22" & Treatment=="K")
comK0720rawseq1 <- prune_taxa(taxa_sums(comK0720rawseq)>0, comK0720rawseq)

# Date: 09-13-2022
# 1. MINERAL
comM09rawseq<- subset_samples(com.bulk.rawseq1, Date=="09-13-22" & Treatment=="M")
comM09rawseq1 <- prune_taxa(taxa_sums(comM09rawseq)>0, comM09rawseq)
# 2. BIODYNAMIC
comD09rawseq<- subset_samples(com.bulk.rawseq1, Date=="09-13-22" & Treatment=="D")
comD09rawseq1 <- prune_taxa(taxa_sums(comD09rawseq)>0, comD09rawseq)
# 3. CONVENTIONAL
comK09rawseq<- subset_samples(com.bulk.rawseq1, Date=="09-13-22" & Treatment=="K")
comK09rawseq1 <- prune_taxa(taxa_sums(comK09rawseq)>0, comK09rawseq)

################################################################################
################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- comM0720rawseq1
physeq.subset #
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Irrigation)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Irrigation)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Irrigation)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Irrigation == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Irrigation == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 80,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 46 taxa

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(comM0720rawseq1)) # using unfiltered data
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson', offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOA_BulkSoil_rare_prev80/')
write.csv(glmT3s.pairwise.global, file = "AOA_K09_130923.csv")
aoa.K09.fil <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aoa.K09.fil, file = "AOA_K09.tab_130923.csv")
################################################################################