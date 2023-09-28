##################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST
#################################################################################################

# GROUP & SEPARATE PHYLOSEQ OBJECT BY TYPE, DATE AND TREATMENT:

##### 1. BULK SOIL - RAREFIED #####

aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS")
aob.physeq_bulk
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1 # 937 taxa, 119 samples

################################################################################
### Subset the data set per irrigation-treatment-date (RAREFIED)

# Date: 04-28-2022
# 1. MINERAL
M04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
M04seq1 <- prune_taxa(taxa_sums(M04seq)>0, M04seq)
sort(rowSums(otu_table(M04seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="D")
D04seq1 <- prune_taxa(taxa_sums(D04seq)>0, D04seq)
# 3. CONVENTIONAL
K04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="K")
K04seq1 <- prune_taxa(taxa_sums(K04seq)>0, K04seq)

# Date: 06-01-2022
# 1. MINERAL
M06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="M")
M06seq1 <- prune_taxa(taxa_sums(M06seq)>0, M06seq)
# 2. BIODYNAMIC
D06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="D")
D06seq1 <- prune_taxa(taxa_sums(D06seq)>0, D06seq)
# 3. CONVENTIONAL
K06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="K")
K06seq1 <- prune_taxa(taxa_sums(K06seq)>0, K06seq)

# Date: 07-05-2022
# 1. MINERAL
M0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="M")
M0705seq1 <- prune_taxa(taxa_sums(M0705seq)>0, M0705seq)
# 2. BIODYNAMIC
D0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="D")
D0705seq1 <- prune_taxa(taxa_sums(D0705seq)>0, D0705seq)
# 3. CONVENTIONAL
K0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="K")
K0705seq1 <- prune_taxa(taxa_sums(K0705seq)>0, K0705seq)

# Date: 07-20-2022
# 1. MINERAL
M0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="M")
M0720seq1 <- prune_taxa(taxa_sums(M0720seq)>0, M0720seq)
# 2. BIODYNAMIC
D0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="D")
D0720seq1 <- prune_taxa(taxa_sums(D0720seq)>0, D0720seq)
# 3. CONVENTIONAL
K0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="K")
K0720seq1 <- prune_taxa(taxa_sums(K0720seq)>0, K0720seq)

# Date: 09-13-2022
# 1. MINERAL
M09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="M")
M09seq1 <- prune_taxa(taxa_sums(M09seq)>0, M09seq)
# 2. BIODYNAMIC
D09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="D")
D09seq1 <- prune_taxa(taxa_sums(D09seq)>0, D09seq)
# 3. CONVENTIONAL
K09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="K")
K09seq1 <- prune_taxa(taxa_sums(K09seq)>0, K09seq)


##### 2. RHIZOSPHERE - RAREFIED #####

aob.physeq_rh <- subset_samples(aob.rare.1282.seq, Type=="RS")
aob.physeq_rh # 1222 taxa 72 samples
aob.physeq_rh1 <- prune_taxa(taxa_sums(aob.physeq_rh)>0, aob.physeq_rh)
aob.physeq_rh1 # 831 taxa

# Date: 04-28-2022
# 1. MINERAL
M04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="M")
M04.rh.seq1 <- prune_taxa(taxa_sums(M04.rh.seq)>0, M04.rh.seq)
sort(rowSums(otu_table(M04.rh.seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="D")
D04.rh.seq1 <- prune_taxa(taxa_sums(D04.rh.seq)>0, D04.rh.seq)
# 3. CONVENTIONAL
K04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="K")
K04.rh.seq1 <- prune_taxa(taxa_sums(K04.rh.seq)>0, K04.rh.seq)

# Date: 06-01-2022
# 1. MINERAL
M06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="M")
M06.rh.seq1 <- prune_taxa(taxa_sums(M06.rh.seq)>0, M06.rh.seq)
# 2. BIODYNAMIC
D06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="D")
D06.rh.seq1 <- prune_taxa(taxa_sums(D06.rh.seq)>0, D06.rh.seq)
# 3. CONVENTIONAL
K06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="K")
K06.rh.seq1 <- prune_taxa(taxa_sums(K06.rh.seq)>0, K06.rh.seq)

# Date: 07-05-2022
# 1. MINERAL
M0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="M")
M0705.rh.seq1 <- prune_taxa(taxa_sums(M0705.rh.seq)>0, M0705.rh.seq)
# 2. BIODYNAMIC
D0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="D")
D0705.rh.seq1 <- prune_taxa(taxa_sums(D0705.rh.seq)>0, D0705.rh.seq)
# 3. CONVENTIONAL
K0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="K")
K0705.rh.seq1 <- prune_taxa(taxa_sums(K0705.rh.seq)>0, K0705.rh.seq)
################################################################################

###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- M04seq1
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
    glmT3s <- glmmTMB(y ~ a + (1|z) , family="poisson", offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
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

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
write.csv(glmT3s.pairwise.global, file = "AOB_K0705.rh_130923.csv")
aob.K0705.rh.fil <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob.K0705.rh.fil, file = "AOB_K0705.rh.tab_130923.csv")
##############################################################################################################################

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare/')
AOB_M04 <- read.csv("AOB_M04_070923.csv")[,-1]
AOB_M04$contrast <- paste("M_042822", AOB_M04$contrast, sep="_")
AOB_D04 <- read.csv("AOB_D04_070923.csv")[,-1]
AOB_D04$contrast <- paste("D_042822", AOB_D04$contrast, sep="_")
AOB_K04 <- read.csv("AOB_K04_070923.csv")[,-1]
AOB_K04$contrast <- paste("K_042822", AOB_K04$contrast, sep="_")

AOB_M06 <- read.csv("AOB_M06_070923.csv")[,-1]
AOB_M06$contrast <- paste("M_060122", AOB_M06$contrast, sep="_")
AOB_D06 <- read.csv("AOB_D06_070923.csv")[,-1]
AOB_D06$contrast <- paste("D_060122", AOB_D06$contrast, sep="_")
AOB_K06 <- read.csv("AOB_K06_070923.csv")[,-1]
AOB_K06$contrast <- paste("K_060122", AOB_K06$contrast, sep="_")

AOB_M0705 <- read.csv("AOB_M0705_070923.csv")[,-1]
AOB_M0705$contrast <- paste("M_070522", AOB_M0705$contrast, sep="_")
AOB_D0705 <- read.csv("AOB_D0705_070923.csv")[,-1]
AOB_D0705$contrast <- paste("D_070522", AOB_D0705$contrast, sep="_")
AOB_K0705 <- read.csv("AOB_K0705_070923.csv")[,-1]
AOB_K0705$contrast <- paste("K_070522", AOB_K0705$contrast, sep="_")

AOB_M0720 <- read.csv("AOB_M0720_070923.csv")[,-1]
AOB_M0720$contrast <- paste("M_072022", AOB_M0720$contrast, sep="_")
AOB_D0720 <- read.csv("AOB_D0720_070923.csv")[,-1]
AOB_D0720$contrast <- paste("D_072022", AOB_D0720$contrast, sep="_")
AOB_K0720 <- read.csv("AOB_K0720_070923.csv")[,-1]
AOB_K0720$contrast <- paste("K_072022", AOB_K0720$contrast, sep="_")

AOB_M09 <- read.csv("AOB_M09_070923.csv")[,-1]
AOB_M09$contrast <- paste("M_091322", AOB_M09$contrast, sep="_")
AOB_D09 <- read.csv("AOB_D09_070923.csv")[,-1]
AOB_D09$contrast <- paste("D_091322", AOB_D09$contrast, sep="_")
AOB_K09 <- read.csv("AOB_K09_070923.csv")[,-1]
AOB_K09$contrast <- paste("K_091322", AOB_K09$contrast, sep="_")

glmT3s.pairwise.global.ALL <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705, AOB_M0720, AOB_D0720, AOB_K0720,
                                    AOB_M09, AOB_D09, AOB_K09)

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
meanotus<-phyloseq_average(aob.physeq_bulk1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
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

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rare/AOB_BulkSoil_rare/')
write.csv(rr, file = "AOB_RR_130923.csv")



####################################################################################################################################
# BULK SOIL-FOR PREVALENCE 80
#############################################################################################################################

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80/')
AOB_M04 <- read.csv("AOB_M04_130923.csv")[,-1]
AOB_M04$contrast <- paste("M_042822", AOB_M04$contrast, sep="_")
AOB_D04 <- read.csv("AOB_D04_130923.csv")[,-1]
AOB_D04$contrast <- paste("D_042822", AOB_D04$contrast, sep="_")
AOB_K04 <- read.csv("AOB_K04_130923.csv")[,-1]
AOB_K04$contrast <- paste("K_042822", AOB_K04$contrast, sep="_")

AOB_M06 <- read.csv("AOB_M06_130923.csv")[,-1]
AOB_M06$contrast <- paste("M_060122", AOB_M06$contrast, sep="_")
AOB_D06 <- read.csv("AOB_D06_130923.csv")[,-1]
AOB_D06$contrast <- paste("D_060122", AOB_D06$contrast, sep="_")
AOB_K06 <- read.csv("AOB_K06_130923.csv")[,-1]
AOB_K06$contrast <- paste("K_060122", AOB_K06$contrast, sep="_")

AOB_M0705 <- read.csv("AOB_M0705_130923.csv")[,-1]
AOB_M0705$contrast <- paste("M_070522", AOB_M0705$contrast, sep="_")
AOB_D0705 <- read.csv("AOB_D0705_130923.csv")[,-1]
AOB_D0705$contrast <- paste("D_070522", AOB_D0705$contrast, sep="_")
AOB_K0705 <- read.csv("AOB_K0705_130923.csv")[,-1]
AOB_K0705$contrast <- paste("K_070522", AOB_K0705$contrast, sep="_")

AOB_M0720 <- read.csv("AOB_M0720_130923.csv")[,-1]
AOB_M0720$contrast <- paste("M_072022", AOB_M0720$contrast, sep="_")
AOB_D0720 <- read.csv("AOB_D0720_130923.csv")[,-1]
AOB_D0720$contrast <- paste("D_072022", AOB_D0720$contrast, sep="_")
AOB_K0720 <- read.csv("AOB_K0720_130923.csv")[,-1]
AOB_K0720$contrast <- paste("K_072022", AOB_K0720$contrast, sep="_")

AOB_M09 <- read.csv("AOB_M09_130923.csv")[,-1]
AOB_M09$contrast <- paste("M_091322", AOB_M09$contrast, sep="_")
AOB_D09 <- read.csv("AOB_D09_130923.csv")[,-1]
AOB_D09$contrast <- paste("D_091322", AOB_D09$contrast, sep="_")
AOB_K09 <- read.csv("AOB_K09_130923.csv")[,-1]
AOB_K09$contrast <- paste("K_091322", AOB_K09$contrast, sep="_")

glmT3s.pairwise.global.ALL <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705, AOB_M0720, AOB_D0720, AOB_K0720,
                                    AOB_M09, AOB_D09, AOB_K09)

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
meanotus<-phyloseq_average(aob.physeq_bulk1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
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
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80/')
#write.csv(rr, file = "AOB_RR_prev80_130923.csv")

################################################################################
#HeatMap
#install.packages("colorRamp2")
library(colorRamp2)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# read the log2fold ratio data files
# bulk soil
#setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/log2fold/')
rr <- read.csv("AOB_RR_Bulk_130923.csv", row.names = 1)
names(rr)=str_sub(names(rr),4)

# rizosphere
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/log2fold/')
rr.rhizo <- read.csv("AOB_RR_Rhizo_130923.csv", row.names = 1)
names(rr.rhizo)=str_sub(names(rr.rhizo),4)
#Set annotation
#setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann <- read.csv("AOB.anno.csv", row.names = 1)

rr.ord <- rr[rownames(ann), ]
rr.rhizo.ord <- rr.rhizo[rownames(ann), ]

colours <- list("Taxonomy"=c("Nitrosospira-sp-17Nsp14_2671457573"="#2C85B2",
                            "Nitrosolobus-multiformis-Nl1_2667636517"="#990F0F",
                             "Nitrosospira-sp_2636913388"="#B2E5FF",
                             "Nitrosomonas-communis-Nm44_2676397764"="#FFB2B2",
                             "Nitrosospira-sp_2630434854"="#7EC3E5"))

colAnn <- rowAnnotation(df=ann,name = "Taxonomy",col=colours,
                            annotation_width=unit(c(1, 4), "cm"), 
                            gap=unit(1, "mm"))

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert <- read.csv("BulkSoil.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#ffcf20FF",
                             "D"="#541352FF",
                             "K"="#2f9aa0FF"))
colFert.Ann <- columnAnnotation(df=ann.fert, col=colours.fert,
                                show_legend =F,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"), 
                                gap=unit(1, "mm"))

# heatmap

col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
#tax = Heatmap(as.matrix(ann), cluster_rows  = F)
aob.bs.hm <- Heatmap(as.matrix(rr.ord),
                     name = "Log2-ratio",
                     column_title = "Bulk Soil",
                     #cluster_columns = F,
                     cluster_rows  = F,
                     column_order = order(colnames(as.matrix(rr.ord))),
                     #row_order = order(rownames(as.matrix(rr))),
                     #column_split = data.frame(rep(c("D", "K", "M"),5,5,5)),
                     #column_split = column_split,
                     #column_names_gp = gpar(fontsize=15, col = c(rep("#ffcf20FF", 5), rep("#541352FF", 5), rep("#2f9aa0FF", 5))),
                     #right_annotation = colAnn,
                     #column_names_gp = gpar(col = c(rep("red", 10), rep("blue", 8)))
                     #column_names_rot = 45,
                     bottom_annotation = colFert.Ann,
                     show_column_dend = F,
                     show_row_dend = F,
                     border_gp = gpar(col = "black", lty = 2),
                     col= col_fun)
aob.bs.hm

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert.rh <- read.csv("Rhizo.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#ffcf20FF",
                                       "D"="#541352FF",
                                       "K"="#2f9aa0FF"))
colFert.Ann.rh <- columnAnnotation(df=ann.fert.rh, col=colours.fert,
                                show_legend =F,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"), 
                                gap=unit(1, "mm"))
aob.rh.hm <- Heatmap(as.matrix(rr.rhizo.ord),
                     name = "Log2-ratio",
                     column_title = "Rhizosphere",
                     cluster_columns = F,
                     cluster_rows  = F,
                     #column_split = column_split,
                     column_order = order(colnames(as.matrix(rr.rhizo.ord))),
                     #column_names_gp = gpar(fontsize=15, col = c(rep("#ffcf20FF", 3), rep("#541352FF", 3), rep("#2f9aa0FF", 3))),
                     right_annotation = colAnn,
                     bottom_annotation = colFert.Ann.rh,
                     #column_names_rot = 45,
                     show_column_dend = F,
                     show_row_dend = F,
                     #heatmap_border = TRUE,
                     border_gp = gpar(col = "black", lty = 2),
                     col= col_fun)
aob.rh.hm
aob.hm <- aob.bs.hm+aob.rh.hm
aob.hm
aob.hm2 <- draw(aob.hm,column_title = "AOB", ht_gap = unit(0.5, "cm"))
                #align_heatmap_legend="heatmap_top", 
                #column_title_gp = gpar(fontsize = 16))
aob.hm2
# save image
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB/')
png("heatm.aob.tiff",width=12,height=5,units="in",res=1200)
aob.hm2
dev.off()

comb.hm <- aob.hm %v% aoa.hm
comb.hm

################################################################################
###compile 3 genes in one heatmap###
################################################################################

#### Bulk Soil #####

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
rr.comp <- read.csv("3genes.RR.csv", row.names = 1)
names(rr.comp)=str_sub(names(rr.comp),4)
#Set annotation
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.comp <- read.csv("3genes.anno.csv", row.names = 1)
#order rownames
rr.comp.ord <- rr.comp[rownames(ann.comp), ]
#remove the character before "_"
#rownames(rr.comp.ord) <- sub('.*_', '', rownames(rr.comp.ord))

#relative abund of bulk soil
# calculate relative abundance 
aob.asv.ra <- transform_sample_counts(aob.rare.1282.seq, function(x) x/sum(x))
aob.asv.ra
#aob.asv.ra.melt <- psmelt(aob.asv.ra)
aob.asv.ra.melt <- psmelt(aob.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
aob.asv.ra.melt$ra.perc <- aob.asv.ra.melt$Mean*100
# aoa
aoa.asv.ra <- transform_sample_counts(aoa.rare.min.physeq, function(x) x/sum(x))
aoa.asv.ra
aoa.asv.ra.melt <- psmelt(aoa.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
aoa.asv.ra.melt$ra.perc <- aoa.asv.ra.melt$Mean*100
# comammox
com.asv.ra <- transform_sample_counts(com.rare.min.physeq, function(x) x/sum(x))
com.asv.ra
com.asv.ra.melt <- psmelt(com.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
com.asv.ra.melt$ra.perc <- com.asv.ra.melt$Mean*100

# save the results in the computer and read the csv file
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
RA.comp <- read.csv("3genes.ra.all.csv", row.names = 1)
RA.comp.ord <- RA.comp[match(rownames(ann.comp), rownames(RA.comp)), ]



#set colors
lgd1 <- Legend(labels = c("Nitrosolobus-multiformis-Nl1_2667636517",
                          "Nitrosomonas-communis-Nm44_2676397764",
                          "Nitrosospira-sp-17Nsp14_2671457573",
                          "Nitrosospira-sp_2630434854",
                          "Nitrosospira-sp_2636913388"),
               legend_gp = gpar(fill=c("#990F0F","#FFB2B2","#2C85B2","#7EC3E5","#B2E5FF")),
               title= "AOB")

lgd2 <- Legend(labels = c("Ca.Nitrosotaleales (NT-Alpha-1.1.2.2)",
                          "Nitrososphaerales (NS-Delta-1.Incertae_sedis)",
                          "Nitrososphaerales (NS-Gamma-1.2)",
                          "Nitrososphaerales (NS-Gamma-2.3.1)"),
               legend_gp = gpar(fill=c("#B22C2C","#A3CC51","#E5FFB2","#B2E5FF")),
               title= "AOA")

lgd3 <- Legend(labels = c("Clade B Nitrospira-sp.GGF-bin22",
                          "Clade B Nitrospira-sp.LM-bin98",
                          "Clade B Nitrospira-sp.LPPL-bin249"),
               legend_gp = gpar(fill=c("#B26F2C","#CC8E51","#E5B17E")),
               title= "COMAMMOX")
pd = packLegend(lgd1, lgd2, lgd3, direction = "vertical")
draw(pd)
              
col.comp.ord <- list("Taxonomy"=c("Nitrosolobus-multiformis-Nl1_2667636517"="#990F0F",
                              "Nitrosomonas-communis-Nm44_2676397764"="#FFB2B2",
                              "Nitrosospira-sp-17Nsp14_2671457573"="#2C85B2",
                              "Nitrosospira-sp_2630434854"="#7EC3E5",
                              "Nitrosospira-sp_2636913388"="#B2E5FF",
                              "Ca.Nitrosotaleales (NT-Alpha-1.1.2.2)"="#B22C2C",
                              "Nitrososphaerales (NS-Delta-1.Incertae_sedis)"="#A3CC51",
                              "Nitrososphaerales (NS-Gamma-1.2)"="#E5FFB2",
                              "Nitrososphaerales (NS-Gamma-2.3.1)"="#B2E5FF",
                              "Clade B Nitrospira-sp.GGF-bin22"="#B26F2C",
                              "Clade B Nitrospira-sp.LM-bin98"="#CC8E51",
                              "Clade B Nitrospira-sp.LPPL-bin249"="#E5B17E"))
col_level <- factor(ann.comp$Taxonomy, levels = c("Nitrosolobus-multiformis-Nl1_2667636517",
                                                  "Nitrosomonas-communis-Nm44_2676397764",
                                                  "Nitrosospira-sp-17Nsp14_2671457573",
                                                  "Nitrosospira-sp_2630434854",
                                                  "Nitrosospira-sp_2636913388",
                                                  "Ca.Nitrosotaleales (NT-Alpha-1.1.2.2)",
                                                  "Nitrososphaerales (NS-Delta-1.Incertae_sedis)",
                                                  "Nitrososphaerales (NS-Gamma-1.2)",
                                                  "Nitrososphaerales (NS-Gamma-2.3.1)",
                                                  "Clade B Nitrospira-sp.GGF-bin22",
                                                  "Clade B Nitrospira-sp.LM-bin98",
                                                  "Clade B Nitrospira-sp.LPPL-bin249"))
tax_level=levels(col_level)

colAnn.comp <- rowAnnotation(df=ann.comp,
                             col=col.comp.ord,
                             show_legend =F,
                             annotation_legend_param = list(Taxonomy = list(
                             title="Taxonomy",
                             ncol=3,
                             at = tax_level)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
colAnn.comp

bar.ann.comp <- rowAnnotation(RelativeAbundance = anno_barplot(RA.comp.ord,
                                                  gp = gpar(fill = c("#990F0F", "#990F0F","#990F0F","#990F0F","#990F0F",
                                                                     "#FFB2B2","#2C85B2","#2C85B2","#2C85B2","#2C85B2",
                                                                     "#7EC3E5","#B2E5FF","#B2E5FF","#B2E5FF","#B2E5FF",
                                                                     "#B2E5FF","#B2E5FF","#B22C2C","#A3CC51","#A3CC51",
                                                                     "#E5FFB2","#B2E5FF","#B26F2C","#B26F2C","#CC8E51",
                                                                     "#E5B17E","#E5B17E","#E5B17E")),
                                                  ylim=c(0,0.18),
                                                  extend = 0.00000000000000001,
                                                  width  = unit(4, "cm"),
                                                  height = unit(6, "cm")))
                                                  #show_annotation_name =F)

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert <- read.csv("BulkSoil.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#ffcf20FF",
                                       "D"="#541352FF",
                                       "K"="#2f9aa0FF"))
colFert.Ann <- columnAnnotation(df=ann.fert, 
                                col=colours.fert,
                                show_legend =F,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"), 
                                gap=unit(1, "mm"))

col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
#tax = Heatmap(as.matrix(ann), cluster_rows  = F)
row_split = rep("AOB", 17)
row_split[18:22] = "AOA"
row_split[23:28] = "COMAMMOX"
row_split.fa = factor(row_split, levels = c("AOB", "AOA", "COMAMMOX"))
comp.bs.hm <- Heatmap(as.matrix(rr.comp.ord),
                     name = "Log2-ratio",
                     column_title = "Bulk Soil",
                     cluster_rows  = F,
                     cluster_row_slices=F,
                     column_order = order(colnames(as.matrix(rr.comp.ord))),
                     row_split = row_split.fa, 
                     left_annotation = bar.ann.comp,
                     bottom_annotation = colFert.Ann,
                     show_column_dend = F,
                     show_row_dend = F,
                     row_gap = unit(0.4, "cm"),
                     border_gp = gpar(col = "black", lty = 2),
                     col= col_fun)
comp.bs.hm
#decorate_annotation("RelativeAbundance", {
  grid.text("Relative Abundance",y = unit(-8.8,"cm"),just = "bottom")
})



#### Rhizosphere #####

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
rr.rhizo.comp <- read.csv("3genes.rhizos.RR.csv", row.names = 1)
names(rr.rhizo.comp)=str_sub(names(rr.rhizo.comp),4)
#Set annotation
setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.comp <- read.csv("3genes.anno.csv", row.names = 1)
#order rownames
rr.rhizo.comp.ord <- rr.rhizo.comp[rownames(ann.comp), ]

colAnn.comp <- rowAnnotation(df=ann.comp,
                             col=col.comp.ord,
                             show_legend =F,
                             annotation_legend_param = list(Taxonomy = list(
                               title="Taxonomy",
                               ncol=3,
                               at = tax_level)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert.rh <- read.csv("Rhizo.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#ffcf20FF",
                                       "D"="#541352FF",
                                       "K"="#2f9aa0FF"))
colFert.Ann.rh <- columnAnnotation(df=ann.fert.rh, col=colours.fert,
                                show_legend =F,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"), 
                                gap=unit(1, "mm"))

col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
row_split = rep("AOB", 17)
row_split[18:22] = "AOA"
row_split[23:28] = "COMAMMOX"
row_split.fa = factor(row_split, levels = c("AOB", "AOA", "COMAMMOX"))
comp.rh.hm <- Heatmap(as.matrix(rr.rhizo.comp.ord),
                      name = "Log2-ratio",
                      column_title = "Rhizosphere",
                      cluster_rows  = F,
                      cluster_row_slices=F,
                      column_order = order(colnames(as.matrix(rr.rhizo.comp.ord))),
                      row_split = row_split.fa, 
                      right_annotation = colAnn.comp,
                      bottom_annotation = colFert.Ann.rh,
                      show_column_dend = F,
                      show_row_dend = F,
                      row_gap = unit(0.4, "cm"),
                      border_gp = gpar(col = "black", lty = 2, width=3),
                      col= col_fun)

comp.rh.hm
comp.heat <- comp.bs.hm + comp.rh.hm
comp.heat2 <- draw(comp.heat,
                   ht_gap = unit(0.4, "cm"),
                    heatmap_legend_list=pd,
                    align_heatmap_legend="heatmap_top")
# save image
setwd('D:/Fina/INRAE_Project/microservices_fig/')
png("heatm.3genes4.tiff",width=14,height=7,units="in",res=1200)
comp.heat2
dev.off()









####################################################################################################################################
# RHIZOSPHERE SOIL-FOR PREVALENCE 80
##############################################################################################################################

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
AOB_M04 <- read.csv("AOB_M04.rh_130923.csv")[,-1]
AOB_M04$contrast <- paste("M_042822", AOB_M04$contrast, sep="_")
AOB_D04 <- read.csv("AOB_D04.rh_130923.csv")[,-1]
AOB_D04$contrast <- paste("D_042822", AOB_D04$contrast, sep="_")
AOB_K04 <- read.csv("AOB_K04.rh_130923.csv")[,-1]
AOB_K04$contrast <- paste("K_042822", AOB_K04$contrast, sep="_")

AOB_M06 <- read.csv("AOB_M06.rh_130923.csv")[,-1]
AOB_M06$contrast <- paste("M_060122", AOB_M06$contrast, sep="_")
AOB_D06 <- read.csv("AOB_D06.rh_130923.csv")[,-1]
AOB_D06$contrast <- paste("D_060122", AOB_D06$contrast, sep="_")
AOB_K06 <- read.csv("AOB_K06.rh_130923.csv")[,-1]
AOB_K06$contrast <- paste("K_060122", AOB_K06$contrast, sep="_")

AOB_M0705 <- read.csv("AOB_M0705.rh_130923.csv")[,-1]
AOB_M0705$contrast <- paste("M_070522", AOB_M0705$contrast, sep="_")
AOB_D0705 <- read.csv("AOB_D0705.rh_130923.csv")[,-1]
AOB_D0705$contrast <- paste("D_070522", AOB_D0705$contrast, sep="_")
AOB_K0705 <- read.csv("AOB_K0705.rh_130923.csv")[,-1]
AOB_K0705$contrast <- paste("K_070522", AOB_K0705$contrast, sep="_")

glmT3s.pairwise.global.ALL <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705)

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
meanotus<-phyloseq_average(aob.physeq_rh1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
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

ctrst.glm.CBFP.T3s.sub.ed <- 
  head(ctrst.glm.CBFP.T3s.sub)

# put the same column order in ctrst.glm.CBFP.T3s.sub and in meanotus
ctrst.glm.CBFP.T3s.sub<-ctrst.glm.CBFP.T3s.sub[,c(7,8,9,
                                                  1,2,3,
                                                  4,5,6)]
ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[rownames(meanotus), ]

# Multiply the matrices to get the RR when it is significant and 0 when it is not significant
rr<-meanotus*ctrst.glm.CBFP.T3s.sub

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
write.csv(rr, file = "AOB_RR_Rhizo_130923.csv")

####################################################################################################################################

##### 1. BULK SOIL - NOT RAREFIED #####

aob.physeq # still contain rhizosphere!, 1338 taxa 192 samples
sort(rowSums(otu_table(aob.physeq), na.rm = FALSE, dims = 1), decreasing = FALSE) # nothing is zero
aob.bulk.rawseq <- subset_samples(aob.physeq, Type=="BS")
aob.bulk.rawseq1 <- prune_taxa(taxa_sums(aob.bulk.rawseq)>0,aob.bulk.rawseq)
aob.bulk.rawseq1 # 1008 taxa

# Date: 04-28-2022

# 1. MINERAL
M04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="M")
M04rawseq1 <- prune_taxa(taxa_sums(M04rawseq)>0, M04rawseq)
M04_table <- as.data.frame(otu_table(M04rawseq1))
M04_table
cond.aldx <- sample_data(M04rawseq1)$Irrigation
# 2. BIODYNAMIC
D04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="D")
D04rawseq1 <- prune_taxa(taxa_sums(D04rawseq)>0, D04rawseq)
D04_table <- as.data.frame(otu_table(D04rawseq1))
D04_table
cond.aldx.D <- sample_data(D04rawseq1)$Irrigation
# 3. CONVENTIONAL
K04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="K")
K04rawseq1 <- prune_taxa(taxa_sums(K04rawseq)>0, K04rawseq)
K04_table <- as.data.frame(otu_table(K04rawseq1))
K04_table
cond.aldx.K <- sample_data(K04rawseq1)$Irrigation

# Date: 06-01-2022

# 1. MINERAL
M06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="M")
M06rawseq1 <- prune_taxa(taxa_sums(M06rawseq)>0, M06rawseq)
M06_table <- as.data.frame(otu_table(M06rawseq1))
M06_table
cond.aldx.M06 <- sample_data(M06rawseq1)$Irrigation
# 2. BIODYNAMIC
D06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="D")
D06rawseq1 <- prune_taxa(taxa_sums(D06rawseq)>0, D06rawseq)
D06_table <- as.data.frame(otu_table(D06rawseq1))
D06_table
# 3. CONVENTIONAL
K06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="K")
K06rawseq1 <- prune_taxa(taxa_sums(K06rawseq)>0, K06rawseq)
K06_table <- as.data.frame(otu_table(K06rawseq1))
K06_table
sample_data(K06rawseq1)$Irrigation

# Date: 07-05-2022

# 1. MINERAL
M0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="M")
M0705rawseq1 <- prune_taxa(taxa_sums(M0705rawseq)>0, M0705rawseq)
M0705_table <- as.data.frame(otu_table(M0705rawseq1))
M0705_table
cond.aldx.M0705 <- sample_data(M0705rawseq1)$Irrigation
# 2. BIODYNAMIC
D0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="D")
D0705rawseq1 <- prune_taxa(taxa_sums(D0705rawseq)>0, D0705rawseq)
D0705_table <- as.data.frame(otu_table(D0705rawseq1))
D0705_table
# 3. CONVENTIONAL
K0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="K")
K0705rawseq1 <- prune_taxa(taxa_sums(K0705rawseq)>0, K0705rawseq)
K0705_table <- as.data.frame(otu_table(K0705rawseq1))
K0705_table
sample_data(K0705rawseq1)$Irrigation

# Date: 07-20-2022

# 1. MINERAL
M0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="M")
M0720rawseq1 <- prune_taxa(taxa_sums(M0720rawseq)>0, M0720rawseq)
M0720_table <- as.data.frame(otu_table(M0720rawseq1))
M0720_table
cond.aldx.M0720 <- sample_data(M0720rawseq1)$Irrigation
# 2. BIODYNAMIC
D0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="D")
D0720rawseq1 <- prune_taxa(taxa_sums(D0720rawseq)>0, D0720rawseq)
D0720_table <- as.data.frame(otu_table(D0720rawseq1))
D0720_table
# 3. CONVENTIONAL
K0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="K")
K0720rawseq1 <- prune_taxa(taxa_sums(K0720rawseq)>0, K0720rawseq)
K0720_table <- as.data.frame(otu_table(K0720rawseq1))
K0720_table
sample_data(K0720rawseq1)$Irrigation

# Date: 09-13-2022

# 1. MINERAL
M09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="M")
M09rawseq1 <- prune_taxa(taxa_sums(M09rawseq)>0, M09rawseq)
M09_table <- as.data.frame(otu_table(M09rawseq1))
M09_table
cond.aldx.M09 <- sample_data(M09rawseq1)$Irrigation
# 2. BIODYNAMIC
D09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="D")
D09rawseq1 <- prune_taxa(taxa_sums(D09rawseq)>0, D09rawseq)
D09_table <- as.data.frame(otu_table(D09rawseq1))
D09_table
# 3. CONVENTIONAL
K09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="K")
K09rawseq1 <- prune_taxa(taxa_sums(K09rawseq)>0, K09rawseq)
K09_table <- as.data.frame(otu_table(K09rawseq1))
K09_table
sample_data(K09rawseq1)$Irrigation

##### 2. RHIZOSPHERE - NOT RAREFIED #####

aob.rh.rawseq <- subset_samples(aob.physeq, Type=="RS")
aob.rh.rawseq1 <- prune_taxa(taxa_sums(aob.rh.rawseq)>0,aob.rh.rawseq)
aob.rh.rawseq1 # 940 taxa 72 samples

# Date: 04-28-2022

# 1. MINERAL
M04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="M")
M04.rh.rawseq1 <- prune_taxa(taxa_sums(M04.rh.rawseq)>0, M04.rh.rawseq)
M04.rh_table <- as.data.frame(otu_table(M04.rh.rawseq1))
M04.rh_table
# 2. BIODYNAMIC
D04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="D")
D04.rh.rawseq1 <- prune_taxa(taxa_sums(D04.rh.rawseq)>0, D04.rh.rawseq)
D04.rh_table <- as.data.frame(otu_table(D04.rh.rawseq1))
D04.rh_table
# 3. CONVENTIONAL
K04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="K")
K04.rh.rawseq1 <- prune_taxa(taxa_sums(K04.rh.rawseq)>0, K04.rh.rawseq)
K04.rh_table <- as.data.frame(otu_table(K04.rh.rawseq1))
K04.rh_table

# Date: 06-01-2022

# 1. MINERAL
M06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="M")
M06.rh.rawseq1 <- prune_taxa(taxa_sums(M06.rh.rawseq)>0, M06.rh.rawseq)
M06.rh_table <- as.data.frame(otu_table(M06.rh.rawseq1))
M06.rh_table
# 2. BIODYNAMIC
D06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="D")
D06.rh.rawseq1 <- prune_taxa(taxa_sums(D06.rh.rawseq)>0, D06.rh.rawseq)
D06.rh_table <- as.data.frame(otu_table(D06.rh.rawseq1))
D06.rh_table
# 3. CONVENTIONAL
K06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="K")
K06.rh.rawseq1 <- prune_taxa(taxa_sums(K06.rh.rawseq)>0, K06.rh.rawseq)
K06.rh_table <- as.data.frame(otu_table(K06.rh.rawseq1))
K06.rh_table

# Date: 07-05-2022

# 1. MINERAL
M0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="M")
M0705.rh.rawseq1 <- prune_taxa(taxa_sums(M0705.rh.rawseq)>0, M0705.rh.rawseq)
M0705.rh_table <- as.data.frame(otu_table(M0705.rh.rawseq1))
M0705.rh_table
# 2. BIODYNAMIC
D0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="D")
D0705.rh.rawseq1 <- prune_taxa(taxa_sums(D0705.rh.rawseq)>0, D0705.rh.rawseq)
D0705.rh_table <- as.data.frame(otu_table(D0705.rh.rawseq1))
D0705.rh_table
# 3. CONVENTIONAL
K0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="K")
K0705.rh.rawseq1 <- prune_taxa(taxa_sums(K0705.rh.rawseq)>0, K0705.rh.rawseq)
K0705.rh_table <- as.data.frame(otu_table(K0705.rh.rawseq1))
K0705.rh_table

###################################################################################################################################
###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- M04rawseq1
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
physeq.subset.75 # 36 taxa
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
o = log(sample_sums(M04rawseq1))
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
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson',offset = o)
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

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_raw/')
AOB_M04 <- read.csv("AOB_M04_120923.csv")[,-1]
AOB_M04$contrast <- paste("M_042822", AOB_M04$contrast, sep="_")
AOB_D04 <- read.csv("AOB_D04_120923.csv")[,-1]
AOB_D04$contrast <- paste("D_042822", AOB_D04$contrast, sep="_")
AOB_K04 <- read.csv("AOB_K04_120923.csv")[,-1]
AOB_K04$contrast <- paste("K_042822", AOB_K04$contrast, sep="_")

AOB_M06 <- read.csv("AOB_M06_120923.csv")[,-1]
AOB_M06$contrast <- paste("M_060122", AOB_M06$contrast, sep="_")
AOB_D06 <- read.csv("AOB_D06_120923.csv")[,-1]
AOB_D06$contrast <- paste("D_060122", AOB_D06$contrast, sep="_")
AOB_K06 <- read.csv("AOB_K06_120923.csv")[,-1]
AOB_K06$contrast <- paste("K_060122", AOB_K06$contrast, sep="_")

AOB_M0705 <- read.csv("AOB_M0705_120923.csv")[,-1]
AOB_M0705$contrast <- paste("M_070522", AOB_M0705$contrast, sep="_")
AOB_D0705 <- read.csv("AOB_D0705_120923.csv")[,-1]
AOB_D0705$contrast <- paste("D_070522", AOB_D0705$contrast, sep="_")
AOB_K0705 <- read.csv("AOB_K0705_120923.csv")[,-1]
AOB_K0705$contrast <- paste("K_070522", AOB_K0705$contrast, sep="_")

AOB_M0720 <- read.csv("AOB_M0720_120923.csv")[,-1]
AOB_M0720$contrast <- paste("M_072022", AOB_M0720$contrast, sep="_")
AOB_D0720 <- read.csv("AOB_D0720_120923.csv")[,-1]
AOB_D0720$contrast <- paste("D_072022", AOB_D0720$contrast, sep="_")
AOB_K0720 <- read.csv("AOB_K0720_120923.csv")[,-1]
AOB_K0720$contrast <- paste("K_072022", AOB_K0720$contrast, sep="_")

AOB_M09 <- read.csv("AOB_M09_120923.csv")[,-1]
AOB_M09$contrast <- paste("M_091322", AOB_M09$contrast, sep="_")
AOB_D09 <- read.csv("AOB_D09_120923.csv")[,-1]
AOB_D09$contrast <- paste("D_091322", AOB_D09$contrast, sep="_")
AOB_K09 <- read.csv("AOB_K09_120923.csv")[,-1]
AOB_K09$contrast <- paste("K_091322", AOB_K09$contrast, sep="_")

glmT3s.pairwise.global.ALL.raw <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705, AOB_M0720, AOB_D0720, AOB_K0720,
                                    AOB_M09, AOB_D09, AOB_K09)

## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.ALL.raw$p.value <= 0.06)
table(glmT3s.pairwise.global.ALL.raw$p.adjust <= 0.06)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global.ALL.raw$OTU[glmT3s.pairwise.global.ALL.raw$p.adjust <= 0.06])
glmT3s.pairwise.global.signif.raw = glmT3s.pairwise.global.ALL.raw[glmT3s.pairwise.global.ALL.raw$p.adjust <=0.06,]

length(tmp_otu3s)
tmp_otu3s























