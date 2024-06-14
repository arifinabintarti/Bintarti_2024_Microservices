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
M04seq<- subset_samples(aob.physeq_bulk1, Date=="Apr 28th" & Treatment=="M")
M04seq1 <- prune_taxa(taxa_sums(M04seq)>0, M04seq) # 32 taxa after filtering (8)
sort(rowSums(otu_table(M04seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04seq<- subset_samples(aob.physeq_bulk1, Date=="Apr 28th" & Treatment=="D")
D04seq1 <- prune_taxa(taxa_sums(D04seq)>0, D04seq) # 20 taxa after filtering (13)
# 3. CONVENTIONAL
K04seq<- subset_samples(aob.physeq_bulk1, Date=="Apr 28th" & Treatment=="K")
K04seq1 <- prune_taxa(taxa_sums(K04seq)>0, K04seq) # 35 taxa after filtering (15)

# Date: 06-01-2022
# 1. MINERAL
M06seq<- subset_samples(aob.physeq_bulk1, Date=="Jun 1st" & Treatment=="M")
M06seq1 <- prune_taxa(taxa_sums(M06seq)>0, M06seq) # 55 taxa after filtering (18)
# 2. BIODYNAMIC
D06seq<- subset_samples(aob.physeq_bulk1, Date=="Jun 1st" & Treatment=="D")
D06seq1 <- prune_taxa(taxa_sums(D06seq)>0, D06seq) # 23 taxa after filtering (14)
# 3. CONVENTIONAL
K06seq<- subset_samples(aob.physeq_bulk1, Date=="Jun 1st" & Treatment=="K")
K06seq1 <- prune_taxa(taxa_sums(K06seq)>0, K06seq) # 61 taxa after filtering (31)

# Date: 07-05-2022
# 1. MINERAL
M0705seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 5th" & Treatment=="M")
M0705seq1 <- prune_taxa(taxa_sums(M0705seq)>0, M0705seq) # 59 taxa after filtering (22)
# 2. BIODYNAMIC
D0705seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 5th" & Treatment=="D")
D0705seq1 <- prune_taxa(taxa_sums(D0705seq)>0, D0705seq)# 36 taxa after filtering (24)
# 3. CONVENTIONAL
K0705seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 5th" & Treatment=="K")
K0705seq1 <- prune_taxa(taxa_sums(K0705seq)>0, K0705seq) # 54 taxa after filtering (26)

# Date: 07-20-2022
# 1. MINERAL
M0720seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 20th" & Treatment=="M")
M0720seq1 <- prune_taxa(taxa_sums(M0720seq)>0, M0720seq) # 40 taxa after filtering (19)
# 2. BIODYNAMIC
D0720seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 20th" & Treatment=="D")
D0720seq1 <- prune_taxa(taxa_sums(D0720seq)>0, D0720seq) # 35 taxa after filtering (21)
# 3. CONVENTIONAL
K0720seq<- subset_samples(aob.physeq_bulk1, Date=="Jul 20th" & Treatment=="K")
K0720seq1 <- prune_taxa(taxa_sums(K0720seq)>0, K0720seq) # 48 taxa after filtering (20)

# Date: 09-13-2022
# 1. MINERAL
M09seq<- subset_samples(aob.physeq_bulk1, Date=="Sept 13th" & Treatment=="M")
M09seq1 <- prune_taxa(taxa_sums(M09seq)>0, M09seq) # 46 taxa after filtering (23)
# 2. BIODYNAMIC
D09seq<- subset_samples(aob.physeq_bulk1, Date=="Sept 13th" & Treatment=="D")
D09seq1 <- prune_taxa(taxa_sums(D09seq)>0, D09seq) # 34 taxa after filtering (22)
# 3. CONVENTIONAL
K09seq<- subset_samples(aob.physeq_bulk1, Date=="Sept 13th" & Treatment=="K")
K09seq1 <- prune_taxa(taxa_sums(K09seq)>0, K09seq) # 47 taxa after filtering(24)
# total 625 taxa (75 prev); total 300 taxa (80 prev)

##### 2. RHIZOSPHERE - RAREFIED #####

aob.physeq_rh <- subset_samples(aob.rare.1282.seq, Type=="RS")
aob.physeq_rh # 1222 taxa 72 samples
aob.physeq_rh1 <- prune_taxa(taxa_sums(aob.physeq_rh)>0, aob.physeq_rh)
aob.physeq_rh1 # 831 taxa

# Date: 04-28-2022
# 1. MINERAL
M04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Apr 28th" & Treatment=="M")
M04.rh.seq1 <- prune_taxa(taxa_sums(M04.rh.seq)>0, M04.rh.seq) # 63 taxa after filtering (28)
sort(rowSums(otu_table(M04.rh.seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Apr 28th" & Treatment=="D")
D04.rh.seq1 <- prune_taxa(taxa_sums(D04.rh.seq)>0, D04.rh.seq) # 43 taxa after filtering (29)
# 3. CONVENTIONAL
K04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Apr 28th" & Treatment=="K")
K04.rh.seq1 <- prune_taxa(taxa_sums(K04.rh.seq)>0, K04.rh.seq) # 70 taxa after filtering (37)

# Date: 06-01-2022
# 1. MINERAL
M06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jun 1st" & Treatment=="M")
M06.rh.seq1 <- prune_taxa(taxa_sums(M06.rh.seq)>0, M06.rh.seq) # 65 taxa after filtering (27)
# 2. BIODYNAMIC
D06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jun 1st" & Treatment=="D")
D06.rh.seq1 <- prune_taxa(taxa_sums(D06.rh.seq)>0, D06.rh.seq) # 47 taxa after filtering (30)
# 3. CONVENTIONAL
K06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jun 1st" & Treatment=="K")
K06.rh.seq1 <- prune_taxa(taxa_sums(K06.rh.seq)>0, K06.rh.seq) # 71 taxa after filtering (30)

# Date: 07-05-2022
# 1. MINERAL
M0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jul 5th" & Treatment=="M")
M0705.rh.seq1 <- prune_taxa(taxa_sums(M0705.rh.seq)>0, M0705.rh.seq) # 60 taxa after filtering (21)
# 2. BIODYNAMIC
D0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jul 5th" & Treatment=="D")
D0705.rh.seq1 <- prune_taxa(taxa_sums(D0705.rh.seq)>0, D0705.rh.seq) # 41 taxa after filtering (24)
# 3. CONVENTIONAL
K0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="Jul 5th" & Treatment=="K")
K0705.rh.seq1 <- prune_taxa(taxa_sums(K0705.rh.seq)>0, K0705.rh.seq) # 67 taxa after filtering (29)
# total 527 taxa (prev 75); total 255 (prev 80)
################################################################################

###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <-aob.physeq_bulk1
physeq.subset 
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 428 BS, 423 RS

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
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$var3)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$var3)))))
view(df_otu_prev_ttt)
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$var3)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$var3 == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$var3 == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)
view(df_otu_prev_ttt)
# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 80,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75  # 68 taxa
#physeq.subset.AOB.BS <- physeq.subset.75
#physeq.subset.AOB.RS <- physeq.subset.75
#physeq.subset.AOB.BS=68 taxa and 119 samples; physeq.subset.AOB.RS=70 taxa and 72 samples

#setwd('D:/Fina/INRAE_Project/microservices/DAA/')
#write.csv(otu_table(physeq.subset.75), file = "aob.filt75.bs.tab.csv")

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
#library(emmeans)
#library(data.table)

otu_table(physeq.subset.75)
sample_data(physeq.subset.75)
tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(physeq.subset.75)) # using filtered data
#Since library sizes are different between groups, accounting for library size results in the gene no longer being DE at the 5% significance level. Not correcting for sequencing depth would thus result in spurious results.
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
    glmT3s.sum = summary(glmT3s)$coefficients$cond
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
#When there is only one thing to test, there is no multiplicity issue, and hence no multiplicity adjustment to the P values.
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
write.csv(glmT3s.pairwise.global, file = "AOB_K0705.rh_130923.csv")
aob.K0705.rh.fil <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob.K0705.rh.fil, file = "AOB_K0705.rh.tab_130923.csv")


#aob.log.dws.emm <- emmeans(aob.log.dws, ~ irrigation | fertilization*sampling.date)
#aob.log.dws.pair <- pairs(aob.log.dws.emm)



##############################################################################################################################
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare/')
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
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80//')
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80/')
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

glmT3s.pairwise.global.AOB <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705, AOB_M0720, AOB_D0720, AOB_K0720,
                                    AOB_M09, AOB_D09, AOB_K09)
dim(glmT3s.pairwise.global.AOB)
view(glmT3s.pairwise.global.AOB)
length(unique(glmT3s.pairwise.global.AOB$OTU)) #68
## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.AOB$p.value <= 0.05) # 40 # not unique
table(glmT3s.pairwise.global.AOB$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
#tmp_otu3s = unique(glmT3s.pairwise.global.AOB$OTU[glmT3s.pairwise.global.AOB$p.adjust <= 0.05])
tmp_otu3s.AOB = unique(glmT3s.pairwise.global.AOB$OTU[glmT3s.pairwise.global.AOB$p.value <= 0.05])
#glmT3s.pairwise.global.signif = glmT3s.pairwise.global.AOB[glmT3s.pairwise.global.AOB$p.adjust <=0.05,]
glmT3s.pairwise.global.signif = glmT3s.pairwise.global.AOB[glmT3s.pairwise.global.AOB$p.value <=0.05,]

length(tmp_otu3s.AOB) # 30 ASVs unique
view(tmp_otu3s)
colnames(glmT3s.pairwise.global.AOB)
# cast pvalues
#contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global.AOB[,c(10,1,2)]
contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global.AOB[,c(7,1,2)]
# numeric variable needs to be named "value" 
colnames(contrasts.glm.CBFP.T3s) <- c("value", "OTU_names", "contrast")
#contrasts.glm.CBFP.T3s <- subset(contrasts.glm.CBFP.T3s, (contrasts.glm.CBFP.T3s$OTU_names %in% BFPOTUs.T3snet.sig))
head(contrasts.glm.CBFP.T3s)
str(contrasts.glm.CBFP.T3s)
contrasts.glm.CBFP.T3s <- data.frame(cast(contrasts.glm.CBFP.T3s, contrast ~ OTU_names, value="value"))
str(contrasts.glm.CBFP.T3s)
rownames(contrasts.glm.CBFP.T3s) <- contrasts.glm.CBFP.T3s$contrast
contrasts.glm.CBFP.T3s$contrast <- NULL
view(contrasts.glm.CBFP.T3s)

# keep OTUs with at least one contrast <0.05 
contrasts.glm.CBFP.T3s.sub <- contrasts.glm.CBFP.T3s[,colSums(contrasts.glm.CBFP.T3s<0.05, na.rm=TRUE) >= 1]
dim(contrasts.glm.CBFP.T3s.sub)
head(contrasts.glm.CBFP.T3s.sub)
str(contrasts.glm.CBFP.T3s.sub)
view(contrasts.glm.CBFP.T3s.sub)
ctrst.glm.CBFP.T3s.sub <- data.frame(t(contrasts.glm.CBFP.T3s.sub))

# replace pvalues to 0 if non significant, or 1 if significant
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub == NA] <- 0
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >0.05] <- 2
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub <0.05] <- 1
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >1] <- 0
ctrst.glm.CBFP.T3s.sub[is.na(ctrst.glm.CBFP.T3s.sub)] <- 0
view(ctrst.glm.CBFP.T3s.sub)

# Calculate the OTU avg per treatment
# CHECK THE OBJECT
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
meanotus<-phyloseq_average(aob.physeq_bulk1,avg_type="arithmetic",acomp_zero_impute = NULL,group="var3")
meanotus<-as.data.frame(otu_table(meanotus));meanotus
view(meanotus)
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


view(meanotus)
# keep only columns containing log2fold ratios (RRs)
colnames(meanotus)
meanotus<-meanotus[,c(31:45)]
head(meanotus)
meanotus[meanotus == "-Inf"] <- 0
meanotus[meanotus == "Inf"] <- 0
meanotus[meanotus == "NaN"] <- 0


ctrst.glm.CBFP.T3s.sub.ed <- 
  head(ctrst.glm.CBFP.T3s.sub)
colnames(ctrst.glm.CBFP.T3s.sub)
colnames(meanotus)



# put the same column order in ctrst.glm.CBFP.T3s.sub and in meanotus
ctrst.glm.CBFP.T3s.sub<-ctrst.glm.CBFP.T3s.sub[,c(11,12,13,14,15,
                                                  1,2,3,4,5,
                                                  6,7,8,9,10)]
ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[rownames(meanotus), ]
view(ctrst.glm.CBFP.T3s.sub)
# Multiply the matrices to get the RR when it is significant and 0 when it is not significant
rr.aob<-meanotus*ctrst.glm.CBFP.T3s.sub
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80/')
#write.csv(rr, file = "AOB_RR_prev80_130923.csv")
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/AOB_BulkSoil_rare_prev80/')
write.csv(rr.aob, file = "AOB_RR_prev80_p.value_270224.csv")
################################################################################
#HeatMap
#install.packages("colorRamp2")
library(colorRamp2)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# read the log2fold ratio data files
# bulk soil
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/log2fold/')
#rr <- read.csv("AOB_RR_Bulk_130923.csv", row.names = 1)
rr <- read.csv("AOB_RR_Bulk_p.val_270224.csv", row.names = 1)
names(rr)=str_sub(names(rr),4)
View(rr)
# rizosphere
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/log2fold/')
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
#rr.rhizo <- read.csv("AOB_RR_Rhizo_130923.csv", row.names = 1)
rr.rhizo <- read.csv("AOB_RR_Rhizo_p.val_270224.csv", row.names = 1)
names(rr.rhizo)=str_sub(names(rr.rhizo),4)
#Set annotation
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#ann <- read.csv("AOB.anno.csv", row.names = 1)
ann <- read.csv("AOB.anno_prev80_pval_270224.csv",row.names = 1)
rr.ord <- rr[rownames(ann), ]
rr.rhizo.ord <- rr.rhizo[rownames(ann), ]

#colours <- list("Taxonomy"=c("Nitrosospira-sp-17Nsp14_2671457573"="#2C85B2",
                            #"Nitrosolobus-multiformis-Nl1_2667636517"="#990F0F",
                             #"Nitrosospira-sp_2636913388"="#B2E5FF",
                            #"Nitrosomonas-communis-Nm44_2676397764"="#FFB2B2",
                             #"Nitrosospira-sp_2630434854"="#7EC3E5",
                            #"Nitrosomonas-europaea-ATCC-19718_637427314"="#99540F"))
colours <- list("Taxonomy"=c("Nitrosolobus-multiformis-Nl1_2667636517"="#990F0F",
                              "Nitrosomonas-communis-Nm44_2676397764"="#FFB2B2",
                              "Nitrosomonas-europaea-ATCC-19718_637427314"="#99540F",
                              "Nitrosospira-sp-17Nsp14_2671457573"="#2C85B2",
                              "Nitrosospira-sp_2630434854"="#7EC3E5",
                              "Nitrosospira-sp_2636913388"="#B2E5FF"))
col_level <- factor(ann$Taxonomy, levels = c("Nitrosolobus-multiformis-Nl1_2667636517",
                                                  "Nitrosomonas-communis-Nm44_2676397764",
                                                  "Nitrosomonas-europaea-ATCC-19718_637427314",
                                                  "Nitrosospira-sp-17Nsp14_2671457573",
                                                  "Nitrosospira-sp_2630434854",
                                                  "Nitrosospira-sp_2636913388"))
tax_level=levels(col_level)
str(tax_level)
tax_level

colAnn <- rowAnnotation(df=ann,
                             col=colours,
                             show_legend =T,
                             annotation_legend_param = list(Taxonomy = list(
                             title="Taxonomy",
                             ncol=1,
                             at = tax_level)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
colAnn
#colAnn <- rowAnnotation(df=ann,name = "Taxonomy",col=colours,
                            #annotation_width=unit(c(1, 4), "cm"), 
                            #gap=unit(1, "mm"))

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
ann.fert <- read.csv("BulkSoil.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#E69F00",
                             "D"="#009E73",
                             "K"="#FF618C"))

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
                     right_annotation = colAnn,
                     #column_names_gp = gpar(col = c(rep("red", 10), rep("blue", 8)))
                     #column_names_rot = 45,
                     bottom_annotation = colFert.Ann,
                     show_column_dend = F,
                     show_row_dend = F,
                     border_gp = gpar(col = "black", lty = 2),
                     col= col_fun)
aob.bs.hm

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert.rh <- read.csv("Rhizo.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("M"="#E69F00",
                             "D"="#009E73",
                             "K"="#FF618C"))
colFert.Ann.rh <- columnAnnotation(df=ann.fert.rh, col=colours.fert,
                                show_legend =F,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"), 
                                gap=unit(1, "mm"))
aob.rh.hm <- Heatmap(as.matrix(rr.rhizo.ord),
                     name = "Log2-ratio",
                     column_title = "Rhizosphere",
                     #cluster_columns = F,
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
#png("heatm.aob.tiff",width=12,height=5,units="in",res=1200)
#aob.hm2
#dev.off()

comb.hm <- aob.hm %v% aoa.hm
comb.hm

################################################################################
###compile 3 genes in one heatmap###
################################################################################

#### Bulk Soil #####

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
#rr.comp <- read.csv("3genes.bulk.RR.csv", row.names = 1)
rr.comp <- read.csv("RR_3Genes_BS_280224.csv", row.names = 1)
head(rr.comp)
#names(rr.comp)=str_sub(names(rr.comp),4)
view(rr.comp)
#Set annotation
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#ann.comp <- read.csv("3genes.anno.csv", row.names = 1)
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB')
ann.comp <- read.csv("Anno_3Genes_280224.csv", row.names = 1)
view(ann.comp)
#order rownames
rr.comp.ord <- rr.comp[rownames(ann.comp), ]
view(rr.comp.ord)
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
View(aob.asv.ra.melt)
aob.asv.ra.melt$ra.perc <- aob.asv.ra.melt$Mean*100
aob.asv.ra.melt <- column_to_rownames(aob.asv.ra.melt, var = "OTU")
rownames(aob.asv.ra.melt)
# make asv and ra
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
ann <- read.csv("AOB.anno_prev80_pval_270224.csv",row.names = 1)
aob.ra <- ann
rownames(aob.ra)
aob.ra$RA <- aob.asv.ra.melt$ra.perc[match(row.names(aob.ra), row.names(aob.asv.ra.melt))]
view(aob.ra)
rownames(aob.ra) <- sub("ASV_", "AOB_ASV ", rownames(aob.ra))

# aoa
aoa.asv.ra <- transform_sample_counts(aoa.rare.min.physeq, function(x) x/sum(x))
aoa.asv.ra
aoa.asv.ra.melt <- psmelt(aoa.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
aoa.asv.ra.melt$ra.perc <- aoa.asv.ra.melt$Mean*100
aoa.asv.ra.melt <- column_to_rownames(aoa.asv.ra.melt, var = "OTU")
# make asv and ra
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
ann.aoa <- read.csv("AOA_anno_prev80_pval_280224.csv", row.names = 1)
aoa.ra <- ann.aoa
rownames(aoa.ra)
aoa.ra$RA <- aoa.asv.ra.melt$ra.perc[match(row.names(aoa.ra), row.names(aoa.asv.ra.melt))]
view(aoa.ra)
rownames(aoa.ra) <- sub("ASV_", "AOA_ASV ", rownames(aoa.ra))

# comammox
com.asv.ra <- transform_sample_counts(com.rare.min.physeq, function(x) x/sum(x))
com.asv.ra
com.asv.ra.melt <- psmelt(com.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
com.asv.ra.melt$ra.perc <- com.asv.ra.melt$Mean*100
view(com.asv.ra.melt)
com.asv.ra.melt <- column_to_rownames(com.asv.ra.melt, var = "OTU")
# make asv and ra
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB')
ann.com <- read.csv("COM_anno_prev80_pval_280224.csv", row.names = 1)
com.ra <- ann.com
rownames(com.ra)
com.ra$RA <- com.asv.ra.melt$ra.perc[match(row.names(com.ra), row.names(com.asv.ra.melt))]
view(com.ra)
rownames(com.ra) <- sub("ASV_", "COM_ASV ", rownames(com.ra))

# Join everything
RA.comp <- rbind(aob.ra, aoa.ra, com.ra)
RA.comp <- rownames_to_column(RA.comp, var="ASV")

RA.comp.df <- as.data.frame(RA.comp)
view(RA.comp.df)
RA.comp.df2 <- RA.comp.df[,-2]
view(RA.comp.df2)
RA.comp.df2 <- column_to_rownames(RA.comp.df2, var="ASV")
str(RA.comp.df2)
# save the results in the computer and read the csv file
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#RA.comp <- read.csv("3genes.ra.all.csv", row.names = 1)
#RA.comp.ord <- RA.comp[match(rownames(ann.comp), rownames(RA.comp)), ]
#view(RA.comp)

library(ComplexHeatmap)
colfunc <- colorRampPalette(c("#FF7F00","white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)
lgd1 <- Legend(labels = c("Nitrosolobus multiformis",
                              "Nitrosomonas communis",
                              "Nitrosomonas europaea",
                              "Nitrosospira sp."),
                              #"Nitrosospira-sp-17Nsp14_2671457573",
                              #"Nitrosospira-sp_2630434854",
                              #"Nitrosospira-sp_2636913388"),
               legend_gp = gpar(fill=c("#33A02C","#FF7F00","#FF9B38","#6A3D9A")),
               title= "AOB", labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20,fontface='bold'))

lgd2 <- Legend(labels = c("NS-Beta",
                          "NS-Delta",
                          #"NS-Delta",
                          "NS-Gamma",
                          #"NS-Gamma",
                          #"NS-Gamma",
                          "NT-Alpha"),
               legend_gp = gpar(fill=c("#DD701E","#E1823A","#E9A672","#7570B3")),
               #legend_gp = gpar(fill=c("#DD701E","#E1823A","#E1823A","#E9A672","#E9A672","#E9A672","#7570B3")),
               title= "AOA",labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20,fontface='bold'))

lgd3 <- Legend(labels = c("Clade A-Nitrospira sp.",
                                 "Clade B-Nitrospira sp."),
                                 #"Clade B-Nitrospira-sp.LM-bin98",
                                 #"Clade B-Nitrospira-sp.LPPL-bin249",
                                 #"Clade B-Nitrospira-sp.Smid-bin44"),
                          legend_gp = gpar(fill=c("#66C2A5","#FC8D62")),
               #legend_gp = gpar(fill=c("#66C2A5","#FC8D62","#FC8D62","#FC8D62","#FC8D62")),
               title= "COMAMMOX",labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20, fontface='bold'))
pd = packLegend(lgd,lgd1, lgd2, lgd3, direction = "horizontal")
draw(pd)
              
col.comp.ord <- list("Taxonomy"=c("Nitrosolobus-multiformis-Nl1_2667636517"="#33A02C",
                              "Nitrosomonas-communis-Nm44_2676397764"="#FF7F00",
                              "Nitrosomonas-europaea-ATCC-19718_637427314"="#FF9B38",
                              "Nitrosospira-sp-17Nsp14_2671457573"="#6A3D9A",
                              "Nitrosospira-sp_2630434854"="#6A3D9A",
                              "Nitrosospira-sp_2636913388"="#6A3D9A",
                                 "NS-Beta-1.3_OTU1_1-EU885554"="#DD701E",
                                 "NS-Delta-1.Incertae_sedis.2_OTU1_1-EU885561"="#E1823A",
                                 "NS-Delta-1.Incertae_sedis.2_OTU2_1-EU885632"="#E1823A",
                                 "NS-Gamma-1.2_OTU1_1-EU671146"="#E9A672",
                                 "NS-Gamma-1.Incertae_sedis_OTU1_1-EU651089"="#E9A672",
                                 "NS-Gamma-2.3.1_OTU2_1-KC469632"="#E9A672",
                                 "NT-Alpha-1.1.2.2-JN179533"="#7570B3",
                                 "Clade A-Nitrospira-sp.CTRL-LIN-TMP1"="#66C2A5",
                                 "Clade B-Nitrospira-sp.GGF-bin22"="#FC8D62",
                                 "Clade B-Nitrospira-sp.LM-bin98"="#FC8D62",
                                 "Clade B-Nitrospira-sp.LPPL-bin249"="#FC8D62",
                                 "Clade B-Nitrospira-sp.Smid-bin44"="#FC8D62"))
col_level <- factor(ann.comp$Taxonomy, levels = c("Nitrosolobus-multiformis-Nl1_2667636517",
                                                  "Nitrosomonas-communis-Nm44_2676397764",
                                                  "Nitrosomonas-europaea-ATCC-19718_637427314",
                                                  "Nitrosospira-sp-17Nsp14_2671457573",
                                                  "Nitrosospira-sp_2630434854",
                                                  "Nitrosospira-sp_2636913388",
                                                  "NS-Beta-1.3_OTU1_1-EU885554",
                                                  "NS-Delta-1.Incertae_sedis.2_OTU1_1-EU885561",
                                                  "NS-Delta-1.Incertae_sedis.2_OTU2_1-EU885632",
                                                  "NS-Gamma-1.2_OTU1_1-EU671146",
                                                  "NS-Gamma-1.Incertae_sedis_OTU1_1-EU651089",
                                                  "NS-Gamma-2.3.1_OTU2_1-KC469632",
                                                  "NT-Alpha-1.1.2.2-JN179533",
                                                  "Clade A-Nitrospira-sp.CTRL-LIN-TMP1",
                                 "Clade B-Nitrospira-sp.GGF-bin22",
                                 "Clade B-Nitrospira-sp.LM-bin98",
                                 "Clade B-Nitrospira-sp.LPPL-bin249",
                                 "Clade B-Nitrospira-sp.Smid-bin44"))

tax_level=levels(col_level)
str(tax_level)
tax_level

colAnn.comp <- rowAnnotation(annotation_name_gp= gpar(fontsize = 15),df=ann.comp,
                             col=col.comp.ord,
                             show_legend =F,
                             annotation_legend_param = list(Taxonomy = list(
                             title="Taxonomy",
                             #annotation_name_gp= gpar(fontsize = 20)
                             #ncol=1,
                             at = tax_level)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
colAnn.comp

bar.ann.comp <- rowAnnotation(annotation_name_gp= gpar(fontsize = 15),'Relative\nAbundance (%)' = anno_barplot(RA.comp.df2,
                                                               axis_param=list(gp=gpar(fontsize = 15)),
                                                  gp = gpar(fill = c("#33A02C", "#33A02C","#33A02C","#33A02C","#33A02C","#33A02C","#33A02C","#33A02C",
                                                                     "#FF7F00","#FF9B38","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     "#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     "#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     "#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     "#DD701E","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A",
                                                                     "#E1823A","#E1823A","#E9A672","#E9A672","#E9A672","#E9A672","#E9A672","#E9A672","#E9A672","#7570B3","#7570B3",
                                                                     "#66C2A5","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62",
                                                                     "#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62",
                                                                     "#FC8D62","#FC8D62","#FC8D62")),
                                                  ylim=c(0,0.18),
                                                  extend = 100,
                                                  width  = unit(4, "cm"),
                                                  height = unit(6, "cm")))
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB')
ann.fert <- read.csv("BulkSoil.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("CONMIN"="#E69F00",
                             "BIODYN"="#009E73",
                             "CONFYM"="#FF618C"))
#label = as.vector(c("group1","group2","group3"))
colFert.Ann <- columnAnnotation(df=ann.fert,
                                col=colours.fert,
                                show_legend =F,
                                #labels=label,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"),
                                annotation_name_gp= gpar(fontsize = 20),
                                gap=unit(1, "mm"))

col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
#tax = Heatmap(as.matrix(ann), cluster_rows  = F)
#row_split = rep("AOB", 17)
row_split = rep("AOB", 42)
#row_split[18:22] = "AOA"
row_split[43:65] = "AOA"
#row_split[23:28] = "COMAMMOX"
row_split[66:92] = "COMAMMOX"
row_split.fa = factor(row_split, levels = c("AOB", "AOA", "COMAMMOX"))
col_split = rep("CONM", 5)
col_split[6:10] = "BIOD"
col_split[11:15] = "CONF"
col_split.fa = factor(col_split, levels = c("BIOD", "CONF", "CONM"))
comp.bs.hm <- Heatmap(as.matrix(rr.comp.ord),
                     name = "Log2-ratio",
                     #column_title = "Bulk Soil",
                     cluster_rows  = F,
                     cluster_row_slices=F,
                     column_order=as.integer(c(6,7,8,9,10,11,12,13,14,15,1,2,3,4,5)),
                     #column_order = order(colnames(as.matrix(rr.comp.ord))),
                     row_split = row_split.fa, 
                     column_split = col_split.fa, 
                     left_annotation = bar.ann.comp,
                     #right_annotation = colAnn.comp,
                     #bottom_annotation = colFert.Ann,
                     show_column_dend = F,
                     show_row_dend = F,
                     row_gap = unit(0.4, "cm"),
                     column_gap = unit(0, "cm"),
                     row_names_gp = gpar(fontsize = 12),
                      row_title_gp = gpar(fontsize = 25),
                      column_title_gp = gpar(fontsize = 18),
                      column_names_gp = gpar(fontsize=20),
                     border_gp = gpar(col = "black", lty = 1),
                     show_heatmap_legend = F,
                    col= col_fun)
comp.bs.hm
#decorate_annotation("RelativeAbundance", {
  #grid.text("Relative Abundance",y = unit(-8.8,"cm"),just = "bottom")
#})

#### Rhizosphere #####
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.rhizo.comp <- read.csv("RR_3Genes_RS_280224.csv", row.names = 1)
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#rr.rhizo.comp <- read.csv("3genes.rhizos.RR.csv", row.names = 1)
#names(rr.rhizo.comp)=str_sub(names(rr.rhizo.comp),4)
#Set annotation
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.comp <- read.csv("Anno_3Genes_280224.csv", row.names = 1)
#order rownames
rr.rhizo.comp.ord <- rr.rhizo.comp[rownames(ann.comp), ]
view(rr.rhizo.comp.ord)
colAnn.comp <- rowAnnotation(annotation_name_gp= gpar(fontsize = 15),df=ann.comp,
                             col=col.comp.ord,
                             show_legend =F,
                             annotation_legend_param = list(Taxonomy = list(
                               title="Taxonomy",
                               ncol=1,
                               at = tax_level)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert.rh <- read.csv("Rhizo.anno.csv", row.names = 1)
colours.fert <- list("Fertilization"=c("CONMIN"="#E69F00",
                             "BIODYN"="#009E73",
                             "CONFYM"="#FF618C"))
colFert.Ann.rh <- columnAnnotation(df=ann.fert.rh, 
                                   #anno_block(labels = ann.fert.rh$Fertilization),
                                   col=colours.fert,
                                   show_legend =F,
                                   show_annotation_name =F,
                                   annotation_width=unit(c(1, 4), "cm"),
                                   annotation_name_gp= gpar(fontsize = 20),
                                   gap=unit(1, "mm"))

library(colorRamp2)
col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
lgd_fun = Legend(title="Log2-ratio",
                 col_fun = col_fun, labels_gp = gpar(fontsize=18), 
                 title_gp = gpar(fontsize=20,fontface='bold'))
lgd = Legend(title="Log2-ratio",col_fun = col_fun, labels_gp = gpar(fontsize=18),
             direction = "horizontal", 
             title_position = "topcenter",
             title_gp = gpar(fontsize=20,fontface='bold'))
#row_split = rep("AOB", 17)
#row_split[18:22] = "AOA"
#row_split[23:28] = "COMAMMOX"
row_split = rep("AOB", 42)
row_split[43:65] = "AOA"
row_split[66:92] = "COMAMMOX"
row_split.fa = factor(row_split, levels = c("AOB", "AOA", "COMAMMOX"))
col_split.rh = rep("CONM", 3)
col_split.rh[4:6] = "BIOD"
col_split.rh[7:9] = "CONF"
col_split.fa.rh = factor(col_split.rh, levels = c("BIOD", "CONF", "CONM"))
comp.rh.hm <- Heatmap(as.matrix(rr.rhizo.comp.ord),
                      #rect_gp = gpar(col = "grey", lwd = 2),
                      name = "Log2-ratio",
                      #column_title = "Rhizosphere",
                      cluster_rows  = F,
                      cluster_row_slices=F,
                      column_order=as.integer(c(4,5,6,7,8,9,1,2,3)),
                      #column_order = order(colnames(as.matrix(rr.rhizo.comp.ord))),
                      row_split = row_split.fa, 
                      column_split = col_split.fa.rh, 
                      right_annotation = colAnn.comp,
                      #bottom_annotation = colFert.Ann.rh,
                      show_column_dend = F,
                      show_row_dend = F,
                      row_gap = unit(0.4, "cm"),
                      column_gap = unit(0, "cm"),
                      row_names_gp = gpar(fontsize = 12),
                      row_title_gp = gpar(fontsize = 25),
                      column_title_gp = gpar(fontsize = 18),
                      column_names_gp = gpar(fontsize=20),
                      border_gp = gpar(col = "black", lty = 1),#width=0.5),
                      show_heatmap_legend = F,
                      col= col_fun)

comp.rh.hm
comp.heat <- comp.bs.hm + comp.rh.hm
comp.heat2 <- draw(comp.heat,
                   ht_gap = unit(1, "cm"),
                    heatmap_legend_list=pd,
                    heatmap_legend_side = "bottom")
                    #align_heatmap_legend="global_center")
comp.heat2

# save image
setwd('D:/Fina/INRAE_Project/microservices_fig/')
png("heatm.3genes4.tiff",width=14,height=7,units="in",res=1200)
comp.heat2
dev.off()

setwd('/Users/arifinabintarti/Documents/France/Figures/')
#png("Fig.5.tiff",width=19.7,height=13,units="in",res=1200)
png("Fig.5.tiff",width=12,height=17,units="in",res=1200)
comp.heat2
dev.off()

setwd('/Users/arifinabintarti/Documents/France/Figures/')
png("Fig.5_2.tiff",width=19,height=16,units="in",res=1200)
#png("heatm.3genes2.tiff",width=25,height=17,units="in",res=1200)
comp.heat2
dev.off()



####################################################################################################################################
# RHIZOSPHERE SOIL-FOR PREVALENCE 80
##############################################################################################################################
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
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

glmT3s.pairwise.global.AOB.RS <- rbind(AOB_M04, AOB_D04, AOB_K04, AOB_M06, AOB_D06, AOB_K06,
                                    AOB_M0705, AOB_D0705, AOB_K0705)
length(unique(glmT3s.pairwise.global.AOB.RS$OTU)) #70

## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.AOB.RS$p.value <= 0.05)
table(glmT3s.pairwise.global.AOB.RS$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
#tmp_otu3s = unique(glmT3s.pairwise.global.AOB.RS$OTU[glmT3s.pairwise.global.AOB.RS$p.adjust <= 0.05])
#glmT3s.pairwise.global.signif = glmT3s.pairwise.global.AOB.RS[glmT3s.pairwise.global.AOB.RS$p.adjust <=0.05,]
tmp_otu3s = unique(glmT3s.pairwise.global.AOB.RS$OTU[glmT3s.pairwise.global.AOB.RS$p.value <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global.AOB.RS[glmT3s.pairwise.global.AOB.RS$p.value <=0.05,]

length(tmp_otu3s) # 25 unique
tmp_otu3s

# cast pvalues
contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global.AOB.RS[,c(7,1,2)]
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
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub == NA] <- 0
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

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
#write.csv(rr, file = "AOB_RR_Rhizo_130923.csv")
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/AOB_Rhizo_rare_prev80/')
write.csv(rr, file = "AOB_RR_Rhizo_prev80_p.value_270224.csv")
####################################################################################################################################

##### 1. BULK SOIL - NOT RAREFIED #####

aob.physeq # still contain rhizosphere!, 1338 taxa 192 samples
sort(rowSums(otu_table(aob.physeq), na.rm = FALSE, dims = 1), decreasing = FALSE) # nothing is zero
aob.bulk.rawseq <- subset_samples(aob.physeq, Type=="BS")
aob.bulk.rawseq1 <- prune_taxa(taxa_sums(aob.bulk.rawseq)>0,aob.bulk.rawseq)
aob.bulk.rawseq1 # 1008 taxa

# Date: 04-28-2022

# 1. MINERAL
M04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Apr 28th" & Treatment=="M")
M04rawseq1 <- prune_taxa(taxa_sums(M04rawseq)>0, M04rawseq)
M04_table <- as.data.frame(otu_table(M04rawseq1))
M04_table
cond.aldx <- sample_data(M04rawseq1)$Irrigation
# 2. BIODYNAMIC
D04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Apr 28th" & Treatment=="D")
D04rawseq1 <- prune_taxa(taxa_sums(D04rawseq)>0, D04rawseq)
D04_table <- as.data.frame(otu_table(D04rawseq1))
D04_table
cond.aldx.D <- sample_data(D04rawseq1)$Irrigation
# 3. CONVENTIONAL
K04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Apr 28th" & Treatment=="K")
K04rawseq1 <- prune_taxa(taxa_sums(K04rawseq)>0, K04rawseq)
K04_table <- as.data.frame(otu_table(K04rawseq1))
K04_table
cond.aldx.K <- sample_data(K04rawseq1)$Irrigation

# Date: 06-01-2022

# 1. MINERAL
M06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jun 1st" & Treatment=="M")
M06rawseq1 <- prune_taxa(taxa_sums(M06rawseq)>0, M06rawseq)
M06_table <- as.data.frame(otu_table(M06rawseq1))
M06_table
cond.aldx.M06 <- sample_data(M06rawseq1)$Irrigation
# 2. BIODYNAMIC
D06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jun 1st" & Treatment=="D")
D06rawseq1 <- prune_taxa(taxa_sums(D06rawseq)>0, D06rawseq)
D06_table <- as.data.frame(otu_table(D06rawseq1))
D06_table
# 3. CONVENTIONAL
K06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jun 1st" & Treatment=="K")
K06rawseq1 <- prune_taxa(taxa_sums(K06rawseq)>0, K06rawseq)
K06_table <- as.data.frame(otu_table(K06rawseq1))
K06_table
sample_data(K06rawseq1)$Irrigation

# Date: 07-05-2022

# 1. MINERAL
M0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 5th" & Treatment=="M")
M0705rawseq1 <- prune_taxa(taxa_sums(M0705rawseq)>0, M0705rawseq)
M0705_table <- as.data.frame(otu_table(M0705rawseq1))
M0705_table
cond.aldx.M0705 <- sample_data(M0705rawseq1)$Irrigation
# 2. BIODYNAMIC
D0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 5th" & Treatment=="D")
D0705rawseq1 <- prune_taxa(taxa_sums(D0705rawseq)>0, D0705rawseq)
D0705_table <- as.data.frame(otu_table(D0705rawseq1))
D0705_table
# 3. CONVENTIONAL
K0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 5th" & Treatment=="K")
K0705rawseq1 <- prune_taxa(taxa_sums(K0705rawseq)>0, K0705rawseq)
K0705_table <- as.data.frame(otu_table(K0705rawseq1))
K0705_table
sample_data(K0705rawseq1)$Irrigation

# Date: 07-20-2022

# 1. MINERAL
M0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 20th" & Treatment=="M")
M0720rawseq1 <- prune_taxa(taxa_sums(M0720rawseq)>0, M0720rawseq)
M0720_table <- as.data.frame(otu_table(M0720rawseq1))
M0720_table
cond.aldx.M0720 <- sample_data(M0720rawseq1)$Irrigation
# 2. BIODYNAMIC
D0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 20th" & Treatment=="D")
D0720rawseq1 <- prune_taxa(taxa_sums(D0720rawseq)>0, D0720rawseq)
D0720_table <- as.data.frame(otu_table(D0720rawseq1))
D0720_table
# 3. CONVENTIONAL
K0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Jul 20th" & Treatment=="K")
K0720rawseq1 <- prune_taxa(taxa_sums(K0720rawseq)>0, K0720rawseq)
K0720_table <- as.data.frame(otu_table(K0720rawseq1))
K0720_table
sample_data(K0720rawseq1)$Irrigation

# Date: 09-13-2022

# 1. MINERAL
M09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Sept 13th" & Treatment=="M")
M09rawseq1 <- prune_taxa(taxa_sums(M09rawseq)>0, M09rawseq)
M09_table <- as.data.frame(otu_table(M09rawseq1))
M09_table
cond.aldx.M09 <- sample_data(M09rawseq1)$Irrigation
# 2. BIODYNAMIC
D09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Sept 13th" & Treatment=="D")
D09rawseq1 <- prune_taxa(taxa_sums(D09rawseq)>0, D09rawseq)
D09_table <- as.data.frame(otu_table(D09rawseq1))
D09_table
# 3. CONVENTIONAL
K09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="Sept 13th" & Treatment=="K")
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
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
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

# combine different 





















