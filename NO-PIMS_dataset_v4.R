library(phyloseq)
library(picante)
# setwd("C:/Users/eperezvaler/OneDrive/GD/WORK/02.PROYECTOS/2023_NO-PIMS/DATA/03.ProcessedData/01.16SData2R/")
# setwd("/home/eduardo/GDrive/WORK/02.PROYECTOS/2023_NO-PIMS/DATA/03.ProcessedData/01.16SData2R/")

# import all data
physeq<-import_biom(BIOMfilename = "sampled_otu_table.biom",treefilename = "phylogeny.tre")
metadata<-read.table("sampled_mapping_file_GENES.tsv",header = T,row.names = 1)
sample_names(physeq)<-paste0("S",sample_names(physeq))
sample_data(physeq)<-metadata

# This step is to make the script working the same was as it was shared
sample_data(physeq)$Treatment<-metadata$Genotype
sample_data(physeq)$mesh_size_um<-metadata$Genotype
sample_data(physeq)$sample_name<-row.names(sample_data(physeq))

sample_data(physeq)$Treatment<-factor(sample_data(physeq)$Treatment,levels=c("No","Col0","Nia1Nia2","Ahb1","Nox1","GSNOR1"))
physeq<-subset_samples(physeq,row.names(sample_data(physeq))!="S5"&row.names(sample_data(physeq))!="S29")

physeq<-subset_taxa(physeq, Rank5!="f__Mitochondria")
physeq<-subset_taxa(physeq, Rank5!="f__Chloroplast")


physeqr<-rarefy_even_depth(physeq,sample.size = 19000,replace = F,rngseed = 200)


##########################################
## BASIC STATS
##########################################
sum(otu_table(physeq))
#  [1] 2329978

mean(colSums(otu_table(physeq)))
# [1] 40172.03

str(otu_table(physeq))
# 58 samples

mean(colSums(otu_table(physeq)))/sqrt(58)
# [1] 5274.847


################################
# Filter low-abundant taxa
###############################
physeq.subset <- physeqr
data.obs <- as.data.frame(otu_table(physeq.subset))

### keeping OTUs with at least 0.02 % relative abundance across all samples
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]

new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incroporate into phyloseq Object


physeq.subset
# # it means only 793 taxa remain in the dataset after filtering. 


# # KEEPING THE SUM OF NON-INCLUDED TAXA AS OTU-104
# physeq.subset <- physeq
# data.obs <- as.data.frame(otu_table(physeq.subset))
# 
# ### keeping OTUs with at least 0.02 % relative abundance across all samples 
# keep.taxa.id <- which((rowSums(data.obs)/sum(data.obs)) > 0.0002)
# removed.taxa.id <- which((rowSums(data.obs)/sum(data.obs)) <= 0.0002)
# removed.taxa.abundances <- colSums(data.obs[removed.taxa.id, ])
# 
# data.F <- data.obs[keep.taxa.id, , drop = FALSE] 
# data.F <- rbind(data.F, removed.taxa.abundances)
# 
# row.names(data.F)[794]<-"OTU-104"
# 
# 
# new.otu <- as.matrix(data.F) # convert it into a matrix.
# new.otu <- otu_table(new.otu, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
# 
# otu_table(physeq.subset) <- new.otu # incroporate into phyloseq Object
# 
# physeq.subset


########################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
########################################################################################

prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset)
)
prevalencedf[1:10,]
#write.table(x=prevelancedf, file="Filtered_OTUtable-prevalence.csv")
dim(prevalencedf)


### calculate prevalence /!\ takes from 30min to 3h /!\

ps = physeq.subset 

df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE

df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Treatment)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Treatment)))))
#attention il ya Sample et sample

for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Treatment)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Treatment == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Treatment == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# write.csv(df_otu_prev_ttt, file = "df_otu_prev_ttt.csv")


#____________________________
### filtre otu par prevalence

ps =  physeq.subset 
df_prev = df_otu_prev_ttt

tmp_otu_F = rownames(df_prev[df_prev$max_prev > 60,])

physeq.subset <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)

# physeqT1F.Bacteria <- physeq.subset #416 #72 samples
#write.csv(as.data.frame(otu_table(physeqT1F)), file = "otu_table_filteredT1.csv")

rm(ps,df_prev,tmp_otu_F)


physeq.subset


####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
library(tidyr)
library(phyloseq)
library(lme4)
library(emmeans)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(MASS)
library(broom)
library(glmmTMB)
library(DHARMa)
library(stringr)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(reshape)
library(pals)
library(rstatix)
library(UpSetR)

# 

tmp_T3s <- physeq.subset

str(tmp_T3s)

#  treatment
a = tibble("sample"= tmp_T3s@sam_data$sample_name,
           "treatment"= as.character(tmp_T3s@sam_data$mesh_size_um))
# force control as intercept
a[a == "Col0"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.factor(tmp_T3s@sam_data$sample_name)
#tmp_T3s@sam_data$block <- paste(c("b"),tmp_T3s@sam_data$block, sep="")
# x <- as.factor(tmp_T3s@sam_data$block)



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
    glmT3s <- glmer(y ~ a + (1 | z), family='poisson', offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    #glmT3s.sum = summary(glmT3s)$coefficients$cond
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
    glmT3s.sum
    
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    
    ### multiple comparaison
    
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
    tmp_df = tmp_df[tmp_df$a == "1a" | tmp_df$b == "1a" ,]
    
    # adjust pval n=nb of contrasts
    #tmp_df[,"p.adjust"] <- p.adjust(tmp_df$p.value,"fdr",n=21)
    #tmp_df[,"p.adjust"] <- p.adjust(tmp_df$p.value,"bonferroni",n=21)
    
    tmp_df = cbind("OTU"=OTU,tmp_df)
    
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  rm(OTU,y,glmT3s,glmT3s.sum)
  
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global

glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value < glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value > glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value == glmT3s.pairwise.global$p.adjust,])



## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global$p.value <= 0.05)
table(glmT3s.pairwise.global$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global$OTU[glmT3s.pairwise.global$p.adjust <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global[glmT3s.pairwise.global$p.adjust <=0.05,]

length(tmp_otu3s) 

######### OJO
# save(glmT3s.pairwise.global, glmT3s.model.global, file = "glm-T3s.RData")
# load("/home/mathilde/Bureau/STAT/R-objects/glm-T3s.RData")



# cast pvalues
contrasts.glm.CBFP.T3s <- glmT3s.pairwise.global[,c(10,1,2)]
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
contrasts.glm.CBFP.T3s.sub <- contrasts.glm.CBFP.T3s[,colSums(contrasts.glm.CBFP.T3s<0.05) >= 1]
dim(contrasts.glm.CBFP.T3s.sub)
head(contrasts.glm.CBFP.T3s.sub)
str(contrasts.glm.CBFP.T3s.sub)

ctrst.glm.CBFP.T3s.sub <- data.frame(t(contrasts.glm.CBFP.T3s.sub))

ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >0.05] <- 1
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub !=1] <- 0
head(ctrst.glm.CBFP.T3s.sub)
colSums(ctrst.glm.CBFP.T3s.sub)


upset.CBFP.T3s <- UpSetR::upset(ctrst.glm.CBFP.T3s.sub, sets = colnames(ctrst.glm.CBFP.T3s.sub),
                                order.by = "freq", keep.order = T)
upset.CBFP.T3s
ctrst.glm.CBFP.T3s.sub[rowSums(ctrst.glm.CBFP.T3s.sub) == 0,]

T3s.sig.OTUs <- rownames(ctrst.glm.CBFP.T3s.sub)




# need to add the sign of the contrast
#load("/home/mathilde/Bureau/STAT/R-objects/BFPOTUs.T3snet.RData") # subset of OTUs used for network : > 2% for bacteria, > 20 occ for all
#BFPOTUs.T3snet.sig <- (intersect(sigOTUs.BFP.T3s, BFPOTUs.T3snet))
#BFPOTUs.T3snet.sig <- (intersect(tmp_otu2, BFPOTUs.T3snet))

tmp0T3s = psmelt(tmp_T3s)
#tmp0T3s <- subset(tmp0T3s, (tmp0T3s$OTU %in% BFPOTUs.T3snet.sig))
tmp0T3s <- tmp0T3s[order(tmp0T3s$OTU),]
# OTU names are like B_OTU-10 with B for bacteria, so you can split the B in another column like this :
tmp0T3s$group <- str_split_fixed(tmp0T3s$OTU, "_", 2)[,1]
tmp = tmp0T3s[,c("OTU","mesh_size_um","Abundance")]  
tmp$mesh_size_um <- as.character(tmp$mesh_size_um)

# load(file="CBFPT3s.mesh.sum.filtered.RData") # from SADs.R
# load(file="C:/Users/eperezvaler/OneDrive/GD/WORK/Scripts/Mathilde_differential_abundance/CBFPT3s.mesh.sum.filtered.RData")


metadata<-as(sample_data(physeq),"data.frame")

CBFPT3s.mesh.sum<-aggregate(rowSums(t(otu_table(physeq)))~Genotype,data=metadata,sum)
colnames(CBFPT3s.mesh.sum)<-c("Genotype","B")
row.names(CBFPT3s.mesh.sum)<-CBFPT3s.mesh.sum$Genotype
CBFPT3s.mesh.sum$mesh_size_um<-CBFPT3s.mesh.sum$Genotype

tmp1T3s <- left_join(tmp, CBFPT3s.mesh.sum) 

# tmp1T3s <- tmp1T3s[order(tmp1T3s$group),]
# tmp1T3s$gp_sum <- c(tmp1T3s[tmp1T3s$group == "B",]$B, tmp1T3s[tmp1T3s$group == "F",]$F, tmp1T3s[tmp1T3s$group == "P",]$P)
tmp1T3s$gp_sum<-tmp1T3s$B

str(tmp1T3s )

################################
# load("C:/Users/eperezvaler/OneDrive/GD/WORK/Scripts/Mathilde_differential_abundance/CBFPT3s.mesh.OTUsum.filtered.RData")

physeq.subset.merge<-merge_samples(physeq.subset, "mesh_size_um")
CBFPT3s.mesh.OTUsum<-as.data.frame(rowSums(decostand(otu_table(physeq.subset.merge),"pa")))
colnames(CBFPT3s.mesh.OTUsum)<-c("B")
CBFPT3s.mesh.OTUsum$mesh_size_um<-row.names(CBFPT3s.mesh.OTUsum)

# this is to merge dataframes based on one column
tmp1T3s0 <- left_join(tmp, CBFPT3s.mesh.OTUsum) 

# tmp1T3s0 <- tmp1T3s0[order(tmp1T3s0$group),]
tmp1T3s0$gp_OTUsum <- tmp1T3s0$B
# tmp1T3s0$B <- NULL
# tmp1T3s0$F <- NULL
# tmp1T3s0$P <- NULL

str(tmp1T3s0 )
# tmp1T3s <- left_join(tmp1T3s, tmp1T3s0)
tmp1T3s$gp_OTUsum<-tmp1T3s0$gp_OTUsum

# str(tmp1T3s)
# sometimes OTU names had a dot instead of a dash in some objects...
tmp1T3s$OTU <- gsub("-", "\\.", tmp1T3s$OTU)
# lo muevo abajo
# tmp1T3s <- tmp1T3s[!duplicated(tmp1T3s),]

# compute RRs -----------------
# relative abundance by treatment

relabund_T3s = as_tibble(tmp1T3s %>% group_by(OTU,mesh_size_um, gp_sum) %>% summarise(sum=sum(Abundance), 
                                                                                             avg=mean(Abundance)) %>% summarise(rel_abund=(sum*100/gp_sum), avg=avg)) 
tmp1T3s <- tmp1T3s[!duplicated(tmp1T3s),]
#relabund_T3s = as_tibble(tmp0 %>% group_by(OTU,mesh_size_um) %>% summarise(avg=mean(Abundance)))



head(relabund_T3s)
relabund_T3s <- as.data.frame(relabund_T3s)
str(relabund_T3s)
relabund_T3s$mesh_size_um <- as.factor(relabund_T3s$mesh_size_um)


# cast is like the opposite of melt
meanT3s_mesh <- data.frame(cast(relabund_T3s, mesh_size_um ~ OTU, value="avg"))
str(meanT3s_mesh)
rownames(meanT3s_mesh) <- meanT3s_mesh$mesh_size_um
meanT3s_mesh$mesh_size_um <- NULL
meanT3s_mesh <- data.frame(t(meanT3s_mesh))
str(meanT3s_mesh)

meanT3s_mesh$RR_No <- log2(meanT3s_mesh$No / meanT3s_mesh$Col0)
meanT3s_mesh$RR_Ahb1 <- log2(meanT3s_mesh$Ahb1 / meanT3s_mesh$Col0)
meanT3s_mesh$RR_Nia1Nia2 <- log2(meanT3s_mesh$Nia1Nia2 / meanT3s_mesh$Col0)
meanT3s_mesh$RR_Nox1 <- log2(meanT3s_mesh$Nox1 / meanT3s_mesh$Col0)
meanT3s_mesh$RR_GSNOR1 <- log2(meanT3s_mesh$GSNOR1 / meanT3s_mesh$Col0)



meanT3s_mesh[meanT3s_mesh == "-Inf"] <- 0
meanT3s_mesh[meanT3s_mesh == "Inf"] <- 0
meanT3s_mesh[meanT3s_mesh == "NaN"] <- 0
dim(meanT3s_mesh)
#meanT3s_mesh <- meanT3s_mesh[BFPOTUs.T3snet.sig,]
head(meanT3s_mesh)
dim(meanT3s_mesh)

# keep OTUs with significant contrasts
meanT3s_mesh <- meanT3s_mesh[T3s.sig.OTUs,]


# remove RR with pval >0.05
dim(ctrst.glm.CBFP.T3s.sub)
head(ctrst.glm.CBFP.T3s.sub)
# just to put the columns in the right order (same as the contrast df):
# ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[,c(4,5,1,3,6,2)]
ctrst.glm.CBFP.T3s.sub <- ctrst.glm.CBFP.T3s.sub[,c(4,1,3,5,2)]
# take only columns with RRs
RRT3s_mesh <- meanT3s_mesh[,7:11]
head(RRT3s_mesh)
# replace RR by zero when contrast is NS
RRT3s_mesh[ctrst.glm.CBFP.T3s.sub == 1]= 0
# check if some OTUs have only zeros
RRT3s_mesh <- RRT3s_mesh[rowSums(RRT3s_mesh) != 0,]

head(RRT3s_mesh)
# RRT3s_mesh$group <- str_split_fixed(rownames(RRT3s_mesh), "_", 2)[,1]
# RRT3s_mesh$OTU_names <- str_split_fixed(rownames(RRT3s_mesh), "_", 2)[,2]
# save(RRT3s_mesh, file="/home/mathilde/Bureau/STAT/R-objects/RRT3s_mesh_sig.RData") # to use for Itol
# load("/home/mathilde/Bureau/STAT/R-objects/RRT3s_mesh_sig.RData")





# save(meanT3s_mesh, file="/home/mathilde/Bureau/STAT/R-objects/meanT3s_mesh.RData")
# load("/home/mathilde/Bureau/STAT/R-objects/meanT3s_mesh.RData")
# nrow(meanT3s_mesh)


tmp_otus_T3s <- data.frame(str_split_fixed(rownames(RRT3s_mesh), "_", 2))
colnames(tmp_otus_T3s)<- c("OTU_names")


RRT3s_mesh <- cbind(tmp_otus_T3s, RRT3s_mesh)
head(RRT3s_mesh)


# nrow(RRT3s_mesh[RRT3s_mesh$group == "B",3:ncol(RRT3s_mesh)]) # check column numbers
# nrow(RRT3s_mesh[RRT3s_mesh$group == "P",3:ncol(RRT3s_mesh)])




RRT3s_meshB <- RRT3s_mesh
# rownames(RRT3s_meshB) <- str_split_fixed(rownames(RRT3s_meshB),"_", 2)[,2]
nrow(RRT3s_meshB)
head(RRT3s_meshB)

# #reorder column
# RRT3s_meshB<-RRT3s_meshB[,c(1,2,6,3,4,5,7)]

# Heatmaps --------------
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(as.matrix(RRT3s_meshB[,3:7]), cluster_columns = FALSE, col= col_fun)
# Heatmap(as.matrix(RRT3s_mesh[RRT3s_mesh$group == "P",3:8]), cluster_columns = FALSE, col= col_fun)
# Heatmap(as.matrix(RRT3s_mesh[RRT3s_mesh$group == "F",3:8]), cluster_columns = FALSE, col= col_fun)



# make a column with the type of response (positive or negative)



meanT3s_mesh <- cbind(tmp_otus_T3s, meanT3s_mesh)
head(meanT3s_mesh)
head(RRT3s_mesh)
# replace the RRs by the the ones already sorted out by significance, check column number
# meanT3s_mesh<-meanT3s_mesh[,c(1,2,7,3,4,6,8,5,12,9,10,11,13)]
meanT3s_mesh[,9:13] <- RRT3s_mesh[,3:7]


meanT3s_mesh0 <- meanT3s_mesh
meanT3s_mesh0$OTUs <- rownames(meanT3s_mesh0)
head(meanT3s_mesh0)
meanT3s_mesh0$RR_Col0 <- NULL
#meanT3s_mesh0$OTU_names <- NULL
#meanT3s_mesh0$X40000 <- NULL
# keep the columns of the means, make a column with RR types as variable, followed value column (of RRs)

meanT3s_mesh0<-meanT3s_mesh0[,-2]

meanT3s_mesh0.melt <- melt(meanT3s_mesh0, id.vars=c("OTU_names", "No", "Col0", "Ahb1", "Nia1Nia2", "Nox1", "GSNOR1", "OTUs"))
# meanT3s_mesh0.melt <- melt(meanT3s_mesh0, id.vars=c("OTU_names", "Ahb1", "Col0", "GSNOR1", "Nia1Nia2", "No", "Nox1","OTUs"))

head(meanT3s_mesh0.melt)
meanT3s_mesh0.melt$response <- c("")
meanT3s_mesh0.melt[meanT3s_mesh0.melt$value > 0 ,]$response <- c("positive")
meanT3s_mesh0.melt[meanT3s_mesh0.melt$value < 0 ,]$response <- c("negative")
meanT3s_mesh0.melt[meanT3s_mesh0.melt$value == 0 ,]$response <- c("NS")

meanT3s_mesh0.melt$response
class(meanT3s_mesh0.melt$response)
head(meanT3s_mesh0.melt)


######################################################################################
#### 
######################################################################################

# Nox1<0
meanT3s_mesh[which(meanT3s_mesh[,12]<0),]

# Nox1>0
meanT3s_mesh[which(meanT3s_mesh[,12]>0),]




# table(meanT3s_mesh0.melt[!duplicated(meanT3s_mesh0.melt$OTU_names),c(1,13)])
# table(meanT3s_mesh0.melt[,c(1,13)])



# cast the responses to get one category of response for each OTU (ie majority of neg = overall negative)----

repT3s_cast <- data.frame(dcast(meanT3s_mesh0.melt, OTUs ~ variable, value.var = "response"))
head(repT3s_cast)
rownames(repT3s_cast)<- repT3s_cast$OTUs
repT3s_cast$OTUs <- NULL


counts_repT3s <- data.frame(count_positive=rowSums(repT3s_cast == "positive"), count_negative= rowSums(repT3s_cast=="negative"))
counts_repT3s$overall <- c("")
counts_repT3s[counts_repT3s$count_positive != 0 & counts_repT3s$count_negative == 0,]$overall <- c("positive")
counts_repT3s[counts_repT3s$count_positive == 0 & counts_repT3s$count_negative != 0,]$overall <- c("negative")
counts_repT3s[counts_repT3s$count_positive != 0 & counts_repT3s$count_negative != 0,]$overall <- c("both")
counts_repT3s[counts_repT3s$count_positive == 0 & counts_repT3s$count_negative == 0,]$overall <- c("NS")

# save(counts_repT3s, file="/home/mathilde/Bureau/STAT/R-objects/counts_repT3s.RData")

head(counts_repT3s)
(table(counts_repT3s[,c(3,6)]))


# mean of all OTU abundance (rel abund within treatment not across all samples) with negative or positive response across mesh sizes ----
# rel abund calculation ---------
head(tmp1T3s)
colnames(tmp1T3s)[1]<- c("OTUs")
head(meanT3s_mesh0.melt) # insignificant contrasts already sorted out


mean_responseT3s <- (left_join(meanT3s_mesh0.melt[,c(1,8:11)], tmp1T3s[,c(1:7)], by=c("OTUs")))

head((mean_responseT3s[1:100,]))
mean_responseT3s$mesh_sign <- paste(mean_responseT3s$variable, mean_responseT3s$response, sep="-")
mean_responseT3s$mesh_rr <- str_split_fixed(mean_responseT3s$variable, "_", 2)[,2]
mean_responseT3s_match <- mean_responseT3s[mean_responseT3s$mesh_rr == mean_responseT3s$mesh_size_um,]
nrow(mean_responseT3s)/nrow(mean_responseT3s_match)

str(mean_responseT3s_match )
mean_by_mesh_responseT3s <- as.data.frame(mean_responseT3s_match %>% group_by(gp_sum,gp_OTUsum, mesh_sign) %>% summarise(sum_bymesh=sum(Abundance)))
head(mean_by_mesh_responseT3s )
mean_by_mesh_responseT3s$rel_abund <- mean_by_mesh_responseT3s$sum_bymesh*100/mean_by_mesh_responseT3s$gp_sum

# for the plot regrouping all meshes, rel abund needs to be calculated with total = total of all mesh, -----
#instead of suming the rel abund within meshes
CBFPT3s.mesh.sum

head(mean_by_mesh_responseT3s)
mean_by_mesh_responseT3s$group_total <- c("")
mean_by_mesh_responseT3s$group_total <- sum(CBFPT3s.mesh.sum$B)
# mean_by_mesh_responseT3s[mean_by_mesh_responseT3s$group == "F",]$group_total <- sum(CBFPT3s.mesh.sum$F)
# mean_by_mesh_responseT3s[mean_by_mesh_responseT3s$group == "P",]$group_total <- sum(CBFPT3s.mesh.sum$P)
mean_by_mesh_responseT3s$group_total <- as.numeric(mean_by_mesh_responseT3s$group_total)
mean_by_mesh_responseT3s$rel_abund_total <- mean_by_mesh_responseT3s$sum_bymesh*100/mean_by_mesh_responseT3s$group_total









# make a table of abundance and number of decreasing and increasing OTUs ----

nrow(meanT3s_mesh0)
head(meanT3s_mesh0)
head(RRT3s_mesh)

# an idea to make it cleaner :
#https://www.r-bloggers.com/2012/06/transforming-subsets-of-data-in-r-with-by-ddply-and-data-table/
RRT3s_meshB<- RRT3s_mesh
head(RRT3s_meshB)

BT3s.shift.counts <- data.frame(RR_No = c(nrow(RRT3s_mesh[RRT3s_mesh$RR_No > 0,]), nrow(RRT3s_mesh[RRT3s_mesh$RR_No < 0,])),
                                RR_Ahb1 = c(nrow(RRT3s_mesh[RRT3s_mesh$RR_Ahb1 > 0,]), nrow(RRT3s_mesh[RRT3s_mesh$RR_Ahb1 < 0,])),
                                RR_Nia1Nia2 = c(nrow(RRT3s_mesh[RRT3s_mesh$RR_Nia1Nia2 > 0,]), nrow(RRT3s_mesh[RRT3s_mesh$RR_Nia1Nia2 < 0,])),
                                RR_Nox1 = c(nrow(RRT3s_mesh[RRT3s_mesh$RR_Nox1 > 0,]), nrow(RRT3s_mesh[RRT3s_mesh$RR_Nox1 < 0,])),
                                RR_GSNOR1 = c(nrow(RRT3s_mesh[RRT3s_mesh$RR_GSNOR1 > 0,]), nrow(RRT3s_mesh[RRT3s_mesh$RR_GSNOR1 < 0,])),
                                row.names = c("B-positive", "B-negative"))


# RRT3s_meshF <- RRT3s_mesh[RRT3s_mesh$group == "F",]
# head(RRT3s_meshF)


# FT3s.shift.counts <- data.frame(RR_31 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_31 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_31 < 0,])),
                                # RR_50 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_50 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_50 < 0,])),
                                # RR_100 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_100 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_100 < 0,])),
                                # RR_250 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_250 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_250 < 0,])),
                                # RR_500 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_500 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_500 < 0,])),
                                # RR_1000 = c(nrow(RRT3s_meshF[RRT3s_meshF$RR_1000 > 0,]), nrow(RRT3s_meshF[RRT3s_meshF$RR_1000 < 0,])),
                                # row.names = c("F-positive", "F-negative"))

# RRT3s_meshP <- RRT3s_mesh[RRT3s_mesh$group == "P",]
# head(RRT3s_meshP)


# PT3s.shift.counts <- data.frame(RR_31 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_31 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_31 < 0,])),
                                # RR_50 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_50 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_50 < 0,])),
                                # RR_100 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_100 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_100 < 0,])),
                                # RR_250 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_250 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_250 < 0,])),
                                # RR_500 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_500 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_500 < 0,])),
                                # RR_1000 = c(nrow(RRT3s_meshP[RRT3s_meshP$RR_1000 > 0,]), nrow(RRT3s_meshP[RRT3s_meshP$RR_1000 < 0,])),
                                # row.names = c("P-positive", "P-negative"))

# BFPT3s.shift.counts <- rbind(BT3s.shift.counts, FT3s.shift.counts, PT3s.shift.counts)

BFPT3s.shift.counts<-BT3s.shift.counts

BFPT3s.shift.counts$reponse <- str_split_fixed(rownames(BFPT3s.shift.counts), "-", 2)[,2]
# BFPT3s.shift.counts$group <- str_split_fixed(rownames(BFPT3s.shift.counts), "-", 2)[,1]
BFPT3s.shift.counts.melt <- melt(BFPT3s.shift.counts)
colnames(BFPT3s.shift.counts.melt)[2]<- c("RR_type")
colnames(BFPT3s.shift.counts.melt)[1]<- c("response")
head(BFPT3s.shift.counts.melt)







mean_by_mesh_responseT3s$RR_type <- str_split_fixed(mean_by_mesh_responseT3s$mesh_sign, "-", 2)[,1]
mean_by_mesh_responseT3s$response <- str_split_fixed(mean_by_mesh_responseT3s$mesh_sign, "-", 2)[,2]
head(mean_by_mesh_responseT3s)

abund_counts_responseT3s <- left_join(BFPT3s.shift.counts.melt, mean_by_mesh_responseT3s[,c(1:9)], by=c("RR_type", "response"))

head(abund_counts_responseT3s)
colnames(abund_counts_responseT3s)[3]<- c("OTU_number")

abund_counts_responseT3s$rel_otu_nb_total <- abund_counts_responseT3s$OTU_number*100/abund_counts_responseT3s$gp_OTUsum
unique(abund_counts_responseT3s$gp_OTUsum)
head(abund_counts_responseT3s)

abund_counts_responseT3s$gp_OTUsum_all <- abund_counts_responseT3s$gp_OTUsum
# # change numbers accordingly
# abund_counts_responseT3s$gp_OTUsum_all <- gsub(703, 703, abund_counts_responseT3s$gp_OTUsum_all)
# abund_counts_responseT3s$gp_OTUsum_all <- gsub(71, 71, abund_counts_responseT3s$gp_OTUsum_all)
# abund_counts_responseT3s$gp_OTUsum_all <- gsub(231, 231, abund_counts_responseT3s$gp_OTUsum_all)
# abund_counts_responseT3s$gp_OTUsum_all <- as.numeric(abund_counts_responseT3s$gp_OTUsum_all )
# abund_counts_responseT3s$rel_otu_nb_total_all <- abund_counts_responseT3s$OTU_number*100/abund_counts_responseT3s$gp_OTUsum_all

abund_counts_responseT3s.melt <- melt(abund_counts_responseT3s)

str(abund_counts_responseT3s.melt)
abund_counts_responseT3s.melt$response <- as.factor(abund_counts_responseT3s.melt$response)
# abund_counts_responseT3s.melt$group <- as.factor(abund_counts_responseT3s.melt$group)
# abund_counts_responseT3s.melt$RR <- factor(abund_counts_responseT3s.melt$RR, levels=c("RR_1000", "RR_500", "RR_250", "RR_100", "RR_50", "RR_31"))
abund_counts_responseT3s.melt$RR <- factor(abund_counts_responseT3s.melt$RR, levels=c("RR_No","RR_Ahb1", "RR_Nia1Nia2", "RR_Nox1", "RR_GSNOR1"))

str(abund_counts_responseT3s.melt)

ggplot(abund_counts_responseT3s.melt, aes(x = RR, y = value)) + geom_bar(stat = "identity")+ facet_wrap(~variable * response, ncol=3, scales = "free")

abund_counts_responseT3s.melt0<- abund_counts_responseT3s.melt
abund_counts_responseT3s.melt0[abund_counts_responseT3s.melt0$response == "negative",]$value <- abund_counts_responseT3s.melt0[abund_counts_responseT3s.melt0$response == "negative",]$value * -1
abund_counts_responseT3s.melt0$value
head(abund_counts_responseT3s.melt0)


# plots------

# abund_counts_responseT3s.melt00 <- abund_counts_responseT3s.melt0[abund_counts_responseT3s.melt0$variable == "rel_abund" | 
                                                                    # abund_counts_responseT3s.melt0$variable == "rel_otu_nb_total",]
abund_counts_responseT3s.melt00 <- abund_counts_responseT3s.melt0[abund_counts_responseT3s.melt0$variable == "rel_otu_nb_total",]

abund_counts_responseT3s.melt00$RR <- gsub("RR_", "", abund_counts_responseT3s.melt00$RR)
abund_counts_responseT3s.melt00$RR <- factor(abund_counts_responseT3s.melt00$RR, levels = c("No","Ahb1", "Nia1Nia2", "Nox1", "GSNOR1"))
# abund_counts_responseT3s.melt00$group <- gsub("B", "Bacteria", abund_counts_responseT3s.melt00$group)
# abund_counts_responseT3s.melt00$group <- gsub("F", "Fungi", abund_counts_responseT3s.melt00$group)
# abund_counts_responseT3s.melt00$group <- gsub("P", "Fauna", abund_counts_responseT3s.melt00$group)
# abund_counts_responseT3s.melt00$group <- factor(abund_counts_responseT3s.melt00$group, levels=c("Bacteria", "Fungi", "Fauna"))

abund_counts_responseT3s.p <- ggplot(abund_counts_responseT3s.melt00, aes(x = RR, y = value, fill = response)) + 
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  theme_bw() + theme(axis.text.x = element_text(angle=90, size = 12, hjust=1, vjust=0.5),
                                                  strip.text.x= element_text(size = 16),
                                                  axis.text.y = element_text(size = 12),
                                                  axis.title.y = element_text(size=16),
                                                  axis.title.x = element_blank(),
                                                  panel.grid = element_blank())+ylab("affected OTUs (%)")

abund_counts_responseT3s.p




head(counts_repT3s)
counts_repT3s$OTUs <- rownames(counts_repT3s)

tmp1T3s$OTUs <- gsub("-", "\\.",tmp1T3s$OTU)

#View(tmp1T3s)
head(tmp1T3s)
tmp1T3s$group_total <- c("")
tmp1T3s$group_total <- sum(CBFPT3s.mesh.sum$B)
# tmp1T3s[tmp1T3s$group == "F",]$group_total <- sum(CBFPT3s.mesh.sum$F)
# tmp1T3s[tmp1T3s$group == "P",]$group_total <- sum(CBFPT3s.mesh.sum$P)

CBFPT3s.mesh.OTUsum
tmp1T3s$group_OTUnb <- c("")
tmp1T3s$group_OTUnb <- max(CBFPT3s.mesh.OTUsum$B)
# tmp1T3s[tmp1T3s$group == "F",]$group_OTUnb <- max(CBFPT3s.mesh.OTUsum$F)
# tmp1T3s[tmp1T3s$group == "P",]$group_OTUnb <- max(CBFPT3s.mesh.OTUsum$P)


tmp1T3s$pa <- tmp1T3s$Abundance
tmp1T3s[tmp1T3s$pa> 0,]$pa <- 1



rep_all_T3s <- left_join(counts_repT3s, tmp1T3s)
head(rep_all_T3s)

rep_all_T3s_sum <- as.data.frame(rep_all_T3s %>% group_by(overall,group_total, group_OTUnb) %>% 
                                   summarise(sum_byresponse=sum(Abundance), prevalence_byresp=length(unique(OTUs))))

rep_all_T3s_sum$group_total <- as.numeric(rep_all_T3s_sum$group_total)
rep_all_T3s_sum$group_OTUnb <- as.numeric(rep_all_T3s_sum$group_OTUnb)


rep_all_T3s_sum$total_rel_abund <- rep_all_T3s_sum$sum_byresponse*100/rep_all_T3s_sum$group_total
rep_all_T3s_sum$total_prevalence <- rep_all_T3s_sum$prevalence_byresp*100/rep_all_T3s_sum$group_OTUnb

rep_all_T3s_sum.melt <- melt(rep_all_T3s_sum[,c(1,2,6,7)])
head(rep_all_T3s_sum.melt)


overall_resp_T3s <- ggplot(subset(rep_all_T3s_sum.melt,rep_all_T3s_sum.melt$variable=="total_prevalence"), aes(x = overall, y = value, fill = overall)) + geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ 
  theme_bw() + theme(axis.text.x = element_text(angle=90, size = 10, hjust=1),
                                                 strip.text.x= element_text(size = 11),
                                                 axis.text.y = element_text(size = 10),
                                                 axis.title.y = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 panel.grid = element_blank())+
  labs(title = "")+ facet_wrap(~ variable)#+ scale_y_continuous(limits = c(-600, 250))


overall_resp_T3s
overall_resp_T1

ggarrange(overall_resp_T3s, abund_counts_responseT3s.p)
library(ggpubr)

ggarrange(overall_resp_T1, overall_resp_T2, overall_resp_T3s, nrow=1, legend = "bottom", common.legend = T)

ggarrange(abund_counts_responseT1.p, abund_counts_responseT2.p, abund_counts_responseT3sm.p, nrow=1, legend = "bottom", common.legend = T)




head(RRT3s_mesh)

RRT3s_mesh.melt <- melt(RRT3s_mesh[,c(1,3:7)])
head(RRT3s_mesh.melt)
RRT3s_mesh.melt$sign <- c("")
RRT3s_mesh.melt[RRT3s_mesh.melt$value < 0,]$sign <- c("negative")
RRT3s_mesh.melt[RRT3s_mesh.melt$value > 0,]$sign <- c("positive")

RRT3s_mesh.melt$sign_RR <- paste(RRT3s_mesh.melt$sign, RRT3s_mesh.melt$variable, sep="_")

RRT3s_mesh.melt <- RRT3s_mesh.melt[RRT3s_mesh.melt$value != 0,]
head(RRT3s_mesh.melt)
RRT3s_mesh.melt$variable <- gsub("RR_", "", RRT3s_mesh.melt$variable)
RRT3s_mesh.melt$variable <- factor(RRT3s_mesh.melt$variable, levels = c("No", "Ahb1", "Nia1Nia2", "Nox1", "GSNOR1"))
# RRT3s_mesh.melt$group <- gsub("B", "Bacteria", RRT3s_mesh.melt$group)
# RRT3s_mesh.melt$group <- gsub("F", "Fungi", RRT3s_mesh.melt$group)
# RRT3s_mesh.melt$group <- gsub("P", "Fauna", RRT3s_mesh.melt$group)
# RRT3s_mesh.melt$group <- factor(RRT3s_mesh.melt$group, levels = c("Bacteria", "Fungi", "Fauna"))

RRT3s.p= ggplot(RRT3s_mesh.melt, aes(x = variable, y = value,group=sign_RR, color=sign)) + 
  geom_point(size = 2, alpha = 0.5)+
  scale_color_manual(values = c("#f51853","#2166AC"))+
  theme_bw()  + ylab("log2 fold change")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5,size = 0.3, color="black")+
  scale_y_continuous(breaks =2*-10:10, limits = c(-10,10))+geom_hline(yintercept=0, color="black")+
  scale_color_manual(values=c("#FF6666","#6666FF"))+
  theme(axis.text.x = element_text(angle=90, size = 12, hjust=1, vjust=0.5),
        strip.text.x= element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        panel.grid = element_blank())+ylab("Log2-fold change")+xlab("mesh size")
RRT3s.p


#RRT3s.p= ggplot(RRT3s_mesh.melt, aes(x = variable, y = value,group=sign_RR, color=sign)) + 
#  geom_jitter(position = position_jitter(width = .20), alpha = 0.7, size = 1.5) + theme_bw() +
#  geom_boxplot(colour = "black", lwd=0.3, alpha=0.1, outlier.colour = NA, position = position_dodge(width = 0)) +
#  xlab("")+ylab("")+ylim(-10,10)+
#  ggtitle("")+
#  theme(plot.title = element_text(size = rel(1.2)),
#        axis.text.x = element_text(angle=90, size = 10, hjust=1, vjust=0.5),
#        strip.text.x= element_text(size = 13),
#        axis.text.y = element_text(size = 10),
#        axis.title.y = element_text(size=12))+facet_wrap(~ group)+scale_color_manual(values=c("#FF6666","#6666FF"))
#RRT3s.p

ggarrange(abund_counts_responseT3s.p, RRT3s.p, ncol=1, common.legend = T, legend = "bottom")


# test quickly RR difference between positive and negative response ------------
shapiro.test(RRT3s_mesh.melt$value_abs)
RRT3s_mesh.melt$value_abs <- abs(RRT3s_mesh.melt$value)
wilcox.test(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Bacteria",]$value_abs ~ RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Bacteria",]$sign)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Bacteria"&RRT3s_mesh.melt$sign== "negative",]$value_abs)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Bacteria"&RRT3s_mesh.melt$sign== "positive",]$value_abs)
wilcox.test(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fungi",]$value_abs ~ RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fungi",]$sign)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fungi"&RRT3s_mesh.melt$sign== "negative",]$value_abs)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fungi"&RRT3s_mesh.melt$sign== "positive",]$value_abs)
wilcox.test(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna",]$value_abs ~ RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna",]$sign)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna"&RRT3s_mesh.melt$sign== "negative",]$value_abs)
mean(RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna"&RRT3s_mesh.melt$sign== "positive",]$value_abs)


dunn_test(value_abs ~ variable, data=RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Bacteria"&RRT3s_mesh.melt$sign== "negative",])
dunn_test(value_abs ~ variable, data=RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fungi",])
dunn_test(value_abs ~ variable, data=RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna"&RRT3s_mesh.melt$sign== "positive",])

#test gradual increase across sizes ----
PposRRT3s_mesh.melt <- RRT3s_mesh.melt[RRT3s_mesh.melt$group == "Fauna"&RRT3s_mesh.melt$sign== "positive",]
ad.test((PposRRT3s_mesh.melt$value))
dunn_test(value ~ sign_RR, data=PposRRT3s_mesh.melt)
