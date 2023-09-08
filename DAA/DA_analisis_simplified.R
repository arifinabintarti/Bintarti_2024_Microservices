############################################################################################################################################################################################ 
# DIFFERENTIAL ABUNDANCE - ENDO
# 
################################
# Filter low-abundant taxa ENDO
###############################
physeq.subset <- endo

data.obs <- as.data.frame(otu_table(physeq.subset))

### keeping OTUs with at least 0.02 % relative abundance across all samples
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.00005)
data.F=data.obs[keep.taxa.id,,drop=FALSE]

new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incroporate into phyloseq Object

physeq.subset


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
dim(prevalencedf)


### calculate prevalence /!\ takes from 30min to 3h /!\

ps = physeq.subset 

df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1

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
o = log(sample_sums(endo))
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
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson', offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    #glmT3s.sum = summary(glmT3s)$coefficients$cond
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
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

# replace pvalues to 0 if non significant, or 1 if significant
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >0.05] <- 2
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub <0.05] <- 1
ctrst.glm.CBFP.T3s.sub[ctrst.glm.CBFP.T3s.sub >1] <- 0
head(ctrst.glm.CBFP.T3s.sub)

# Calculate the OTU avg per treatment
# CHECK THE OBJECT
meanotus<-phyloseq_average(endo,avg_type="arithmetic",acomp_zero_impute = NULL,group="Genotype")
meanotus<-as.data.frame(otu_table(meanotus));meanotus
# same order for both meanotus and tmp_otu3s
meanotus<-meanotus[tmp_otu3s,]

# Calculate log2fold ratios for all OTUs in the filtered table
meanotus$RR_Ahb1 <- log2(meanotus$Ahb1 / meanotus$Col0)
meanotus$RR_Nia1Nia2 <- log2(meanotus$Nia1Nia2 / meanotus$Col0)
meanotus$RR_Nox1 <- log2(meanotus$Nox1 / meanotus$Col0)
meanotus$RR_GSNOR1 <- log2(meanotus$GSNOR1 / meanotus$Col0)

head(meanotus)
# keep only columns containing log2fold ratios (RRs)
meanotus<-meanotus[,c(6:9)]
head(meanotus)

head(ctrst.glm.CBFP.T3s.sub)

# put the same column order in ctrst.glm.CBFP.T3s.sub and in meanotus
ctrst.glm.CBFP.T3s.sub<-ctrst.glm.CBFP.T3s.sub[,c(1,3,4,2)]
# replace "-" by "." to be able to compare both datasets. Also, put OTUs in the same order in both cases
row.names(meanotus)<-gsub("-", ".", row.names(meanotus))
meanotus<-meanotus[row.names(ctrst.glm.CBFP.T3s.sub),]
head(meanotus)

# Multiply the matrices to get the RR when it is significant and 0 when it is not significant
rr<-meanotus*ctrst.glm.CBFP.T3s.sub

# Initialize a list to store the results
rm(results_list,result)
results_list <- list()

# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr[[treatment]] > 0) / length(taxa_names(physeq.subset))
  decreased <- sum(rr[[treatment]] < 0) / length(taxa_names(physeq.subset))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}

# Combine the results into a single dataframe and sort the treatments according to the desired order
result <- do.call(rbind, results_list)
result$treatment<-factor(result$treatment, levels=c("RR_Ahb1","RR_GSNOR1","RR_Nia1Nia2","RR_Nox1"))

# Check that the object is in proper shape
print(result)

# Save the table for later
write.table(result,"result_endo.txt")
result<-read.table("result_endo.txt")

## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_responseT3s.p <- ggplot(result, aes(x = treatment, y = percentage*100, fill = response)) + 
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  theme_bw() + theme(axis.text.x = element_text(angle=90, size = 12, hjust=1, vjust=0.5),
                     strip.text.x= element_text(size = 16),
                     axis.text.y = element_text(size = 12),
                     axis.title.y = element_text(size=16),
                     axis.title.x = element_blank(),
                     panel.grid = element_blank())+ylab("affected OTUs (%)")+
  scale_y_continuous(limits = c(-35,20),breaks = c(-30, -20, -10, 0, 10, 20))


abund_counts_responseT3s.p
saveRDS(abund_counts_responseT3s.p,"abund_counts_responseT3s.p.endo.RDS")
svglite::svglite("abund_counts_responseT3s.p.endo.svg",fix_text_size = F,width = 5,height = 4)
abund_counts_responseT3s.p
dev.off()

# Make the product again for the graph involving log2fold ratios
rr<-meanotus*ctrst.glm.CBFP.T3s.sub
# write.table(rr,"rr_endo.txt")
# rr<-read.table("rr_endo.txt")
rr<-melt(rr)
rr <- rr %>% mutate(response = ifelse(value >= 0, "positive", "negative"))

# # Count the number of positive and negative values for each level
# rr<-rr %>%
#   group_by(variable) %>%
#   summarize(
#     positive_count = sum(value > 0),
#     negative_count = sum(value < 0)
#   )
# rr<-data.frame(rr)

#Remove RRs when 0
rr <- rr[rr$value != 0, ]

#Order
rr$variable<-factor(rr$variable, levels=c("RR_Ahb1","RR_GSNOR1","RR_Nia1Nia2","RR_Nox1"))

# Plot the log2fold ratios
RRT3s.p= ggplot(rr, aes(x = variable, y = value,group=response, color=response)) + 
  geom_point(size = 2, alpha = 0.5)+
  scale_color_manual(values = c("#f51853","#2166AC"))+
  theme_bw()  + ylab("log2 fold change")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5,size = 0.3, color="black")+
  scale_y_continuous(breaks =2*-10:10, limits = c(-7.5,7.5))+geom_hline(yintercept=0, color="black")+
  scale_color_manual(values=c("#FF6666","#6666FF"))+
  theme(axis.text.x = element_text(angle=90, size = 12, hjust=1, vjust=0.5),
        strip.text.x= element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        panel.grid = element_blank())+ylab("Log2-fold change")+xlab("mesh size")
RRT3s.p
# saveRDS(RRT3s.p,"RRT3s.p.endo.RDS")
RRT3s.p<-readRDS("RRT3s.p.endo.RDS")

svglite::svglite("RRT3s.p.endo.svg",fix_text_size = F,width = 5,height = 4)
RRT3s.p
dev.off()
