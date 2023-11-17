###################################################################
#######################Unifrac distance between treatments

BCstep1 = distance(step1.rarefied,  "wunifrac")
BCstep1
d <- as.dist(BCstep1)
#BCstep1 <- as.matrix(BCstep1)
#write.table(BCstep1, file="weighted UniFrac distance matrix.csv")

item_groups <- sample_data(step1.rarefied)
item_groups <- item_groups$Treatment
#calculate dist between groups
d.calcul <- dist_groups(d, item_groups)


#Control
tab.distance = as_tibble(d.calcul) 
tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)


tab.distance.C <- subset(tab.distance, Label%in% c("BetweenCd1andCd2","BetweenCd1andCd3", "BetweenCd1andCd4","BetweenCd2andCd3","BetweenCd2andCd4","BetweenCd3andCd4"))



kruskal.test(Distance ~ Label,
             data = tab.distance.C)

### NemenyiTes
tab.distance.C$Label <- as.factor(tab.distance.C$Label)
PT = PMCMRplus::kwAllPairsNemenyiTest(x=tab.distance.C$Distance, g=tab.distance.C$Label,
                                      dist="Tukey")

PT = PT$p.value
PT
PT1 = fullPTable(PT)

Tukey <- multcompLetters(PT1,  
                         compare="<",  
                         threshold=0.05,
                         Letters=letters,  
                         reversed = FALSE)
Tukey <- as.data.frame(Tukey$Letters)
colnames(Tukey) <- c("Letters")
Tukey$Label <- rownames(Tukey)

### Tuckey Groups

treatOrderC <- factor(c("BetweenCd1andCd2","BetweenCd1andCd3", "BetweenCd1andCd4","BetweenCd2andCd3","BetweenCd2andCd4","BetweenCd3andCd4"))

tuckeyGroups <- data.frame(Label = Tukey$Label,Groups = Tukey$Letters)
tuckeyGroups <- tuckeyGroups[match(treatOrderC , tuckeyGroups$Label),]
tuckeyGroups$Label <- reorder.factor(tuckeyGroups$Label, new.order = treatOrderC )

sumData <- ddply(tab.distance.C, "Label", summarise,
                 N    = length(Distance),
                 Mean = mean(Distance),
                 Sd   = sd(Distance),
                 Se   = Sd / sqrt(N)
)



CONTROL<- ggplot(sumData, aes(x = Label, y = Mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_text(data=tuckeyGroups,aes(x=Label, y = sumData$Mean + sumData$Sd, label=Groups), vjust=-1) +
  
  ylab("Weighted Unifrac Distances") +
  theme_bw() + 
  
  theme(plot.title = element_text(hjust = 0.5))

################################################################################################################################################
###################################################################
#######################Unifrac distance between treatments
#devtools::install_github("kylebittinger/usedist")
#install.packages("PMCMRplus")
#install.packages("gdata")
#install.packages("PMCMR")
library(PMCMR)
library(usedist)
library(PMCMRplus)
library(gdata)
#BCstep1 = distance(aoa.physeq_bulk1,  "wunifrac")
#BCstep1
#d <- as.dist(BCstep1)
#BCstep1 <- as.matrix(BCstep1)
#write.table(BCstep1, file="weighted UniFrac distance matrix.csv")

item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.bulk_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)

#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)


#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))


kruskal.test(Distance ~ Label,
             data = tab.distance.C)
str(tab.distance.C)
tab.distance.C


dunn=dunnTest(Distance ~ Label,
         data = tab.distance.C, method = "bh")
library(rcompanion)
CLD = cldList(P.adj ~ Comparison, data=dunn$res)
rownames(CLD) = sumData$Label
CLD = rownames_to_column(CLD, var="Label")





### NemenyiTes
tab.distance$Label <- as.factor(tab.distance$Label)
PT = PMCMRplus::kwAllPairsNemenyiTest(x=tab.distance.C$Distance, g=tab.distance.C$Label,
                                      dist="Tukey")

PT = PT$p.value
PT
PT1 = fullPTable(PT)

Tukey <- multcompLetters(PT1,  
                         compare="<",  
                         threshold=0.05,
                         Letters=letters,  
                         reversed = FALSE)
Tukey <- as.data.frame(Tukey$Letters)
colnames(Tukey) <- c("Letters")
Tukey$Label <- rownames(Tukey)
Tukey

### Tuckey Groups

treatOrderC <- factor(c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))

tuckeyGroups <- data.frame(Label = Tukey$Label,Groups = Tukey$Letters)
tuckeyGroups <- tuckeyGroups[match(treatOrderC , tuckeyGroups$Label),]
tuckeyGroups$Label <- reorder.factor(tuckeyGroups$Label, new.order = treatOrderC )

sumData <- ddply(tab.distance.C, "Label", summarise,
                 N    = length(Distance),
                 Mean = mean(Distance),
                 Sd   = sd(Distance),
                 Se   = Sd / sqrt(N)
)



CONTROL<- ggplot(sumData, aes(x = Label, y = Mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           linewidth=.3) +      # Thinner lines
  geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se),
                linewidth=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_text(data=CLD,aes(x=Label, y = sumData$Mean + sumData$Sd, label=Letter), vjust=0) +
  
  ylab("Bray-Curtis Distances") +
  theme_bw() + 
  
  theme(plot.title = element_text(hjust = 0.5))
CONTROL


library(reshape2)

# wrangle distance matrix into a longer dataframe
tmp_dist_matrix = melt(as.matrix(aoa.bulk_dist_bc))
# remove self-comparisons
tmp_dist_matrix = tmp_dist_matrix[tmp_dist_matrix$Var1 != tmp_dist_matrix$Var2,]
# select sample data
tmp_sam_data = tibble("sample"=rownames(aoa.meta.bulk),
                      "treatment"=aoa.meta.bulk$x)
# combined distance matrix with sample data
colnames(tmp_sam_data) = c("Var1", "treatment1")
tmp_data <- left_join(tmp_dist_matrix, tmp_sam_data, by = "Var1")
colnames(tmp_sam_data) = c("Var2", "treatment2")
tmp_data <- left_join(tmp_data, tmp_sam_data, by = "Var2")
# select distances from the control
tmp_data <- tmp_data[tmp_data$treatment1 == "",]


library(agricolae)

# linear model
tmp_lm = lm(tmp_data$value ~ tmp_data$treatment2, na.action = na.omit )
summary(tmp_lm)
# anova
anova(tmp_lm)
tmp_aov = aov(tmp_lm)
summary(tmp_aov)
# post hoc Tukey test
tmp_comp <- HSD.test(tmp_aov,'tmp_data$treatment',alpha = 0.05,group = T)


tmp_lm = kruskal.test(tmp_data$value ~ tmp_data$treatment2, na.action = na.omit)
tmp_lm

dunnTest(tmp_data$value ~ as.factor(tmp_data$treatment2), method = "bh")





# tibble with statistical groups
tmp_stat = tibble("treatment"=rownames(tmp_comp[["groups"]]),
                  "mean"=tmp_comp[["groups"]][["tmp_data$value"]],
                  "stat_groups"=tmp_comp[["groups"]][["groups"]])

tmp_stat






###############################
library(usedist)
dist <- aoa.bulk_dist_bc

aoa.physeq_bulk_D <- subset_samples(aoa.physeq_bulk1, Treatment=="D")
aoa.physeq_bulk_D1 <- prune_taxa(taxa_sums(aoa.physeq_bulk_D)>0, aoa.physeq_bulk_D)
sort(taxa_sums(aoa.physeq_bulk_D1), decreasing =F)
aoa.bulk.D <- as.data.frame(otu_table(aoa.physeq_bulk_D1))

aoa_bulk_D_dist <- vegdist(t(aoa.bulk.D), method = "bray")  
  
  

samples<- c(1:40)
control <- c(1,3,5,7,9,11,13,17,19,21,23,25,27,29,31,33,35,37,39)


results<-c("")
for (i in seq_along(samples)) {
  results[i] <- dist_between_centroids(d = aoa_bulk_D_dist, idx1 = samples[i], 
                                       idx2 = control)
}
results

aoa.t.bulk.D <- sample_data(aoa.physeq_bulk_D1)
aoa.t.bulk.D.irri <- aoa.t.bulk.D$Irrigation
str(aoa.t.bulk.D.irri)

tmp_kw = kruskal.test(as.numeric(results) ~ as.factor(aoa.t.bulk.M.irri), na.action = na.omit)
tmp_kw
result.t <- as.numeric(results)

c.d <- result.t[c(1,3,5,7,9,11,13,17,19,21,23,25,27,29,31,33,35,37,39)]
d <- result.t[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

res <- wilcox.test(c.d, d)
res


dunnTest(result.t ~ as.factor(aoa.t.bulk.M.irri), method = "bh")
wilcox.test(result.t ~ aoa.t.bulk.M.irri)


aoa.physeq_bulk_K <- subset_samples(aoa.physeq_bulk1, Treatment=="K")
aoa.physeq_bulk_K1 <- prune_taxa(taxa_sums(aoa.physeq_bulk_K)>0, aoa.physeq_bulk_K)
sort(taxa_sums(aoa.physeq_bulk_K1), decreasing =F)
aoa.bulk.K <- as.data.frame(otu_table(aoa.physeq_bulk_K1))

aoa_bulk_K_dist <- vegdist(t(aoa.bulk.K), method = "bray") 

################################################################################
#microeco
library(microeco)
# Preparing the microtable class

# rarefied AOA ASV table of bulk soil
aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS")
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# taxonomy table AOA
aoa.asv.tax <- as.data.frame(tax_table(aoa.physeq_bulk1))
dim(aoa.asv.tax)
aoa.asv.tax$Group <- 'AOA'#adding a new column with AOA as value
# enviromental data
str(aoa.meta.bulk)
aoa.asv.env <- aoa.meta.bulk %>% 
  dplyr::rename("AOA_Richness" = "Richness",
                "AOA_Shannon" = "Shannon",
                "AOA_InvSimpson" = "InvSimpson") # change column names
aoa.asv.env <- aoa.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','AOA_Richness'), as.numeric)

# create  a microtable
aoa.microdata <- microtable$new(sample_table = aoa.asv.env, otu_table = aoa.asv.tab, tax_table = aoa.asv.tax)
# calculate beta diversity
aoa.microdata$tidy_dataset()
aoa.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = aoa.microdata, group = "x", measure = "bray")
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean")  
# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "KW_dunn")
t1$plot_group_distance(boxplot_add = "mean")

