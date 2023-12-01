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
library(magrittr)
library(purrr)
#Distance matrix

# 1. AOA Bulk Soil - CAP Distance - CONFYM (K)
aoa.cap.bulk.dist
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.bulk.dist, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.K <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.K$Label <- factor(tab.distance.K$Label)
str(tab.distance.K) 
# One-way ANOVA
set.seed(13)
K.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.K)
summary(K.dist.cap.aov) # significant, p-val = 0.0036, f-val= 5.65
# Post-Hoc Test
K.dist.cap.tuk <- TukeyHSD(K.dist.cap.aov)
Tukey.K <- as.data.frame(K.dist.cap.tuk$Label)
Tukey.K <- rownames_to_column(Tukey.K, var = "Comparison")
colnames(Tukey.K)[5] <- "p.adj"
# Make the significance letter
CLD.k = cldList(p.adj ~ Comparison, data=Tukey.K)
CLD.k
# re-order
CLD.k.ed = CLD.k[c(3,1,2),] # re-order
CLD.k.ed
#_______________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.K)
str(tab.distance.K)
tab.distance.K
# Post Hoc Dunn Test
set.seed(1333)
dunn.K=dunnTest(Distance ~ Label,
         data = tab.distance.K, method = "bh")
library(rcompanion)
CLD.k = cldList(P.adj ~ Comparison, data=dunn.K$res)
CLD.k
#_______________________________________________________________________________
# Plot
sumData.K <- ddply(tab.distance.K, "Label", summarise,
                 Max = max(Distance),
                 N    = length(Distance),
                 Mean = mean(Distance),
                 Sd   = sd(Distance),
                 Se   = Sd / sqrt(N))
sumData.K
rownames(CLD.k.ed) = sumData.K$Label
CLD.k.ed = rownames_to_column(CLD.k.ed, var="Label")
CLD.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

dist.CAP.bulk.k<- ggplot(tab.distance.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                     labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=CLD.k.ed,aes(x=Label, y = sumData.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
              legend.title = element_text(size=13, face='bold'),
              plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x =element_text(size=14),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              #axis.title.y=element_text(size=13,face="bold"),
              axis.title.x=element_blank(),
              legend.text=element_text(size=13),
              legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.003", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
dist.CAP.bulk.k


# 2. AOA Bulk Soil - CAP Distance - CONMIN
aoa.cap.bulk.dist
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.bulk.dist, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.M <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.M 
# One-way ANOVA
set.seed(13)
M.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.M)
summary(M.dist.cap.aov) # not significant, p-val = 0.129, f-val= 2.055
# Post-Hoc Test
M.dist.cap.tuk <- TukeyHSD(M.dist.cap.aov)
Tukey.M <- as.data.frame(M.dist.cap.tuk$Label)
Tukey.M <- rownames_to_column(Tukey.M, var = "Comparison")
colnames(Tukey.M)[5] <- "p.adj"
# Make the significance letter
CLD.m = cldList(p.adj ~ Comparison, data=Tukey.M)
CLD.m
# re-order
CLD.m.ed = CLD.m[c(3,1,2),] # re-order
CLD.m.ed
#_______________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.M)
str(tab.distance.M)
tab.distance.M
# Post Hoc Dunn Test
set.seed(1333)
dunn.M=dunnTest(Distance ~ Label,
                data = tab.distance.M, method = "bh")
library(rcompanion)
CLD.m.ed = cldList(P.adj ~ Comparison, data=dunn.M$res)
CLD.m.ed
#_______________________________________________________________________________
# Plot
sumData.M <- ddply(tab.distance.M, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.m.ed) = sumData.M$Label
CLD.m.ed = rownames_to_column(CLD.m.ed, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

dist.CAP.bulk.m<- ggplot(tab.distance.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.m.ed,aes(x=Label, y = sumData.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.129", hjust = 0, size = 6, fontface='italic')+
  ggtitle("CONMIN (M)")
dist.CAP.bulk.m

# 3. AOA Bulk Soil - Bray Curtis Distance - BIODYN (D)
aoa.cap.bulk.dist
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.bulk.dist, item_groups)
#Control
tab.distance = as_tibble(d.calcul) 
tab.distance.D <- subset(tab.distance, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.D$Label <- factor(tab.distance.D$Label)
tab.distance.D$Label
# One-way ANOVA
set.seed(13)
D.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.D)
summary(D.dist.cap.aov) # not significant, p-val < 0.0001 f-val= 11.88
# Post-Hoc Test
D.dist.cap.tuk <- TukeyHSD(D.dist.cap.aov)
D.dist.cap.tuk
Tukey.D <- as.data.frame(D.dist.cap.tuk$Label)
Tukey.D <- rownames_to_column(Tukey.D, var = "Comparison")
colnames(Tukey.D)[5] <- "p.adj"
# Make the significance letter
CLD.d = cldList(p.adj ~ Comparison, data=Tukey.D)
CLD.d
# re-order
CLD.d.ed = CLD.d[c(3,1,2),] # re-order
CLD.d.ed
#_______________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.D)
str(tab.distance.D)
tab.distance.D
# Post Hoc Dunn Test
set.seed(1333)
dunn.D=dunnTest(Distance ~ Label,
                data = tab.distance.D, method = "bh")
library(rcompanion)
CLD.d = cldList(P.adj ~ Comparison, data=dunn.D$res)
CLD.d
#_______________________________________________________________________________
# Plot
sumData.D <- ddply(tab.distance.D, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.d.ed) = sumData.D$Label
CLD.d.ed = rownames_to_column(CLD.d.ed, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

dist.CAP.bulk.d<- ggplot(tab.distance.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.d.ed,aes(x=Label, y = sumData.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value < 0.001", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
dist.CAP.bulk.d
# Combine figures
library(patchwork)
dist.CAP.AOA <- dist.CAP.bulk.d | dist.CAP.bulk.k | dist.CAP.bulk.m
dist.CAP.AOA
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist.CAP.Bulk.tiff",
       dist.CAP.AOA, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

#_______________________________________________________________________________
#_______________________________________________________________________________


# 1. AOA Rhizosphere - CAP Distance - CONFYM (K)
aoa.cap.rh.dist
#Sample group
item_groups.rh <- sample_data(aoa.physeq_rh1)
item_groups.rh <- item_groups.rh$x
#calculate dist between groups
d.calcul.rh <- dist_groups(aoa.cap.rh.dist, item_groups.rh)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.rh = as_tibble(d.calcul.rh) 
tab.distance.rh.K <- subset(tab.distance.rh, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.rh.K$Label <- factor(tab.distance.rh.K$Label)
str(tab.distance.rh.K) 
# One-way ANOVA
set.seed(13)
rh.K.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.rh.K)
summary(rh.K.dist.cap.aov) # significant, p-val = 0.021, f-val= 3.89
# Post-Hoc Test
rh.K.dist.cap.tuk <- TukeyHSD(rh.K.dist.cap.aov)
Tukey.rh.K <- as.data.frame(rh.K.dist.cap.tuk$Label)
Tukey.rh.K <- rownames_to_column(Tukey.rh.K, var = "Comparison")
colnames(Tukey.rh.K)[5] <- "p.adj"
# Make the significance letter
CLD.rh.k = cldList(p.adj ~ Comparison, data=Tukey.rh.K)
CLD.rh.k
# re-order
CLD.rh.k.ed = CLD.rh.k[c(3,1,2),] # re-order
CLD.rh.k.ed
# Plot
sumData.rh.K <- ddply(tab.distance.rh.K, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
sumData.rh.K
rownames(CLD.rh.k.ed) = sumData.rh.K$Label
CLD.rh.k.ed = rownames_to_column(CLD.rh.k.ed, var="Label")
CLD.rh.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

dist.CAP.rhizo.k<- ggplot(tab.distance.rh.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=CLD.rh.k.ed,aes(x=Label, y = sumData.rh.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.021", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
dist.CAP.rhizo.k


# 2. AOA Rhizosphere - CAP Distance - CONMIN (M)
aoa.cap.rh.dist
#Sample group
item_groups.rh <- sample_data(aoa.physeq_rh1)
item_groups.rh <- item_groups.rh$x
#calculate dist between groups
d.calcul.rh <- dist_groups(aoa.cap.rh.dist, item_groups.rh)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.rh = as_tibble(d.calcul.rh) 
tab.distance.rh.M <- subset(tab.distance.rh, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.rh.M$Label <- factor(tab.distance.rh.M$Label)
str(tab.distance.rh.M) 
# One-way ANOVA
set.seed(13)
rh.M.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.rh.M)
summary(rh.M.dist.cap.aov) # not significant, p-val = 0.43, f-val= 0.84
# Post-Hoc Test
rh.M.dist.cap.tuk <- TukeyHSD(rh.M.dist.cap.aov)
Tukey.rh.M <- as.data.frame(rh.M.dist.cap.tuk$Label)
Tukey.rh.M <- rownames_to_column(Tukey.rh.M, var = "Comparison")
colnames(Tukey.rh.M)[5] <- "p.adj"
# Make the significance letter
CLD.rh.m = cldList(p.adj ~ Comparison, data=Tukey.rh.M)
CLD.rh.m
# re-order
CLD.rh.m.ed = CLD.rh.m[c(3,1,2),] # re-order
CLD.rh.m.ed
# Plot
sumData.rh.M <- ddply(tab.distance.rh.M, "Label", summarise,
                      Max = max(Distance),
                      N    = length(Distance),
                      Mean = mean(Distance),
                      Sd   = sd(Distance),
                      Se   = Sd / sqrt(N))
sumData.rh.M
rownames(CLD.rh.m.ed) = sumData.rh.M$Label
CLD.rh.m.ed = rownames_to_column(CLD.rh.m.ed, var="Label")
CLD.rh.m.ed
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

dist.CAP.rhizo.m<- ggplot(tab.distance.rh.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=CLD.rh.m.ed,aes(x=Label, y = sumData.rh.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.43", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONMIN (M)")
dist.CAP.rhizo.m

# 3. AOA Rhizosphere - CAP Distance - BIODYN (D)
aoa.cap.rh.dist
#Sample group
item_groups.rh <- sample_data(aoa.physeq_rh1)
item_groups.rh <- item_groups.rh$x
#calculate dist between groups
d.calcul.rh <- dist_groups(aoa.cap.rh.dist, item_groups.rh)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.rh = as_tibble(d.calcul.rh) 
tab.distance.rh.D <- subset(tab.distance.rh, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.rh.D$Label <- factor(tab.distance.rh.D$Label)
str(tab.distance.rh.D) 
# One-way ANOVA
set.seed(13)
rh.D.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.rh.D)
summary(rh.D.dist.cap.aov) # significant, p-val = 0.034, f-val= 3.4
# Post-Hoc Test
rh.D.dist.cap.tuk <- TukeyHSD(rh.D.dist.cap.aov)
Tukey.rh.D <- as.data.frame(rh.D.dist.cap.tuk$Label)
Tukey.rh.D <- rownames_to_column(Tukey.rh.D, var = "Comparison")
colnames(Tukey.rh.D)[5] <- "p.adj"
# Make the significance letter
CLD.rh.d = cldList(p.adj ~ Comparison, data=Tukey.rh.D)
CLD.rh.d
# re-order
CLD.rh.d.ed = CLD.rh.d[c(3,1,2),] # re-order
CLD.rh.d.ed
# Plot
sumData.rh.D <- ddply(tab.distance.rh.D, "Label", summarise,
                      Max = max(Distance),
                      N    = length(Distance),
                      Mean = mean(Distance),
                      Sd   = sd(Distance),
                      Se   = Sd / sqrt(N))
sumData.rh.D
rownames(CLD.rh.d.ed) = sumData.rh.D$Label
CLD.rh.d.ed = rownames_to_column(CLD.rh.d.ed, var="Label")
CLD.rh.d.ed
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

dist.CAP.rhizo.d<- ggplot(tab.distance.rh.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.rh.d.ed,aes(x=Label, y = sumData.rh.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.034", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
dist.CAP.rhizo.d
# Combine figures
library(patchwork)
dist.CAP.AOA.rhizo <- dist.CAP.rhizo.d | dist.CAP.rhizo.k | dist.CAP.rhizo.m
dist.CAP.AOA.rhizo
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist.CAP.Rhizo.tiff",
       dist.CAP.AOA.rhizo, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

################################################################################
################################################################################

# 1. AOB Bulk Soil - CAP Distance - CONFYM (K)
aob.cap.bulk.dist
#Sample group
aob.item_groups <- sample_data(aob.physeq_bulk1)
aob.item_groups <- aob.item_groups$x
#calculate dist between groups
d.calcul.aob <- dist_groups(aob.cap.bulk.dist, aob.item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.aob = as_tibble(d.calcul.aob) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.aob.K <- subset(tab.distance.aob, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.aob.K$Label <- factor(tab.distance.aob.K$Label)
str(tab.distance.aob.K) 
# One-way ANOVA
set.seed(13)
aob.K.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.aob.K)
summary(aob.K.dist.cap.aov) # significant, p-val = 0.01, f-val= 3.9
# Post-Hoc Test
aob.K.dist.cap.tuk <- TukeyHSD(aob.K.dist.cap.aov)
Tukey.aob.K <- as.data.frame(aob.K.dist.cap.tuk$Label)
Tukey.aob.K <- rownames_to_column(Tukey.aob.K, var = "Comparison")
colnames(Tukey.aob.K)[5] <- "p.adj"
# Make the significance letter
CLD.aob.k = cldList(p.adj ~ Comparison, data=Tukey.aob.K)
CLD.aob.k
# re-order
CLD.aob.k.ed = CLD.aob.k[c(3,1,2),] # re-order
CLD.aob.k.ed
# Plot
sumData.aob.K <- ddply(tab.distance.aob.K, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
sumData.aob.K
rownames(CLD.aob.k.ed) = sumData.aob.K$Label
CLD.aob.k.ed = rownames_to_column(CLD.aob.k.ed, var="Label")
CLD.aob.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

dist.aob.CAP.bulk.k<- ggplot(tab.distance.aob.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=CLD.aob.k.ed,aes(x=Label, y = sumData.aob.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.01", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
dist.aob.CAP.bulk.k


# 2. AOB Bulk Soil - CAP Distance - CONMIN (M)

aob.cap.bulk.dist
#Sample group
aob.item_groups <- sample_data(aob.physeq_bulk1)
aob.item_groups <- aob.item_groups$x
#calculate dist between groups
d.calcul.aob <- dist_groups(aob.cap.bulk.dist, aob.item_groups)
#Control
tab.distance.aob = as_tibble(d.calcul.aob) 
tab.distance.aob.M <- subset(tab.distance.aob, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.aob.M$Label <- factor(tab.distance.aob.M$Label)
str(tab.distance.aob.M) 
# One-way ANOVA
set.seed(13)
aob.M.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.aob.M)
summary(aob.M.dist.cap.aov) # not significant, p-val = 0.41, f-val= 0.89
# Post-Hoc Test
aob.M.dist.cap.tuk <- TukeyHSD(aob.M.dist.cap.aov)
Tukey.aob.M <- as.data.frame(aob.M.dist.cap.tuk$Label)
Tukey.aob.M <- rownames_to_column(Tukey.aob.M, var = "Comparison")
colnames(Tukey.aob.M)[5] <- "p.adj"
# Make the significance letter
CLD.aob.m = cldList(p.adj ~ Comparison, data=Tukey.aob.M)
CLD.aob.m
# re-order
CLD.aob.m.ed = CLD.aob.m[c(3,1,2),] # re-order
CLD.aob.m.ed
# Plot
sumData.aob.M <- ddply(tab.distance.aob.M, "Label", summarise,
                       Max = max(Distance),
                       N    = length(Distance),
                       Mean = mean(Distance),
                       Sd   = sd(Distance),
                       Se   = Sd / sqrt(N))
sumData.aob.M
rownames(CLD.aob.m.ed) = sumData.aob.M$Label
CLD.aob.m.ed = rownames_to_column(CLD.aob.m.ed, var="Label")
CLD.aob.m.ed
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

dist.aob.CAP.bulk.m<- ggplot(tab.distance.aob.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.aob.m.ed,aes(x=Label, y = sumData.aob.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.41", hjust = 0, size = 6, fontface='italic')+
  ggtitle("CONMIN (M)")
dist.aob.CAP.bulk.m

# 3. AOA Bulk Soil - Bray Curtis Distance - BIODYN (D)

aob.cap.bulk.dist
#Sample group
aob.item_groups <- sample_data(aob.physeq_bulk1)
aob.item_groups <- aob.item_groups$x
#calculate dist between groups
d.calcul.aob <- dist_groups(aob.cap.bulk.dist, aob.item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.aob = as_tibble(d.calcul.aob) 
tab.distance.aob.D <- subset(tab.distance.aob, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.aob.D$Label <- factor(tab.distance.aob.D$Label)
str(tab.distance.aob.D) 
# One-way ANOVA
set.seed(13)
aob.D.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.aob.D)
summary(aob.D.dist.cap.aov) # not significant, p-val = 0.1, f-val= 2.24
# Post-Hoc Test
aob.D.dist.cap.tuk <- TukeyHSD(aob.D.dist.cap.aov)
Tukey.aob.D <- as.data.frame(aob.D.dist.cap.tuk$Label)
Tukey.aob.D <- rownames_to_column(Tukey.aob.D, var = "Comparison")
colnames(Tukey.aob.D)[5] <- "p.adj"
# Make the significance letter
CLD.aob.d = cldList(p.adj ~ Comparison, data=Tukey.aob.D)
CLD.aob.d
# re-order
CLD.aob.d.ed = CLD.aob.d[c(3,1,2),] # re-order
CLD.aob.d.ed
# Plot
sumData.aob.D <- ddply(tab.distance.aob.D, "Label", summarise,
                       Max = max(Distance),
                       N    = length(Distance),
                       Mean = mean(Distance),
                       Sd   = sd(Distance),
                       Se   = Sd / sqrt(N))
sumData.aob.D
rownames(CLD.aob.d.ed) = sumData.aob.D$Label
CLD.aob.d.ed = rownames_to_column(CLD.aob.d.ed, var="Label")
CLD.aob.d.ed
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

dist.aob.CAP.bulk.d<- ggplot(tab.distance.aob.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.aob.d.ed,aes(x=Label, y = sumData.aob.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.1", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
dist.aob.CAP.bulk.d
# Combine figures
library(patchwork)
dist.CAP.AOB <- dist.aob.CAP.bulk.d | dist.aob.CAP.bulk.k | dist.aob.CAP.bulk.m
dist.CAP.AOB
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist.CAP.Bulk.tiff",
       dist.CAP.AOB, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

#_______________________________________________________________________________
#_______________________________________________________________________________


# 1. AOB Rhizosphere - CAP Distance - CONFYM (K)
aob.cap.rh.dist
#Sample group
aob.item_groups.rh <- sample_data(aob.physeq_rh1)
aob.item_groups.rh <- aob.item_groups.rh$x
#calculate dist between groups
aob.d.calcul.rh <- dist_groups(aob.cap.rh.dist, aob.item_groups.rh)
#Control
aob.tab.distance.rh = as_tibble(aob.d.calcul.rh) 
aob.tab.distance.rh.K <- subset(aob.tab.distance.rh, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
aob.tab.distance.rh.K$Label <- factor(aob.tab.distance.rh.K$Label)
str(aob.tab.distance.rh.K) 
# One-way ANOVA
set.seed(13)
aob.rh.K.dist.cap.aov <- aov(Distance ~ Label, data = aob.tab.distance.rh.K)
summary(aob.rh.K.dist.cap.aov) # not significant, p-val = 0.057, f-val= 2.89
# Post-Hoc Test
aob.rh.K.dist.cap.tuk <- TukeyHSD(aob.rh.K.dist.cap.aov)
aob.Tukey.rh.K <- as.data.frame(aob.rh.K.dist.cap.tuk$Label)
aob.Tukey.rh.K <- rownames_to_column(aob.Tukey.rh.K, var = "Comparison")
colnames(aob.Tukey.rh.K)[5] <- "p.adj"
# Make the significance letter
aob.CLD.rh.k = cldList(p.adj ~ Comparison, data=aob.Tukey.rh.K)
aob.CLD.rh.k
# re-order
aob.CLD.rh.k.ed = aob.CLD.rh.k[c(3,1,2),] # re-order
aob.CLD.rh.k.ed
# Plot
aob.sumData.rh.K <- ddply(aob.tab.distance.rh.K, "Label", summarise,
                      Max = max(Distance),
                      N    = length(Distance),
                      Mean = mean(Distance),
                      Sd   = sd(Distance),
                      Se   = Sd / sqrt(N))
aob.sumData.rh.K
rownames(aob.CLD.rh.k.ed) = aob.sumData.rh.K$Label
aob.CLD.rh.k.ed = rownames_to_column(aob.CLD.rh.k.ed, var="Label")
aob.CLD.rh.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

aob.dist.CAP.rhizo.k<- ggplot(aob.tab.distance.rh.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=aob.CLD.rh.k.ed,aes(x=Label, y = aob.sumData.rh.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.057", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
aob.dist.CAP.rhizo.k


# 2. AOB Rhizosphere - CAP Distance - CONMIN (M)
aob.cap.rh.dist
#Sample group
aob.item_groups.rh <- sample_data(aob.physeq_rh1)
aob.item_groups.rh <- aob.item_groups.rh$x
#calculate dist between groups
aob.d.calcul.rh <- dist_groups(aob.cap.rh.dist, aob.item_groups.rh)
#Control
aob.tab.distance.rh = as_tibble(aob.d.calcul.rh) 
aob.tab.distance.rh.M <- subset(aob.tab.distance.rh, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
aob.tab.distance.rh.M$Label <- factor(aob.tab.distance.rh.M$Label)
str(aob.tab.distance.rh.M) 
# One-way ANOVA
set.seed(13)
aob.rh.M.dist.cap.aov <- aov(Distance ~ Label, data = aob.tab.distance.rh.M)
summary(aob.rh.M.dist.cap.aov) # not significant, p-val = 0.31, f-val= 1.175
# Post-Hoc Test
aob.rh.M.dist.cap.tuk <- TukeyHSD(aob.rh.M.dist.cap.aov)
aob.Tukey.rh.M <- as.data.frame(aob.rh.M.dist.cap.tuk$Label)
aob.Tukey.rh.M <- rownames_to_column(aob.Tukey.rh.M, var = "Comparison")
colnames(aob.Tukey.rh.M)[5] <- "p.adj"
# Make the significance letter
aob.CLD.rh.m = cldList(p.adj ~ Comparison, data=aob.Tukey.rh.M)
aob.CLD.rh.m
# re-order
aob.CLD.rh.m.ed = aob.CLD.rh.m[c(3,1,2),] # re-order
aob.CLD.rh.m.ed
# Plot
aob.sumData.rh.M <- ddply(aob.tab.distance.rh.M, "Label", summarise,
                          Max = max(Distance),
                          N    = length(Distance),
                          Mean = mean(Distance),
                          Sd   = sd(Distance),
                          Se   = Sd / sqrt(N))
aob.sumData.rh.M
rownames(aob.CLD.rh.m.ed) = aob.sumData.rh.M$Label
aob.CLD.rh.m.ed = rownames_to_column(aob.CLD.rh.m.ed, var="Label")
aob.CLD.rh.m.ed
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

aob.dist.CAP.rhizo.m<- ggplot(aob.tab.distance.rh.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.10)+
  geom_text(data=aob.CLD.rh.m.ed,aes(x=Label, y = aob.sumData.rh.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.31", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONMIN (M)")
aob.dist.CAP.rhizo.m

# 3. AOB Rhizosphere - CAP Distance - BIODYN (D)

aob.cap.rh.dist
#Sample group
aob.item_groups.rh <- sample_data(aob.physeq_rh1)
aob.item_groups.rh <- aob.item_groups.rh$x
#calculate dist between groups
aob.d.calcul.rh <- dist_groups(aob.cap.rh.dist, aob.item_groups.rh)
#Control
aob.tab.distance.rh = as_tibble(aob.d.calcul.rh) 
aob.tab.distance.rh.D <- subset(aob.tab.distance.rh, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
aob.tab.distance.rh.D$Label <- factor(aob.tab.distance.rh.D$Label)
str(aob.tab.distance.rh.D) 
# One-way ANOVA
set.seed(13)
aob.rh.D.dist.cap.aov <- aov(Distance ~ Label, data = aob.tab.distance.rh.D)
summary(aob.rh.D.dist.cap.aov) # significant, p-val = 0.0002, f-val= 8.57
# Post-Hoc Test
aob.rh.D.dist.cap.tuk <- TukeyHSD(aob.rh.D.dist.cap.aov)
aob.Tukey.rh.D <- as.data.frame(aob.rh.D.dist.cap.tuk$Label)
aob.Tukey.rh.D <- rownames_to_column(aob.Tukey.rh.D, var = "Comparison")
colnames(aob.Tukey.rh.D)[5] <- "p.adj"
# Make the significance letter
aob.CLD.rh.d = cldList(p.adj ~ Comparison, data=aob.Tukey.rh.D)
aob.CLD.rh.d
# re-order
aob.CLD.rh.d.ed = aob.CLD.rh.d[c(3,1,2),] # re-order
aob.CLD.rh.d.ed
# Plot
aob.sumData.rh.D <- ddply(aob.tab.distance.rh.D, "Label", summarise,
                          Max = max(Distance),
                          N    = length(Distance),
                          Mean = mean(Distance),
                          Sd   = sd(Distance),
                          Se   = Sd / sqrt(N))
aob.sumData.rh.D
rownames(aob.CLD.rh.d.ed) = aob.sumData.rh.D$Label
aob.CLD.rh.d.ed = rownames_to_column(aob.CLD.rh.d.ed, var="Label")
aob.CLD.rh.d.ed
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

aob.dist.CAP.rhizo.d<- ggplot(aob.tab.distance.rh.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=aob.CLD.rh.d.ed,aes(x=Label, y = aob.sumData.rh.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.10)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.10,label= "P-value = 0.0002", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
aob.dist.CAP.rhizo.d

# Combine figures
library(patchwork)
dist.CAP.AOB.rhizo <- aob.dist.CAP.rhizo.d | aob.dist.CAP.rhizo.k | aob.dist.CAP.rhizo.m
dist.CAP.AOB.rhizo
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist.CAP.Rhizo.tiff",
       dist.CAP.AOB.rhizo, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

################################################################################
################################################################################

# 1. COMAMMOX Bulk Soil - CAP Distance - CONFYM (K)

com.cap.bulk.dist
#Sample group
com.item_groups <- sample_data(com.physeq_bulk1)
com.item_groups <- com.item_groups$x
#calculate dist between groups
d.calcul.com <- dist_groups(com.cap.bulk.dist, com.item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.com = as_tibble(d.calcul.com) 
tab.distance.com.K <- subset(tab.distance.com, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.com.K$Label <- factor(tab.distance.com.K$Label)
str(tab.distance.com.K) 
# One-way ANOVA
set.seed(13)
com.K.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.com.K)
summary(com.K.dist.cap.aov) # not significant, p-val = 0.79, f-val= 0.22
# Post-Hoc Test
com.K.dist.cap.tuk <- TukeyHSD(com.K.dist.cap.aov)
Tukey.com.K <- as.data.frame(com.K.dist.cap.tuk$Label)
Tukey.com.K <- rownames_to_column(Tukey.com.K, var = "Comparison")
colnames(Tukey.com.K)[5] <- "p.adj"
# Make the significance letter
CLD.com.k = cldList(p.adj ~ Comparison, data=Tukey.com.K)
CLD.com.k
# re-order
CLD.com.k.ed = CLD.com.k[c(3,1,2),] # re-order
CLD.com.k.ed
# Plot
sumData.com.K <- ddply(tab.distance.com.K, "Label", summarise,
                       Max = max(Distance),
                       N    = length(Distance),
                       Mean = mean(Distance),
                       Sd   = sd(Distance),
                       Se   = Sd / sqrt(N))
sumData.com.K
rownames(CLD.com.k.ed) = sumData.com.K$Label
CLD.com.k.ed = rownames_to_column(CLD.com.k.ed, var="Label")
CLD.com.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

dist.com.CAP.bulk.k<- ggplot(tab.distance.com.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.4)+
  geom_text(data=CLD.com.k.ed,aes(x=Label, y = sumData.com.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value = 0.79", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
dist.com.CAP.bulk.k


# 2. COMAMMOX Bulk Soil - CAP Distance - CONMIN (M)

com.cap.bulk.dist
#Sample group
com.item_groups <- sample_data(com.physeq_bulk1)
com.item_groups <- com.item_groups$x
#calculate dist between groups
d.calcul.com <- dist_groups(com.cap.bulk.dist, com.item_groups)
#Control
tab.distance.com = as_tibble(d.calcul.com) 
tab.distance.com.M <- subset(tab.distance.com, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.com.M$Label <- factor(tab.distance.com.M$Label)
str(tab.distance.com.M) 
# One-way ANOVA
set.seed(13)
com.M.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.com.M)
summary(com.M.dist.cap.aov) # not significant, p-val = 0.39, f-val= 0.93
# Post-Hoc Test
com.M.dist.cap.tuk <- TukeyHSD(com.M.dist.cap.aov)
Tukey.com.M <- as.data.frame(com.M.dist.cap.tuk$Label)
Tukey.com.M <- rownames_to_column(Tukey.com.M, var = "Comparison")
colnames(Tukey.com.M)[5] <- "p.adj"
# Make the significance letter
CLD.com.m = cldList(p.adj ~ Comparison, data=Tukey.com.M)
CLD.com.m
# re-order
CLD.com.m.ed = CLD.com.m[c(3,1,2),] # re-order
CLD.com.m.ed
# Plot
sumData.com.M <- ddply(tab.distance.com.M, "Label", summarise,
                       Max = max(Distance),
                       N    = length(Distance),
                       Mean = mean(Distance),
                       Sd   = sd(Distance),
                       Se   = Sd / sqrt(N))
sumData.com.M
rownames(CLD.com.m.ed) = sumData.com.M$Label
CLD.com.m.ed = rownames_to_column(CLD.com.m.ed, var="Label")
CLD.com.m.ed
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

dist.com.CAP.bulk.m<- ggplot(tab.distance.com.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.com.m.ed,aes(x=Label, y = sumData.com.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  ylim(0,1.4)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value = 0.39", hjust = 0, size = 6, fontface='italic')+
  ggtitle("CONMIN (M)")
dist.com.CAP.bulk.m

# 3. COMAMMOX Bulk Soil - Bray Curtis Distance - BIODYN (D)

com.cap.bulk.dist
#Sample group
com.item_groups <- sample_data(com.physeq_bulk1)
com.item_groups <- com.item_groups$x
#calculate dist between groups
d.calcul.com <- dist_groups(com.cap.bulk.dist, com.item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance.com = as_tibble(d.calcul.com) 
tab.distance.com.D <- subset(tab.distance.com, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.com.D$Label <- factor(tab.distance.com.D$Label)
str(tab.distance.com.D) 
# One-way ANOVA
set.seed(13)
com.D.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.com.D)
summary(com.D.dist.cap.aov) # significant, p-val < 0.0001, f-val= 12.56
# Post-Hoc Test
com.D.dist.cap.tuk <- TukeyHSD(com.D.dist.cap.aov)
Tukey.com.D <- as.data.frame(com.D.dist.cap.tuk$Label)
Tukey.com.D <- rownames_to_column(Tukey.com.D, var = "Comparison")
colnames(Tukey.com.D)[5] <- "p.adj"
# Make the significance letter
CLD.com.d = cldList(p.adj ~ Comparison, data=Tukey.com.D)
CLD.com.d
# re-order
CLD.com.d.ed = CLD.com.d[c(3,1,2),] # re-order
CLD.com.d.ed
# Plot
sumData.com.D <- ddply(tab.distance.com.D, "Label", summarise,
                       Max = max(Distance),
                       N    = length(Distance),
                       Mean = mean(Distance),
                       Sd   = sd(Distance),
                       Se   = Sd / sqrt(N))
sumData.com.D
rownames(CLD.com.d.ed) = sumData.com.D$Label
CLD.com.d.ed = rownames_to_column(CLD.com.d.ed, var="Label")
CLD.com.d.ed
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

dist.com.CAP.bulk.d<- ggplot(tab.distance.com.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.com.d.ed,aes(x=Label, y = sumData.com.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.4)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value < 0.0001", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
dist.com.CAP.bulk.d

# Combine figures
library(patchwork)
dist.CAP.COM <- dist.com.CAP.bulk.d | dist.com.CAP.bulk.k | dist.com.CAP.bulk.m
dist.CAP.COM
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist.CAP.Bulk.tiff",
       dist.CAP.COM, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

#_______________________________________________________________________________
#_______________________________________________________________________________

# 1. COMAMMOX Rhizosphere - CAP Distance - CONFYM (K)
com.cap.rh.dist
#Sample group
com.item_groups.rh <- sample_data(com.physeq_rh1)
com.item_groups.rh <- com.item_groups.rh$x
#calculate dist between groups
com.d.calcul.rh <- dist_groups(com.cap.rh.dist, com.item_groups.rh)
#Control
com.tab.distance.rh = as_tibble(com.d.calcul.rh) 
com.tab.distance.rh.K <- subset(com.tab.distance.rh, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
com.tab.distance.rh.K$Label <- factor(com.tab.distance.rh.K$Label)
str(com.tab.distance.rh.K) 
# One-way ANOVA
set.seed(13)
com.rh.K.dist.cap.aov <- aov(Distance ~ Label, data = com.tab.distance.rh.K)
summary(com.rh.K.dist.cap.aov) # not significant, p-val = 0.01, f-val= 4.635
# Post-Hoc Test
com.rh.K.dist.cap.tuk <- TukeyHSD(com.rh.K.dist.cap.aov)
com.Tukey.rh.K <- as.data.frame(com.rh.K.dist.cap.tuk$Label)
com.Tukey.rh.K <- rownames_to_column(com.Tukey.rh.K, var = "Comparison")
colnames(com.Tukey.rh.K)[5] <- "p.adj"
# Make the significance letter
com.CLD.rh.k = cldList(p.adj ~ Comparison, data=com.Tukey.rh.K)
com.CLD.rh.k
# re-order
com.CLD.rh.k.ed = com.CLD.rh.k[c(3,1,2),] # re-order
com.CLD.rh.k.ed
# Plot
com.sumData.rh.K <- ddply(com.tab.distance.rh.K, "Label", summarise,
                          Max = max(Distance),
                          N    = length(Distance),
                          Mean = mean(Distance),
                          Sd   = sd(Distance),
                          Se   = Sd / sqrt(N))
com.sumData.rh.K
rownames(com.CLD.rh.k.ed) = com.sumData.rh.K$Label
com.CLD.rh.k.ed = rownames_to_column(com.CLD.rh.k.ed, var="Label")
com.CLD.rh.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

com.dist.CAP.rhizo.k<- ggplot(com.tab.distance.rh.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.4)+
  geom_text(data=com.CLD.rh.k.ed,aes(x=Label, y = com.sumData.rh.K$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value = 0.01", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
com.dist.CAP.rhizo.k


# 2. COMAMMOX Rhizosphere - CAP Distance - CONMIN (M)
com.cap.rh.dist
#Sample group
com.item_groups.rh <- sample_data(com.physeq_rh1)
com.item_groups.rh <- com.item_groups.rh$x
#calculate dist between groups
com.d.calcul.rh <- dist_groups(com.cap.rh.dist, com.item_groups.rh)
#Control
com.tab.distance.rh = as_tibble(com.d.calcul.rh) 
com.tab.distance.rh.M <- subset(com.tab.distance.rh, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
com.tab.distance.rh.M$Label <- factor(com.tab.distance.rh.M$Label)
str(com.tab.distance.rh.M) 
# One-way ANOVA
set.seed(13)
com.rh.M.dist.cap.aov <- aov(Distance ~ Label, data = com.tab.distance.rh.M)
summary(com.rh.M.dist.cap.aov) # not significant, p-val = 0.25, f-val= 1.37
# Post-Hoc Test
com.rh.M.dist.cap.tuk <- TukeyHSD(com.rh.M.dist.cap.aov)
com.Tukey.rh.M <- as.data.frame(com.rh.M.dist.cap.tuk$Label)
com.Tukey.rh.M <- rownames_to_column(com.Tukey.rh.M, var = "Comparison")
colnames(com.Tukey.rh.M)[5] <- "p.adj"
# Make the significance letter
com.CLD.rh.m = cldList(p.adj ~ Comparison, data=com.Tukey.rh.M)
com.CLD.rh.m
# re-order
com.CLD.rh.m.ed = com.CLD.rh.m[c(3,1,2),] # re-order
com.CLD.rh.m.ed
# Plot
com.sumData.rh.M <- ddply(com.tab.distance.rh.M, "Label", summarise,
                          Max = max(Distance),
                          N    = length(Distance),
                          Mean = mean(Distance),
                          Sd   = sd(Distance),
                          Se   = Sd / sqrt(N))
com.sumData.rh.M
rownames(com.CLD.rh.m.ed) = com.sumData.rh.M$Label
com.CLD.rh.m.ed = rownames_to_column(com.CLD.rh.m.ed, var="Label")
com.CLD.rh.m.ed
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

com.dist.CAP.rhizo.m<- ggplot(com.tab.distance.rh.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,1.4)+
  geom_text(data=com.CLD.rh.m.ed,aes(x=Label, y = com.sumData.rh.M$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value = 0.25", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONMIN (M)")
com.dist.CAP.rhizo.m

# 3. COMAMMOX Rhizosphere - CAP Distance - BIODYN (D)

com.cap.rh.dist
#Sample group
com.item_groups.rh <- sample_data(com.physeq_rh1)
com.item_groups.rh <- com.item_groups.rh$x
#calculate dist between groups
com.d.calcul.rh <- dist_groups(com.cap.rh.dist, com.item_groups.rh)
#Control
com.tab.distance.rh = as_tibble(com.d.calcul.rh) 
com.tab.distance.rh.D <- subset(com.tab.distance.rh, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
com.tab.distance.rh.D$Label <- factor(com.tab.distance.rh.D$Label)
str(com.tab.distance.rh.D) 
# One-way ANOVA
set.seed(13)
com.rh.D.dist.cap.aov <- aov(Distance ~ Label, data = com.tab.distance.rh.D)
summary(com.rh.D.dist.cap.aov) # not significant, p-val = 0.11, f-val= 2.14
# Post-Hoc Test
com.rh.D.dist.cap.tuk <- TukeyHSD(com.rh.D.dist.cap.aov)
com.Tukey.rh.D <- as.data.frame(com.rh.D.dist.cap.tuk$Label)
com.Tukey.rh.D <- rownames_to_column(com.Tukey.rh.D, var = "Comparison")
colnames(com.Tukey.rh.D)[5] <- "p.adj"
# Make the significance letter
com.CLD.rh.d = cldList(p.adj ~ Comparison, data=com.Tukey.rh.D)
com.CLD.rh.d
# re-order
com.CLD.rh.d.ed = com.CLD.rh.d[c(3,1,2),] # re-order
com.CLD.rh.d.ed
# Plot
com.sumData.rh.D <- ddply(com.tab.distance.rh.D, "Label", summarise,
                          Max = max(Distance),
                          N    = length(Distance),
                          Mean = mean(Distance),
                          Sd   = sd(Distance),
                          Se   = Sd / sqrt(N))
com.sumData.rh.D
rownames(com.CLD.rh.d.ed) = com.sumData.rh.D$Label
com.CLD.rh.d.ed = rownames_to_column(com.CLD.rh.d.ed, var="Label")
com.CLD.rh.d.ed
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

com.dist.CAP.rhizo.d<- ggplot(com.tab.distance.rh.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=com.CLD.rh.d.ed,aes(x=Label, y = com.sumData.rh.D$Max + 0.04, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,1.4)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=1.4,label= "P-value = 0.11", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
com.dist.CAP.rhizo.d

# Combine figures
library(patchwork)
dist.CAP.COM.rhizo <- com.dist.CAP.rhizo.d | com.dist.CAP.rhizo.k | com.dist.CAP.rhizo.m
dist.CAP.COM.rhizo
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist.CAP.Rhizo.tiff",
       dist.CAP.COM.rhizo, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)





















################################################################################
# AOA distance - everything in one plot

# 1. AOA Bulk Soil - Bray Curtis Distance 
aoa.bulk_dist_bc
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.bulk_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.C$Group2 <-  as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.C.ed)
tab.distance.C.ed2 <- tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.C.ed2)
str(tab.distance.C.ed2)
tab.distance.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C=dunnTest(Distance ~ Label,
                data = tab.distance.C.ed2, method = "bh")
dunn.C$res
str(t)
library(rcompanion)
CLD.C = cldList(P.adj ~ Comparison, data=dunn.C$res)
CLD.C.ed = CLD.C[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.ed
# Plot
sumData.C <- tab.distance.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C
rownames(CLD.C.ed) = sumData.C$Label
CLD.C.ed = rownames_to_column(CLD.C.ed, var="Label")
CLD.C.ed

#treatOrderC <- factor(c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                        #"Between cont.K and rain.K","Within cont.K","Within rain.K",
                        #"Between cont.M and rain.M","Within cont.M","Within rain.M"))
#CLD.C <- CLD.C[match(treatOrderC , CLD.C$Label),]
#CLD.C$Label <- reorder.factor(CLD.C$Label, new.order = treatOrderC)
#str(CLD.C)
CLD.C.ed2 <- CLD.C.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

#sumData.C.ed <- sumData.C %>% 
  #mutate(Label = factor(Label, 
                        #levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   #"Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   #"Between cont.M and rain.M","Within cont.M","Within rain.M")))
tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
dist.bulk.c<- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
                    #labels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                # "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                 #"Between cont.M and rain.M","Within cont.M","Within rain.M")) +
 facet_wrap(~ Treatment,scales = "free_x")+
 scale_x_discrete(labels = c("Between","Control","Drought",
                             "Between","Control","Drought",
                             "Between","Control","Drought"))+
  geom_text(data=CLD.C.ed2,aes(x=Label, y = sumData.C$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.bulk.c
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_bulk_bray.tiff",
       dist.bulk.c, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)

# 2. AOA Rhizosphere - Bray Curtis Distance 
aoa.rh_dist_bc
#Sample group
aoa.item_groups.rh <- sample_data(aoa.physeq_rh1)
aoa.item_groups.rh <- aoa.item_groups.rh$x
#calculate dist between groups
d.calcul.aoa.rh <- dist_groups(aoa.rh_dist_bc, aoa.item_groups.rh)
#Control
tab.distance.aoa.rh = as_tibble(d.calcul.aoa.rh) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.aoa.rh.C <- subset(tab.distance.aoa.rh, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.aoa.rh.C$Label <- factor(tab.distance.aoa.rh.C$Label)
tab.distance.aoa.rh.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.aoa.rh.C$Group2 <-  as.character(tab.distance.aoa.rh.C$Group2) 
tab.distance.aoa.rh.C.ed <- tab.distance.aoa.rh.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.aoa.rh.C.ed)
tab.distance.aoa.rh.C.ed2 <- tab.distance.aoa.rh.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.aoa.rh.C.ed2)
str(tab.distance.aoa.rh.C.ed2)
tab.distance.aoa.rh.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C.aoa.rh=dunnTest(Distance ~ Label,
                data = tab.distance.aoa.rh.C.ed2, method = "bh")
dunn.C.aoa.rh$res
library(rcompanion)
CLD.C.aoa.rh = cldList(P.adj ~ Comparison, data=dunn.C.aoa.rh$res)
CLD.C.aoa.rh
CLD.C.aoa.rh.ed = CLD.C.aoa.rh[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.aoa.rh.ed
# Plot
sumData.C.aoa.rh <- tab.distance.aoa.rh.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C.aoa.rh
rownames(CLD.C.aoa.rh.ed) = sumData.C.aoa.rh$Label
CLD.C.aoa.rh.ed = rownames_to_column(CLD.C.aoa.rh.ed, var="Label")
CLD.C.aoa.rh.ed

CLD.C.aoa.rh.ed2 <- CLD.C.aoa.rh.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.aoa.rh.C.ed2$Treatment <- as.factor(tab.distance.aoa.rh.C.ed2$Treatment)

# Plot
dist.rh.aoa<- ggplot(tab.distance.aoa.rh.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.aoa.rh.ed2,aes(x=Label, y = sumData.C.aoa.rh$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.rh.aoa
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_rhizo_bray.tiff",
       dist.rh.aoa, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)

################################################################################
# AOB distance - everything in one plot

# 1. AOB Bulk Soil - Bray Curtis Distance 
aob.bulk_dist_bc
#Sample group
item_groups <- sample_data(aob.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aob.bulk_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.C$Group2 <-  as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.C.ed)
tab.distance.C.ed2 <- tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(tab.distance.C.ed2)
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.C.ed2)
str(tab.distance.C.ed2)
tab.distance.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C=dunnTest(Distance ~ Label,
                data = tab.distance.C.ed2, method = "bh")
dunn.C$res
library(rcompanion)
CLD.C = cldList(P.adj ~ Comparison, data=dunn.C$res)
CLD.C.ed = CLD.C[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.ed
# Plot
sumData.C <- tab.distance.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C
rownames(CLD.C.ed) = sumData.C$Label
CLD.C.ed = rownames_to_column(CLD.C.ed, var="Label")
CLD.C.ed

CLD.C.ed2 <- CLD.C.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
aob.dist.bulk.c<- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.ed2,aes(x=Label, y = sumData.C$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
aob.dist.bulk.c
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist_bulk_bray.tiff",
       aob.dist.bulk.c, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)


# 2. AOB Rhizosphere - Bray Curtis Distance 
aob.rh_dist_bc
#Sample group
aob.item_groups.rh <- sample_data(aob.physeq_rh1)
aob.item_groups.rh <- aob.item_groups.rh$x
#calculate dist between groups
d.calcul.aob.rh <- dist_groups(aob.rh_dist_bc, aob.item_groups.rh)
#Control
tab.distance.aob.rh = as_tibble(d.calcul.aob.rh) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.aob.rh.C <- subset(tab.distance.aob.rh, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.aob.rh.C$Label <- factor(tab.distance.aob.rh.C$Label)
tab.distance.aob.rh.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.aob.rh.C$Group2 <-  as.character(tab.distance.aob.rh.C$Group2) 
tab.distance.aob.rh.C.ed <- tab.distance.aob.rh.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.aob.rh.C.ed)
tab.distance.aob.rh.C.ed2 <- tab.distance.aob.rh.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.aob.rh.C.ed2)
str(tab.distance.aob.rh.C.ed2)
tab.distance.aob.rh.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C.aob.rh=dunnTest(Distance ~ Label,
                       data = tab.distance.aob.rh.C.ed2, method = "bh")
dunn.C.aob.rh$res
library(rcompanion)
CLD.C.aob.rh = cldList(P.adj ~ Comparison, data=dunn.C.aob.rh$res)
CLD.C.aob.rh
CLD.C.aob.rh.ed = CLD.C.aob.rh[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.aob.rh.ed
# Plot
sumData.C.aob.rh <- tab.distance.aob.rh.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C.aob.rh
rownames(CLD.C.aob.rh.ed) = sumData.C.aob.rh$Label
CLD.C.aob.rh.ed = rownames_to_column(CLD.C.aob.rh.ed, var="Label")
CLD.C.aob.rh.ed

CLD.C.aob.rh.ed2 <- CLD.C.aob.rh.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))
CLD.C.aob.rh.ed2
tab.distance.aob.rh.C.ed2$Treatment <- as.factor(tab.distance.aob.rh.C.ed2$Treatment)

# Plot
dist.rh.aob<- ggplot(tab.distance.aob.rh.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.aob.rh.ed2,aes(x=Label, y = sumData.C.aob.rh$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.rh.aob
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist_rhizo_bray.tiff",
       dist.rh.aob, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)

################################################################################

# COMAMMOX distance - everything in one plot

# 1. COM Bulk Soil - Bray Curtis Distance 
com.bulk_dist_bc
#Sample group
item_groups <- sample_data(com.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(com.bulk_dist_bc, item_groups)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.C$Group2 <-  as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.C.ed)
tab.distance.C.ed2 <- tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(tab.distance.C.ed2)
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.C.ed2)
str(tab.distance.C.ed2)
tab.distance.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C=dunnTest(Distance ~ Label,
                data = tab.distance.C.ed2, method = "bh")
dunn.C$res
library(rcompanion)
CLD.C = cldList(P.adj ~ Comparison, data=dunn.C$res)
CLD.C.ed = CLD.C[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.ed
# Plot
sumData.C <- tab.distance.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C
rownames(CLD.C.ed) = sumData.C$Label
CLD.C.ed = rownames_to_column(CLD.C.ed, var="Label")
CLD.C.ed

CLD.C.ed2 <- CLD.C.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
com.dist.bulk.c<- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.ed2,aes(x=Label, y = sumData.C$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
com.dist.bulk.c
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist_bulk_bray.tiff",
       com.dist.bulk.c, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)


# 2. COM Rhizosphere - Bray Curtis Distance 
com.rh_dist_bc
#Sample group
com.item_groups.rh <- sample_data(com.physeq_rh1)
com.item_groups.rh <- com.item_groups.rh$x
#calculate dist between groups
d.calcul.com.rh <- dist_groups(com.rh_dist_bc, com.item_groups.rh)
#Control
tab.distance.com.rh = as_tibble(d.calcul.com.rh) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.com.rh.C <- subset(tab.distance.com.rh, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.com.rh.C$Label <- factor(tab.distance.com.rh.C$Label)
tab.distance.com.rh.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.com.rh.C$Group2 <-  as.character(tab.distance.com.rh.C$Group2) 
tab.distance.com.rh.C.ed <- tab.distance.com.rh.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.com.rh.C.ed)
tab.distance.com.rh.C.ed2 <- tab.distance.com.rh.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.com.rh.C.ed2)
str(tab.distance.com.rh.C.ed2)
tab.distance.com.rh.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C.com.rh=dunnTest(Distance ~ Label,
                       data = tab.distance.com.rh.C.ed2, method = "bh")
dunn.C.com.rh$res
library(rcompanion)
CLD.C.com.rh = cldList(P.adj ~ Comparison, data=dunn.C.com.rh$res)
CLD.C.com.rh
CLD.C.com.rh.ed = CLD.C.com.rh[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.com.rh.ed
# Plot
sumData.C.com.rh <- tab.distance.com.rh.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C.com.rh
rownames(CLD.C.com.rh.ed) = sumData.C.com.rh$Label
CLD.C.com.rh.ed = rownames_to_column(CLD.C.com.rh.ed, var="Label")
CLD.C.com.rh.ed

CLD.C.com.rh.ed2 <- CLD.C.com.rh.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))
CLD.C.com.rh.ed2
tab.distance.com.rh.C.ed2$Treatment <- as.factor(tab.distance.com.rh.C.ed2$Treatment)

# Plot
dist.rh.com<- ggplot(tab.distance.com.rh.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.com.rh.ed2,aes(x=Label, y = sumData.C.com.rh$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.rh.com
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist_rhizo_bray.tiff",
       dist.rh.com, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)





































































##################################################################################
# 1. AOA Bulk Soil - CAP Distance 
aoa.cap.bulk.dist
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.bulk.dist, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.C$Group2 <-  as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.C.ed)
tab.distance.C.ed2 <- tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.C.ed2)
str(tab.distance.C.ed2)
tab.distance.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C=dunnTest(Distance ~ Label,
                data = tab.distance.C.ed2, method = "bh")
dunn.C$res
str(t)
library(rcompanion)
CLD.C = cldList(P.adj ~ Comparison, data=dunn.C$res)
CLD.C.ed = CLD.C[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.ed
# Plot
sumData.C <- tab.distance.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C
rownames(CLD.C.ed) = sumData.C$Label
CLD.C.ed = rownames_to_column(CLD.C.ed, var="Label")
CLD.C.ed

CLD.C.ed2 <- CLD.C.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
dist.bulk.CAP<- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.ed2,aes(x=Label, y = sumData.C$Max + 0.04, label=Letter), vjust=0) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.bulk.CAP
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_bulk_CAP.tiff",
       dist.bulk.CAP, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)

# 2. AOA Rhizosphere - CAP Distance - CONFYM
aoa.cap.rh.dist
#Sample group
aoa.item_groups.rh <- sample_data(aoa.physeq_rh1)
aoa.item_groups.rh <- aoa.item_groups.rh$x
#calculate dist between groups
d.calcul.aoa.rh <- dist_groups(aoa.cap.rh.dist, aoa.item_groups.rh)
#Control
tab.distance.aoa.rh = as_tibble(d.calcul.aoa.rh) 
#tab.distance.C = tab.distance
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
tab.distance.aoa.rh.C <- subset(tab.distance.aoa.rh, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.aoa.rh.C$Label <- factor(tab.distance.aoa.rh.C$Label)
tab.distance.aoa.rh.C$Label
#tab.distance.A <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.aoa.rh.C$Group2 <-  as.character(tab.distance.aoa.rh.C$Group2) 
tab.distance.aoa.rh.C.ed <- tab.distance.aoa.rh.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.aoa.rh.C.ed)
tab.distance.aoa.rh.C.ed2 <- tab.distance.aoa.rh.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.aoa.rh.C.ed2)
str(tab.distance.aoa.rh.C.ed2)
tab.distance.aoa.rh.C.ed2
# Post Hoc Dunn Test
set.seed(1333)
dunn.C.aoa.rh=dunnTest(Distance ~ Label,
                       data = tab.distance.aoa.rh.C.ed2, method = "bh")
dunn.C.aoa.rh$res
library(rcompanion)
CLD.C.aoa.rh = cldList(P.adj ~ Comparison, data=dunn.C.aoa.rh$res)
CLD.C.aoa.rh
CLD.C.aoa.rh.ed = CLD.C.aoa.rh[c(1,4,7,2,5,8,3,6,9),] # re-order
CLD.C.aoa.rh.ed
# Plot
sumData.C.aoa.rh <- tab.distance.aoa.rh.C.ed2 %>%
  group_by(Label,Treatment) %>%
  summarise(Max = max(Distance),
            N = length(Distance),
            Mean = mean(Distance),
            Sd = sd(Distance),
            Se = Sd/sqrt(N))
sumData.C.aoa.rh
rownames(CLD.C.aoa.rh.ed) = sumData.C.aoa.rh$Label
CLD.C.aoa.rh.ed = rownames_to_column(CLD.C.aoa.rh.ed, var="Label")
CLD.C.aoa.rh.ed

CLD.C.aoa.rh.ed2 <- CLD.C.aoa.rh.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.aoa.rh.C.ed2$Treatment <- as.factor(tab.distance.aoa.rh.C.ed2$Treatment)

# Plot
dist.rh.aoa.CAP <- ggplot(tab.distance.aoa.rh.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7","#FF618C","#FF8EAC","#FFBBCD","#E69F00","#EDBA48","#F4D591"))+
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_x_discrete(labels = c("Between","Control","Drought",
                              "Between","Control","Drought",
                              "Between","Control","Drought"))+
  geom_text(data=CLD.C.aoa.rh.ed2,aes(x=Label, y = sumData.C.aoa.rh$Max + 0.04, label=Letter), vjust=0) +
  ylab("Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.rh.aoa.CAP
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_rhizo_CAP.tiff",
       dist.rh.aoa.CAP, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)






























































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

########## TEST

aoa.rh_dist_bc
#Sample group
item_groups <- sample_data(aoa.physeq_rh1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.rh_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.D <- subset(tab.distance, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.D$Label <- factor(tab.distance.D$Label)
str(tab.distance.D)
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.D)
str(tab.distance.D)
tab.distance.D
# Post Hoc Dunn Test
set.seed(1333)
dunn.D=dunnTest(Distance ~ Label,
                data = tab.distance.D, method = "bh")
library(rcompanion)
CLD.d = cldList(P.adj ~ Comparison, data=dunn.D$res)
CLD.d
# Plot
sumData.D <- ddply(tab.distance.D, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.d) = sumData.D$Label
CLD.d = rownames_to_column(CLD.d, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

dist.bulk.d<- ggplot(tab.distance.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.6,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#B6E3D7"),
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.d,aes(x=Label, y = sumData.D$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), 
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.bulk.d

####

aoa.cap.rh.dist
#Sample group
item_groups <- sample_data(aoa.physeq_rh1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.rh.dist, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.M <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.M 
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.M)
str(tab.distance.M)
tab.distance.M
# Post Hoc Dunn Test
set.seed(1333)
dunn.M=dunnTest(Distance ~ Label,
                data = tab.distance.M, method = "bh")
library(rcompanion)
CLD.m = cldList(P.adj ~ Comparison, data=dunn.M$res)
CLD.m
# Plot
sumData.M <- ddply(tab.distance.M, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.m) = sumData.M$Label
CLD.m = rownames_to_column(CLD.m, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)

dist.rh.m<- ggplot(tab.distance.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.6,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#E69F00","#EDBA48","#F4D591"),
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.m,aes(x=Label, y = sumData.M$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), 
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.rh.m

#####

aoa.cap.rh.dist
#Sample group
item_groups <- sample_data(aoa.physeq_rh1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.cap.rh.dist, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.K <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.K 
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = tab.distance.K)
str(tab.distance.K)
tab.distance.K
# Post Hoc Dunn Test
set.seed(1333)
dunn.K=dunnTest(Distance ~ Label,
                data = tab.distance.K, method = "bh")
library(rcompanion)
CLD.k = cldList(P.adj ~ Comparison, data=dunn.K$res)
CLD.k
# Plot
sumData.K <- ddply(tab.distance.K, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.k) = sumData.K$Label
CLD.k = rownames_to_column(CLD.k, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)

dist.bulk.k<- ggplot(tab.distance.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.6,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#FF618C","#FF8EAC","#FFBBCD"),
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  #geom_jitter(width=0.15) +
  geom_text(data=CLD.k,aes(x=Label, y = sumData.K$Max + 0.04, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  #geom_textdata=CLD.k,aes(label = cld, y = w + sd), vjust = -0.5)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), 
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.bulk.k

############################################################################################################################

#Distance matrix

# 1. AOA Bulk Soil - CAP Distance - CONFYM (K)
dist_matrix.aoa
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(dist_matrix.aoa, item_groups)
d.calcul
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.K <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.K$Label <- factor(tab.distance.K$Label)
str(tab.distance.K) 
# One-way ANOVA
set.seed(13)
K.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.K)
summary(K.dist.cap.aov) # significant, p-val < 0.0001, f-val= 939.1
# Post-Hoc Test
K.dist.cap.tuk <- TukeyHSD(K.dist.cap.aov)
Tukey.K <- as.data.frame(K.dist.cap.tuk$Label)
Tukey.K <- rownames_to_column(Tukey.K, var = "Comparison")
colnames(Tukey.K)[5] <- "p.adj"
# Make the significance letter
CLD.k = cldList(p.adj ~ Comparison, data=Tukey.K)
CLD.k
# re-order
CLD.k.ed = CLD.k[c(3,1,2),] # re-order
CLD.k.ed
# Plot
sumData.K <- ddply(tab.distance.K, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
sumData.K
rownames(CLD.k.ed) = sumData.K$Label
CLD.k.ed = rownames_to_column(CLD.k.ed, var="Label")
CLD.k.ed
# find color gradient
colfunc <- colorRampPalette(c("#FF618C", "white"))
colfunc(8)
col.k <- c("#FF618C","#FF8EAC","#FFBBCD")

dist.CAP.bulk.k<- ggplot(tab.distance.K, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.k ,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  ylim(0,15)+
  geom_text(data=CLD.k.ed,aes(x=Label, y = sumData.K$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5,size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  #annotate("text",x=0.5,y=1.10,label= "P-value = 0.003", hjust = 0, size = 6, fontface='italic') +
  ggtitle("CONFYM (K)")
dist.CAP.bulk.k

# 2. AOA Bulk Soil - CAP Distance - CONMIN
dist_matrix.aoa
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(dist_matrix.aoa, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.M <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Within cont.M","Within rain.M"))
tab.distance.M 
# One-way ANOVA
set.seed(13)
M.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.M)
summary(M.dist.cap.aov) # not significant, p-val < 0.0001, f-val= 452.1
# Post-Hoc Test
M.dist.cap.tuk <- TukeyHSD(M.dist.cap.aov)
Tukey.M <- as.data.frame(M.dist.cap.tuk$Label)
Tukey.M <- rownames_to_column(Tukey.M, var = "Comparison")
colnames(Tukey.M)[5] <- "p.adj"
# Make the significance letter
CLD.m = cldList(p.adj ~ Comparison, data=Tukey.M)
CLD.m
# re-order
CLD.m.ed = CLD.m[c(3,1,2),] # re-order
CLD.m.ed
# Plot
sumData.M <- ddply(tab.distance.M, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.m.ed) = sumData.M$Label
CLD.m.ed = rownames_to_column(CLD.m.ed, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#E69F00", "white"))
colfunc(8)
col.m <- c("#E69F00","#EDBA48","#F4D591")

dist.CAP.bulk.m<- ggplot(tab.distance.M, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = col.m,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.m.ed,aes(x=Label, y = sumData.M$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Bray-Curtis Distances") +
  ylim(0,15)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  #annotate("text",x=0.5,y=1.10,label= "P-value = 0.129", hjust = 0, size = 6, fontface='italic')+
  ggtitle("CONMIN (M)")
dist.CAP.bulk.m

# 3. AOA Bulk Soil - Bray Curtis Distance - BIODYN (D)
dist_matrix.aoa
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(dist_matrix.aoa, item_groups)
#Control
tab.distance = as_tibble(d.calcul) 
tab.distance.D <- subset(tab.distance, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.D$Label <- factor(tab.distance.D$Label)
tab.distance.D$Label
# One-way ANOVA
set.seed(13)
D.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.D)
summary(D.dist.cap.aov) # significant, p-val < 0.0001 f-val= 1813
# Post-Hoc Test
D.dist.cap.tuk <- TukeyHSD(D.dist.cap.aov)
D.dist.cap.tuk
Tukey.D <- as.data.frame(D.dist.cap.tuk$Label)
Tukey.D <- rownames_to_column(Tukey.D, var = "Comparison")
colnames(Tukey.D)[5] <- "p.adj"
# Make the significance letter
CLD.d = cldList(p.adj ~ Comparison, data=Tukey.D)
CLD.d
# re-order
CLD.d.ed = CLD.d[c(3,1,2),] # re-order
CLD.d.ed
# Plot
sumData.D <- ddply(tab.distance.D, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.d.ed) = sumData.D$Label
CLD.d.ed = rownames_to_column(CLD.d.ed, var="Label")
# find color gradient
colfunc <- colorRampPalette(c("#009E73", "white"))
colfunc(8)

col.d <- c("#009E73","#6DC7AF","#B6E3D7")

dist.CAP.bulk.d<- ggplot(tab.distance.D, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label),outlier.colour = NULL) + 
  scale_fill_manual(values = col.d,
                    labels = c("between", "within control", "within drought")) +
  scale_x_discrete(labels=c("Between", "Control", "Drought"))+
  geom_text(data=CLD.d.ed,aes(x=Label, y = sumData.D$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  ylim(0,15)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  #annotate("text",x=0.5,y=1.10,label= "P-value < 0.001", hjust = 0, size = 6, fontface='italic')+
  ggtitle("BIODYN (D)")
dist.CAP.bulk.d
# Combine figures
library(patchwork)
dist.CAP.AOA <- dist.CAP.bulk.d | dist.CAP.bulk.k | dist.CAP.bulk.m
dist.CAP.AOA
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist.NEWCAP.tiff",
       dist.CAP.AOA, device = "tiff",
       width = 9, height =4.5, 
       units= "in", dpi = 600)

#_______________________________________________________________________________

# AOA distance - everything in one plot

# 1. AOA Bulk Soil - Bray Curtis Distance 
dist_matrix.aoa
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(dist_matrix.aoa, item_groups)
#Control
tab.distance = as_tibble(d.calcul) 
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Label
str(tab.distance.C)
tab.distance.C$Group2 <- as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(tab.distance.C.ed)
tab.distance.C.ed2 <- tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D","Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K","Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(tab.distance.C.ed2)
tab.distance.C.ed2$Label <- factor(tab.distance.C.ed2$Label)
# One-way ANOVA
set.seed(13)
all.dist.cap.aov <- aov(Distance ~ Label, data = tab.distance.C.ed2)
summary(all.dist.cap.aov) # significant, p-val < 0.0001 f-val= 248.8
# Post-Hoc Test
all.dist.cap.tuk <- TukeyHSD(all.dist.cap.aov)
all.dist.cap.tuk
Tukey.all <- as.data.frame(all.dist.cap.tuk$Label)
Tukey.all <- rownames_to_column(Tukey.all, var = "Comparison")
colnames(Tukey.all)[5] <- "p.adj"
# Make the significance letter
CLD.all = cldList(p.adj ~ Comparison, data=Tukey.all)
CLD.all
# re-order
CLD.all.ed = CLD.all[c(3,1,2),] # re-order
CLD.all.ed
# Plot
sumData.all <- ddply(tab.distance.C.ed2, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(CLD.all.ed) = sumData.all$Label
CLD.all.ed = rownames_to_column(CLD.all.ed, var="Label")
CLD.all.ed

CLD.all.ed2 <- CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
dist.bulk.all<- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"),
  labels = c("Between cont.D and rain.D","Between cont.K and rain.K","Between cont.M and rain.M"))+
  # "Between cont.K and rain.K","Within cont.K","Within rain.K",
  #"Between cont.M and rain.M","Within cont.M","Within rain.M")) +
  facet_wrap(~ Treatment,scales = "free_x")+
  #scale_x_discrete(labels = c("Between cont.D and rain.D","Between cont.K and rain.K","Between cont.M and rain.M"))+
  geom_text(data=CLD.all.ed2,aes(x=Label, y = sumData.all$Max + 1, label=Letter), vjust=0) +
  ylab("Bray-Curtis Distances") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),#angle=45,hjust=1
        axis.title.y=element_text(size=13,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
dist.bulk.all
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_bulk_bray.tiff",
       dist.bulk.c, device = "tiff",
       width = 5.6, height =3, 
       units= "in", dpi = 600)

###########################################################################################################################






