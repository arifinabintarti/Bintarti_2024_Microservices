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
#Distance matrix

# 1. AOA Bulk Soil - Bray Curtis Distance - CONFYM
aoa.bulk_dist_bc
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.bulk_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.K <- subset(tab.distance, Label%in% c("Between cont.K and rain.K","Within cont.K","Within rain.K"))
tab.distance.K$Label <- factor(tab.distance.K$Label)
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


# 2. AOA Bulk Soil - Bray Curtis Distance - CONMIN
aoa.bulk_dist_bc
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.bulk_dist_bc, item_groups)
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

dist.bulk.m<- ggplot(tab.distance.M, aes(x = Label, y = Distance)) + 
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
dist.bulk.m

# 3. AOA Bulk Soil - Bray Curtis Distance - BIODYN
aoa.bulk_dist_bc
#Sample group
item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(aoa.bulk_dist_bc, item_groups)
#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
tab.distance.D <- subset(tab.distance, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D"))
tab.distance.D$Label <- factor(tab.distance.D$Label)
tab.distance.D$Label
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


