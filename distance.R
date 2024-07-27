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

#############################################################################################################################################

# AOA distance - everything in one plot

# 1. AOA Bulk Soil - CAP Distance using LDA results

# calculate distance from the CAP analysis
#dist_matrix.aoa <- dist(aoa.cap.bulk$x)
dist_matrix.aoa <- dist(aoa.cap.bulk$x[,c(1,2)])
#Sample group
aoa.BS.item_groups <- sample_data(aoa.physeq_bulk1)
aoa.BS.item_groups <- aoa.BS.item_groups$x
#calculate dist between groups
aoa.BS.d.calcul <- dist_groups(dist_matrix.aoa, aoa.BS.item_groups)
#Control
aoa.BS.tab.distance = as_tibble(aoa.BS.d.calcul) 
aoa.BS.tab.distance.C.all <- subset(aoa.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
aoa.BS.tab.distance.C <- subset(aoa.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
aoa.BS.tab.distance.C$Label <- factor(aoa.BS.tab.distance.C$Label)
aoa.BS.tab.distance.C$Label
str(aoa.BS.tab.distance.C)
aoa.BS.tab.distance.C$Group2 <- as.character(aoa.BS.tab.distance.C$Group2) 
aoa.BS.tab.distance.C.ed <- aoa.BS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
str(aoa.BS.tab.distance.C.ed)
aoa.BS.tab.distance.C.ed2 <- aoa.BS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
str(aoa.BS.tab.distance.C.ed2)
aoa.BS.tab.distance.C.ed2$Label <- factor(aoa.BS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aoa.BS.tab.distance.C.ed2.out <- aoa.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
aoa.BS.tab.distance.C.ed2.SW <- aoa.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(aoa.BS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.BS.tab.distance.C.ed2.Lave <- aoa.BS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)
aoa.BS.tab.distance.C.ed2.Lave

### LAVENE TEST FOR HOMOGENEITY IS VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.BS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 445.02, df = 2, p-value < 2.2e-16
# Post Hoc Dunn Test
set.seed(1333)
aoa.BS.dunn <- dunnTest(Distance ~ Label,
                data = aoa.BS.tab.distance.C.ed2, method = "bh")
aoa.BS.dunn

# Make the significance letter
aoa.BS.CLD.all = cldList(P.adj ~ Comparison, data=aoa.BS.dunn$res)
aoa.BS.CLD.all
# Plot
aoa.BS.sumData.all <- ddply(aoa.BS.tab.distance.C.ed2, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(aoa.BS.CLD.all) = aoa.BS.sumData.all$Label
aoa.BS.CLD.all.ed = rownames_to_column(aoa.BS.CLD.all, var="Label")
aoa.BS.CLD.all.ed

aoa.BS.CLD.all.ed2 <- aoa.BS.CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

aoa.BS.tab.distance.C.ed2$Treatment <- as.factor(aoa.BS.tab.distance.C.ed2$Treatment)
# Plot
AOA.BS.CAP.dist<- ggplot(aoa.BS.tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  scale_x_discrete(labels = c("BIOD","CONFY","CONM"))+
  geom_text(data=aoa.BS.CLD.all.ed2,aes(x=Label, y = aoa.BS.sumData.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Between treatment distance") +
  ylim(0,11)+
  theme_classic() +
  labs(subtitle = "C. AOA")+
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
AOA.BS.CAP.dist
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_dist_bulk_CAP2.tiff",
       dist.bulk.all, device = "tiff",
       width = 6.5, height =4, 
       units= "in", dpi = 600)
#################################################################################

# With all data between and within treatment

aoa.BS.tab.distance.C.all <- subset(aoa.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
aoa.BS.tab.distance.C.all$Label <- factor(aoa.BS.tab.distance.C.all$Label)
aoa.BS.tab.distance.C.all$Label
view(aoa.BS.tab.distance.C.all)
aoa.BS.tab.distance.C.all$Group2 <- as.character(aoa.BS.tab.distance.C.all$Group2) 
aoa.BS.tab.distance.C.all.ed <- aoa.BS.tab.distance.C.all %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN",
    endsWith(Group2, "K") ~ "CONFYM",
    endsWith(Group2, "M") ~ "CONMIN"))
str(aoa.BS.tab.distance.C.all.ed)
aoa.BS.tab.distance.C.all.ed2 <- aoa.BS.tab.distance.C.all.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D", "Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K", "Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(aoa.BS.tab.distance.C.all.ed2)
aoa.BS.tab.distance.C.all.ed2$Label <- factor(aoa.BS.tab.distance.C.all.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aoa.BS.tab.distance.C.all.ed2.out <- aoa.BS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
aoa.BS.tab.distance.C.all.ed2.out
# Saphiro-Wilk for normality
aoa.BS.tab.distance.C.all.ed2.SW <- aoa.BS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
aoa.BS.tab.distance.C.all.ed2.SW
ggqqplot(aoa.BS.tab.distance.C.all.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data

# Lavene test
aoa.BS.tab.distance.C.all.ed2.Lave <- aoa.BS.tab.distance.C.all.ed2 %>%
  levene_test(Distance ~ Label)
aoa.BS.tab.distance.C.all.ed2.Lave

### LAVENE TEST FOR HOMOGENEITY IS VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.BS.tab.distance.C.all.ed2) #Kruskal-Wallis chi-squared = 996.42, df = 8, p-value < 2.2e-16

aoa.BS.tab.distance.C.all.ed2.D <- subset(aoa.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.D and rain.D","Within cont.D", "Within rain.D"))
View(aoa.BS.tab.distance.C.all.ed2.D)
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.BS.tab.distance.C.all.ed2.D)

aoa.BS.tab.distance.C.all.ed2.K <- subset(aoa.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.K and rain.K","Within cont.K", "Within rain.K"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.BS.tab.distance.C.all.ed2.K)

aoa.BS.tab.distance.C.all.ed2.M <- subset(aoa.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.M and rain.M","Within cont.M", "Within rain.M"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.BS.tab.distance.C.all.ed2.M)
# Post Hoc Dunn Test
set.seed(1333)
aoa.BS.dunn.all <- dunnTest(Distance ~ Label,
                data = aoa.BS.tab.distance.C.all.ed2, method = "bh")
aoa.BS.dunn.all

set.seed(1333)
aoa.BS.dunn.D <- dunnTest(Distance ~ Label,
                data = aoa.BS.tab.distance.C.all.ed2.D, method = "bh")
aoa.BS.dunn.D

set.seed(1333)
aoa.BS.dunn.K <- dunnTest(Distance ~ Label,
                data = aoa.BS.tab.distance.C.all.ed2.K, method = "bh")
aoa.BS.dunn.K

set.seed(1333)
aoa.BS.dunn.M <- dunnTest(Distance ~ Label,
                data = aoa.BS.tab.distance.C.all.ed2.M, method = "bh")
aoa.BS.dunn.M

# Make the significance letter
aoa.BS.CLD.all.all = cldList(P.adj ~ Comparison, data=aoa.BS.dunn.all$res)
aoa.BS.CLD.all.all

aoa.BS.CLD.all.D = cldList(P.adj ~ Comparison, data=aoa.BS.dunn.D$res)
aoa.BS.CLD.all.D

aoa.BS.CLD.all.K = cldList(P.adj ~ Comparison, data=aoa.BS.dunn.K$res)
aoa.BS.CLD.all.K

aoa.BS.CLD.all.M = cldList(P.adj ~ Comparison, data=aoa.BS.dunn.M$res)
aoa.BS.CLD.all.M

aoa.BS.CLD.all.all.ord <- aoa.BS.CLD.all.all[c(1,4,7,2,5,8,3,6,9),] # re-order
aoa.BS.CLD.all.all.ord

# Plot
aoa.BS.sumData.all.all <- ddply(aoa.BS.tab.distance.C.all.ed2, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(aoa.BS.CLD.all.all.ord) = aoa.BS.sumData.all.all$Label
aoa.BS.CLD.all.all.ed = rownames_to_column(aoa.BS.CLD.all.all.ord, var="Label")
aoa.BS.CLD.all.all.ed

aoa.BS.CLD.all.all.ed2 <- aoa.BS.CLD.all.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN",
    endsWith(Group, "K") ~ "CONFYM",
    endsWith(Group, "M") ~ "CONMIN"))

aoa.BS.tab.distance.C.all.ed2$Treatment <- as.factor(aoa.BS.tab.distance.C.all.ed2$Treatment)
# Plot
AOA.BS.CAP.dist.all <- ggplot(aoa.BS.tab.distance.C.all.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_fill_manual(values = c("#009E73","#6DC7AF","#DAF1EB","#FF618C","#FF8EAC","#FFE8EE","#E69F00","#EDBA48","#FBF1DA"))+
  scale_x_discrete(labels = c("Between","within \n control","within \n drought","Between","within \n control","within \n drought","Between","within \n control","within \n drought"))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  #scale_x_discrete(labels = c("BIOD","CONFY","CONM"))+
  geom_text(data=aoa.BS.CLD.all.all.ed2,aes(x=Label, y = aoa.BS.sumData.all.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Distance") +
  ylim(0,11)+
  theme_bw() +
  labs(subtitle = "C. AOA")+
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        strip.text = element_text(size=18),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
AOA.BS.CAP.dist.all
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("AOA_distance_CAP1et2_betweenetwithin_BS.tiff",
       AOA.BS.CAP.dist.all, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600, compression="lzw",bg = 'white')



#_______________________________________________________________________________
# 2. AOA Rhizosphere- CAP Distance using LDA Results

# calculate distance from the CAP analysis
#dist_matrix.aoa.rh <- dist(aoa.cap.rh$x)
dist_matrix.aoa.rh <- dist(aoa.cap.rh$x[,c(1,2)])
dist_matrix.aoa.rh
#Sample group
aoa.RS.item_groups <- sample_data(aoa.physeq_rh1)
aoa.RS.item_groups <- aoa.RS.item_groups$x
#calculate dist between groups
aoa.RS.d.calcul <- dist_groups(dist_matrix.aoa.rh, aoa.RS.item_groups)
#Control
aoa.RS.tab.distance = as_tibble(aoa.RS.d.calcul) 
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
aoa.RS.tab.distance.C <- subset(aoa.RS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
aoa.RS.tab.distance.C$Label <- factor(aoa.RS.tab.distance.C$Label)
aoa.RS.tab.distance.C$Label
str(aoa.RS.tab.distance.C)
aoa.RS.tab.distance.C$Group2 <- as.character(aoa.RS.tab.distance.C$Group2) 
aoa.RS.tab.distance.C.ed <- aoa.RS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN",
    endsWith(Group2, "K") ~ "CONFYM",
    endsWith(Group2, "M") ~ "CONMIN"))
str(aoa.RS.tab.distance.C.ed)
aoa.RS.tab.distance.C.ed2 <- aoa.RS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
str(aoa.RS.tab.distance.C.ed2)
aoa.RS.tab.distance.C.ed2$Label <- factor(aoa.RS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aoa.RS.tab.distance.C.ed2.out <- aoa.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
aoa.RS.tab.distance.C.ed2.SW <- aoa.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(aoa.RS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.RS.tab.distance.C.ed2.Lave <- aoa.RS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)

### SAPHIRO-WILK TEST FOR NORMALITY IS VIOLATED --> Use Kruskal-Wallis

#_______________________________________________________________________________
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
#_______________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.RS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 100.2, df = 2, p-value < 2.2e-16
# Post Hoc Dunn Test
set.seed(1333)
aoa.RS.dunn <- dunnTest(Distance ~ Label,
                        data = aoa.RS.tab.distance.C.ed2, method = "bh")
aoa.RS.dunn
#_______________________________________________________________________________

#################################################################################

# With all data between and within treatment

aoa.RS.tab.distance.C.all <- subset(aoa.RS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
aoa.RS.tab.distance.C.all$Label <- factor(aoa.RS.tab.distance.C.all$Label)
aoa.RS.tab.distance.C.all$Label
view(aoa.RS.tab.distance.C.all)
aoa.RS.tab.distance.C.all$Group2 <- as.character(aoa.RS.tab.distance.C.all$Group2) 
aoa.RS.tab.distance.C.all.ed <- aoa.RS.tab.distance.C.all %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN",
    endsWith(Group2, "K") ~ "CONFYM",
    endsWith(Group2, "M") ~ "CONMIN"))
str(aoa.RS.tab.distance.C.all.ed)
aoa.RS.tab.distance.C.all.ed2 <- aoa.RS.tab.distance.C.all.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D", "Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K", "Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(aoa.RS.tab.distance.C.all.ed2)
aoa.RS.tab.distance.C.all.ed2$Label <- factor(aoa.RS.tab.distance.C.all.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aoa.RS.tab.distance.C.all.ed2.out <- aoa.RS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
aoa.RS.tab.distance.C.all.ed2.out
# Saphiro-Wilk for normality
aoa.RS.tab.distance.C.all.ed2.SW <- aoa.RS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
aoa.RS.tab.distance.C.all.ed2.SW
ggqqplot(aoa.RS.tab.distance.C.all.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data

# Lavene test
aoa.RS.tab.distance.C.all.ed2.Lave <- aoa.RS.tab.distance.C.all.ed2 %>%
  levene_test(Distance ~ Label)
aoa.RS.tab.distance.C.all.ed2.Lave

### LAVENE TEST FOR HOMOGENEITY IS VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.RS.tab.distance.C.all.ed2) #Kruskal-Wallis chi-squared = 414.13, df = 8, p-value < 2.2e-16

aoa.RS.tab.distance.C.all.ed2.D <- subset(aoa.RS.tab.distance.C.all.ed2, Label%in% c("Between cont.D and rain.D","Within cont.D", "Within rain.D"))
View(aoa.RS.tab.distance.C.all.ed2.D)
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.RS.tab.distance.C.all.ed2.D)

aoa.RS.tab.distance.C.all.ed2.K <- subset(aoa.RS.tab.distance.C.all.ed2, Label%in% c("Between cont.K and rain.K","Within cont.K", "Within rain.K"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.RS.tab.distance.C.all.ed2.K)

aoa.RS.tab.distance.C.all.ed2.M <- subset(aoa.RS.tab.distance.C.all.ed2, Label%in% c("Between cont.M and rain.M","Within cont.M", "Within rain.M"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aoa.RS.tab.distance.C.all.ed2.M)
# Post Hoc Dunn Test
set.seed(1333)
aoa.RS.dunn.all <- dunnTest(Distance ~ Label,
                data = aoa.RS.tab.distance.C.all.ed2, method = "bh")
aoa.RS.dunn.all

set.seed(1333)
aoa.RS.dunn.D <- dunnTest(Distance ~ Label,
                data = aoa.RS.tab.distance.C.all.ed2.D, method = "bh")
aoa.RS.dunn.D

set.seed(1333)
aoa.RS.dunn.K <- dunnTest(Distance ~ Label,
                data = aoa.RS.tab.distance.C.all.ed2.K, method = "bh")
aoa.RS.dunn.K

set.seed(1333)
aoa.RS.dunn.M <- dunnTest(Distance ~ Label,
                data = aoa.RS.tab.distance.C.all.ed2.M, method = "bh")
aoa.RS.dunn.M

# Make the significance letter
aoa.RS.CLD.all.all = cldList(P.adj ~ Comparison, data=aoa.RS.dunn.all$res)
aoa.RS.CLD.all.all

aoa.RS.CLD.all.D = cldList(P.adj ~ Comparison, data=aoa.RS.dunn.D$res)
aoa.RS.CLD.all.D

aoa.RS.CLD.all.K = cldList(P.adj ~ Comparison, data=aoa.RS.dunn.K$res)
aoa.RS.CLD.all.K

aoa.RS.CLD.all.M = cldList(P.adj ~ Comparison, data=aoa.RS.dunn.M$res)
aoa.RS.CLD.all.M

# Make the significance letter
aoa.RS.CLD.all.all = cldList(P.adj ~ Comparison, data=aoa.RS.dunn.all$res)
aoa.RS.CLD.all.all

aoa.RS.CLD.all.all.ord <- aoa.RS.CLD.all.all[c(1,4,7,2,5,8,3,6,9),] # re-order
aoa.RS.CLD.all.all.ord

# Plot
aoa.RS.sumData.all.all <- ddply(aoa.RS.tab.distance.C.all.ed2, "Label", summarise,
                     Max = max(Distance),
                     N    = length(Distance),
                     Mean = mean(Distance),
                     Sd   = sd(Distance),
                     Se   = Sd / sqrt(N))
rownames(aoa.RS.CLD.all.all.ord) = aoa.RS.sumData.all.all$Label
aoa.RS.CLD.all.all.ed = rownames_to_column(aoa.RS.CLD.all.all.ord, var="Label")
aoa.RS.CLD.all.all.ed

aoa.RS.CLD.all.all.ed2 <- aoa.RS.CLD.all.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN",
    endsWith(Group, "K") ~ "CONFYM",
    endsWith(Group, "M") ~ "CONMIN"))

aoa.RS.tab.distance.C.all.ed2$Treatment <- as.factor(aoa.RS.tab.distance.C.all.ed2$Treatment)

# Plot
AOA.RS.CAP.dist.all <- ggplot(aoa.RS.tab.distance.C.all.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_fill_manual(values = c("#009E73","#6DC7AF","#DAF1EB","#FF618C","#FF8EAC","#FFE8EE","#E69F00","#EDBA48","#FBF1DA"))+
  scale_x_discrete(labels = c("Between","within \n control","within \n drought","Between","within \n control","within \n drought","Between","within \n control","within \n drought"))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  #scale_x_discrete(labels = c("BIOD","CONFY","CONM"))+
  geom_text(data=aoa.RS.CLD.all.all.ed2,aes(x=Label, y = aoa.RS.sumData.all.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Distance") +
  ylim(0,11)+
  theme_bw() +
  labs(subtitle = "C. AOA")+
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        strip.text = element_text(size=18),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
AOA.RS.CAP.dist.all
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("AOA_distance_CAP1et2_betweenetwithin_RS.tiff",
       AOA.RS.CAP.dist.all, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600, compression="lzw",bg = 'white')



###########################################################################################################################

# AOB distance - everything in one plot

# 1. AOB Bulk Soil - CAP Distance using LDA results

# calculate distance from the CAP analysis

#dist_matrix.aob <- dist(aob.cap.bulk$x)
dist_matrix.aob <- dist(aob.cap.bulk$x[,c(1,2)])
dist_matrix.aob
#Sample group
aob.BS.item_groups <- sample_data(aob.physeq_bulk1)
aob.BS.item_groups <- aob.BS.item_groups$x
#calculate dist between groups
aob.BS.d.calcul <- dist_groups(dist_matrix.aob, aob.BS.item_groups)
#Control
aob.BS.tab.distance = as_tibble(aob.BS.d.calcul) 
aob.BS.tab.distance.C <- subset(aob.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
aob.BS.tab.distance.C$Label <- factor(aob.BS.tab.distance.C$Label)
aob.BS.tab.distance.C$Label
aob.BS.tab.distance.C$Group2 <- as.character(aob.BS.tab.distance.C$Group2) 
aob.BS.tab.distance.C.ed <- aob.BS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN",
    endsWith(Group2, "K") ~ "CONFYM",
    endsWith(Group2, "M") ~ "CONMIN"))
aob.BS.tab.distance.C.ed2 <- aob.BS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
aob.BS.tab.distance.C.ed2$Label <- factor(aob.BS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aob.BS.tab.distance.C.ed2.out <- aob.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
aob.BS.tab.distance.C.ed2.SW <- aob.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(aob.BS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.BS.tab.distance.C.ed2.Lave <- aob.BS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)

### LAVENE TEST FOR HOMOGENEITY IS VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.BS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 65.768, df = 2, p-value = 5.231e-15
# Post Hoc Dunn Test
set.seed(1333)
aob.BS.dunn <- dunnTest(Distance ~ Label,
                        data = aob.BS.tab.distance.C.ed2, method = "bh")
aob.BS.dunn

# Make the significance letter
aob.BS.CLD.all = cldList(P.adj ~ Comparison, data=aob.BS.dunn$res)
aob.BS.CLD.all
# Plot
aob.BS.sumData.all <- ddply(aob.BS.tab.distance.C.ed2, "Label", summarise,
                     Max = max(Distance),
                     N    = length(Distance),
                     Mean = mean(Distance),
                     Sd   = sd(Distance),
                     Se   = Sd / sqrt(N))
rownames(aob.BS.CLD.all) = aob.BS.sumData.all$Label
aob.BS.CLD.all.ed = rownames_to_column(aob.BS.CLD.all, var="Label")
aob.BS.CLD.all.ed

aob.BS.CLD.all.ed2 <- aob.BS.CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

aob.BS.tab.distance.C.ed2$Treatment <- as.factor(aob.BS.tab.distance.C.ed2$Treatment)
# Plot
AOB.BS.CAP.dist<- ggplot(aob.BS.tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  scale_x_discrete(labels = c("BIOD","CONF","CONM"))+
  geom_text(data = aob.BS.CLD.all.ed2,aes(x = Label, y = aob.BS.sumData.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Between treatment distance") + 
  ylim(0,11)+
  labs(subtitle="A. AOB", title="Bulk Soil")+
  theme_classic() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0, size = 25, face='bold'),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
AOB.BS.CAP.dist
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist_bulk_CAP2.tiff",
       aob.dist.bulk.all, device = "tiff",
       width = 6.5, height =4, 
       units= "in", dpi = 600)

#################################################################################

# With all data between and within treatment

aob.BS.tab.distance.C.all <- subset(aob.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K","Within cont.D","Within cont.M","Within cont.K", "Within rain.D","Within rain.M","Within rain.K"))
aob.BS.tab.distance.C.all$Label <- factor(aob.BS.tab.distance.C.all$Label)
aob.BS.tab.distance.C.all$Label
view(aob.BS.tab.distance.C.all)
aob.BS.tab.distance.C.all$Group2 <- as.character(aob.BS.tab.distance.C.all$Group2) 
aob.BS.tab.distance.C.all.ed <- aob.BS.tab.distance.C.all %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN",
    endsWith(Group2, "K") ~ "CONFYM",
    endsWith(Group2, "M") ~ "CONMIN"))
str(aob.BS.tab.distance.C.all.ed)
aob.BS.tab.distance.C.all.ed2 <- aob.BS.tab.distance.C.all.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D", "Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K", "Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M","Within cont.M","Within rain.M")))
str(aob.BS.tab.distance.C.all.ed2)
aob.BS.tab.distance.C.all.ed2$Label <- factor(aob.BS.tab.distance.C.all.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aob.BS.tab.distance.C.all.ed2.out <- aob.BS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
aob.BS.tab.distance.C.all.ed2.out
# Saphiro-Wilk for normality
aob.BS.tab.distance.C.all.ed2.SW <- aob.BS.tab.distance.C.all.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
aob.BS.tab.distance.C.all.ed2.SW
ggqqplot(aob.BS.tab.distance.C.all.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data

# Lavene test
aob.BS.tab.distance.C.all.ed2.Lave <- aob.BS.tab.distance.C.all.ed2 %>%
  levene_test(Distance ~ Label)
aob.BS.tab.distance.C.all.ed2.Lave

### LAVENE TEST FOR HOMOGENEITY IS VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.BS.tab.distance.C.all.ed2) #Kruskal-Wallis chi-squared = 176.89, df = 8, p-value < 2.2e-16

aob.BS.tab.distance.C.all.ed2.D <- subset(aob.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.D and rain.D","Within cont.D", "Within rain.D"))
View(aob.BS.tab.distance.C.all.ed2.D)
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.BS.tab.distance.C.all.ed2.D)

aob.BS.tab.distance.C.all.ed2.K <- subset(aob.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.K and rain.K","Within cont.K", "Within rain.K"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.BS.tab.distance.C.all.ed2.K)

aob.BS.tab.distance.C.all.ed2.M <- subset(aob.BS.tab.distance.C.all.ed2, Label%in% c("Between cont.M and rain.M","Within cont.M", "Within rain.M"))
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.BS.tab.distance.C.all.ed2.M)
# Post Hoc Dunn Test
set.seed(1333)
aob.BS.dunn.all <- dunnTest(Distance ~ Label,
                data = aob.BS.tab.distance.C.all.ed2, method = "bh")
aob.BS.dunn.all

set.seed(1333)
aob.BS.dunn.D <- dunnTest(Distance ~ Label,
                data = aob.BS.tab.distance.C.all.ed2.D, method = "bh")
aob.BS.dunn.D

set.seed(1333)
aob.BS.dunn.K <- dunnTest(Distance ~ Label,
                data = aob.BS.tab.distance.C.all.ed2.K, method = "bh")
aob.BS.dunn.K

set.seed(1333)
aob.BS.dunn.M <- dunnTest(Distance ~ Label,
                data = aob.BS.tab.distance.C.all.ed2.M, method = "bh")
aob.BS.dunn.M

# Make the significance letter
aob.BS.CLD.all.all = cldList(P.adj ~ Comparison, data=aob.BS.dunn.all$res)
aob.BS.CLD.all.all

aob.BS.CLD.all.D = cldList(P.adj ~ Comparison, data=aob.BS.dunn.D$res)
aob.BS.CLD.all.D

aob.BS.CLD.all.K = cldList(P.adj ~ Comparison, data=aob.BS.dunn.K$res)
aob.BS.CLD.all.K

aob.BS.CLD.all.M = cldList(P.adj ~ Comparison, data=aob.BS.dunn.M$res)
aob.BS.CLD.all.M

aob.BS.CLD.all.all.ord <- aob.BS.CLD.all.all[c(1,4,7,2,5,8,3,6,9),] # re-order
aob.BS.CLD.all.all.ord

# Plot
aob.BS.sumData.all.all <- ddply(aob.BS.tab.distance.C.all.ed2, "Label", summarise,
                   Max = max(Distance),
                   N    = length(Distance),
                   Mean = mean(Distance),
                   Sd   = sd(Distance),
                   Se   = Sd / sqrt(N))
rownames(aob.BS.CLD.all.all.ord) = aob.BS.sumData.all.all$Label
aob.BS.CLD.all.all.ed = rownames_to_column(aob.BS.CLD.all.all.ord, var="Label")
aob.BS.CLD.all.all.ed

aob.BS.CLD.all.all.ed2 <- aob.BS.CLD.all.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN",
    endsWith(Group, "K") ~ "CONFYM",
    endsWith(Group, "M") ~ "CONMIN"))

aob.BS.tab.distance.C.all.ed2$Treatment <- as.factor(aob.BS.tab.distance.C.all.ed2$Treatment)
# Plot
AOB.BS.CAP.dist.all <- ggplot(aob.BS.tab.distance.C.all.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  facet_wrap(~ Treatment,scales = "free_x")+
  scale_fill_manual(values = c("#009E73","#6DC7AF","#DAF1EB","#FF618C","#FF8EAC","#FFE8EE","#E69F00","#EDBA48","#FBF1DA"))+
  scale_x_discrete(labels = c("Between","within \n control","within \n drought","Between","within \n control","within \n drought","Between","within \n control","within \n drought"))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  #scale_x_discrete(labels = c("BIOD","CONFY","CONM"))+
  geom_text(data=aob.BS.CLD.all.all.ed2,aes(x=Label, y = aob.BS.sumData.all.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Distance") +
  ylim(0,11)+
  theme_bw() +
  labs(subtitle = "A. AOB")+
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        strip.text = element_text(size=18),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
  #annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
AOB.BS.CAP.dist.all
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("AOB_distance_CAP1et2_betweenetwithin_BS.tiff",
       AOB.BS.CAP.dist.all, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600, compression="lzw",bg = 'white')





#_______________________________________________________________________________
# 2. AOB Rhizosphere- CAP Distance using LDA Results

# calculate distance from the CAP analysis
dist_matrix.aob.rh <- dist(aob.cap.rh$x[,c(1,2)])
dist_matrix.aob.rh
#Sample group
aob.RS.item_groups <- sample_data(aob.physeq_rh1)
aob.RS.item_groups <- aob.RS.item_groups$x
#calculate dist between groups
aob.RS.d.calcul <- dist_groups(dist_matrix.aob.rh, aob.RS.item_groups)
#Control
aob.RS.tab.distance = as_tibble(aob.RS.d.calcul) 
aob.RS.tab.distance.C <- subset(aob.RS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
aob.RS.tab.distance.C$Label <- factor(aob.RS.tab.distance.C$Label)
aob.RS.tab.distance.C$Label
aob.RS.tab.distance.C$Group2 <- as.character(aob.RS.tab.distance.C$Group2) 
aob.RS.tab.distance.C.ed <- aob.RS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
aob.RS.tab.distance.C.ed2 <- aob.RS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
aob.RS.tab.distance.C.ed2$Label <- factor(aob.RS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
aob.RS.tab.distance.C.ed2.out <- aob.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
aob.RS.tab.distance.C.ed2.SW <- aob.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(aob.RS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.RS.tab.distance.C.ed2.Lave <- aob.RS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)

### NO ASSUMPTIONS ARE VIOLATED --> Use ANOVA

#_______________________________________________________________________________
# One-way ANOVA
set.seed(13)
aob.dist.cap.rh.aov <- aov(Distance ~ Label, data = aob.tab.distance.C.rh.ed2)
summary(aob.dist.cap.rh.aov) # significant, p-val < 0.0001 f-val= 248.8

# Post-Hoc Test
aob.dist.cap.rh.tuk <- TukeyHSD(aob.dist.cap.rh.aov)
aob.dist.cap.rh.tuk
Tukey.all <- as.data.frame(aob.dist.cap.rh.tuk$Label)
Tukey.all <- rownames_to_column(Tukey.all, var = "Comparison")
colnames(Tukey.all)[5] <- "p.adj"

#_____________________________________________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = aob.RS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 4.8072, df = 2, p-value = 0.09039
# Post Hoc Dunn Test
set.seed(1333)
aob.RS.dunn <- dunnTest(Distance ~ Label,
                        data = aob.RS.tab.distance.C.ed2, method = "bh")
aob.RS.dunn
#______________________________________________________________________________________________________________
# Make the significance letter
aob.RS.CLD.all = cldList(P.adj ~ Comparison, data=aob.RS.dunn$res)
aob.RS.CLD.all
# Plot
aob.RS.sumData.all <- ddply(aob.RS.tab.distance.C.ed2, "Label", summarise,
                            Max = max(Distance),
                            N    = length(Distance),
                            Mean = mean(Distance),
                            Sd   = sd(Distance),
                            Se   = Sd / sqrt(N))
rownames(aob.RS.CLD.all) = aob.RS.sumData.all$Label
aob.RS.CLD.all.ed = rownames_to_column(aob.RS.CLD.all, var="Label")
aob.RS.CLD.all.ed

aob.RS.CLD.all.ed2 <- aob.RS.CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

aob.RS.tab.distance.C.ed2$Treatment <- as.factor(aob.RS.tab.distance.C.ed2$Treatment)
# Plot
AOB.RS.CAP.dist<- ggplot(aob.RS.tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  scale_x_discrete(labels = c("BIOD","CONF","CONM"))+
  geom_text(data = aob.RS.CLD.all.ed2,aes(x = Label, y = aob.RS.sumData.all$Max + 1, label = Letter), vjust=0, size=6) +
  labs(title = "Rhizosphere", subtitle = "B. AOB")+
  ylab("Between treatment distance") +
  ylim(0,11)+
  theme_classic() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0, size = 25, face='bold'),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15),
        #axis.ticks.y = element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP=0.09", hjust = 0, size = 4.5, fontface='italic')
AOB.RS.CAP.dist
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_dist_rhizo_CAP2.tiff",
       aob.dist.rhizo.all, device = "tiff",
       width = 6.5, height =4, 
       units= "in", dpi = 600)


###########################################################################################################################

# COMAMMOX distance - everything in one plot

# 1. COMAMMOX Bulk Soil - CAP Distance using LDA results

# calculate distance from the CAP analysis
#dist_matrix.com <- dist(com.cap.bulk$x)
dist_matrix.com <- dist(com.cap.bulk$x[,c(1,2)])
dist_matrix.com
#Sample group
com.BS.item_groups <- sample_data(com.physeq_bulk1)
com.BS.item_groups <- com.BS.item_groups$x
#calculate dist between groups
com.BS.d.calcul <- dist_groups(dist_matrix.com, com.BS.item_groups)
#Control
com.BS.tab.distance = as_tibble(com.BS.d.calcul) 
com.BS.tab.distance.C <- subset(com.BS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
com.BS.tab.distance.C$Label <- factor(com.BS.tab.distance.C$Label)
com.BS.tab.distance.C$Label
com.BS.tab.distance.C$Group2 <- as.character(com.BS.tab.distance.C$Group2) 
com.BS.tab.distance.C.ed <- com.BS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
com.BS.tab.distance.C.ed2 <- com.BS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
com.BS.tab.distance.C.ed2$Label <- factor(com.BS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
com.BS.tab.distance.C.ed2.out <- com.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
com.BS.tab.distance.C.ed2.SW <- com.BS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(com.BS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.BS.tab.distance.C.ed2.Lave <- com.BS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)

### SAPHIRO-WILK TEST FOR NORMALITY AND LAVENE TEST FOR HOMOGENEITY ARE VIOLATED --> Use Kruskal-Wallis

# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = com.BS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 423.38, df = 2, p-value < 2.2e-16
# Post Hoc Dunn Test
set.seed(1333)
com.BS.dunn <- dunnTest(Distance ~ Label,
                        data = com.BS.tab.distance.C.ed2, method = "bh")
com.BS.dunn

# Make the significance letter
com.BS.CLD.all = cldList(P.adj ~ Comparison, data=com.BS.dunn$res)
com.BS.CLD.all
# Plot
com.BS.sumData.all <- ddply(com.BS.tab.distance.C.ed2, "Label", summarise,
                            Max = max(Distance),
                            N    = length(Distance),
                            Mean = mean(Distance),
                            Sd   = sd(Distance),
                            Se   = Sd / sqrt(N))
rownames(com.BS.CLD.all) = com.BS.sumData.all$Label
com.BS.CLD.all.ed = rownames_to_column(com.BS.CLD.all, var="Label")
com.BS.CLD.all.ed

com.BS.CLD.all.ed2 <- com.BS.CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

com.BS.tab.distance.C.ed2$Treatment <- as.factor(com.BS.tab.distance.C.ed2$Treatment)
# Plot
COM.BS.CAP.dist<- ggplot(com.BS.tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  scale_x_discrete(labels = c("BIOD","CONF","CONM"))+
  geom_text(data = com.BS.CLD.all.ed2,aes(x = Label, y = com.BS.sumData.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Between treatment distance") +
  labs(subtitle = "E. Comammox")+
  ylim(0,11)+
  theme_classic() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
COM.BS.CAP.dist
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist_bulk_CAP2.tiff",
       com.dist.bulk.all, device = "tiff",
       width = 6.5, height =4, 
       units= "in", dpi = 600)

#_______________________________________________________________________________
# 2. COMAMMOX Rhizosphere- CAP Distance using LDA Results

# calculate distance from the CAP analysis
#dist_matrix.com.rh <- dist(com.cap.rh$x)
dist_matrix.com.rh <- dist(com.cap.rh$x[,c(1,2)])
#Sample group
com.RS.item_groups <- sample_data(com.physeq_rh1)
com.RS.item_groups <- com.RS.item_groups$x
#calculate dist between groups
com.RS.d.calcul <- dist_groups(dist_matrix.com.rh, com.RS.item_groups)
#Control
com.RS.tab.distance = as_tibble(com.RS.d.calcul) 
com.RS.tab.distance.C <- subset(com.RS.tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
com.RS.tab.distance.C$Label <- factor(com.RS.tab.distance.C$Label)
com.RS.tab.distance.C$Label
com.RS.tab.distance.C$Group2 <- as.character(com.RS.tab.distance.C$Group2) 
com.RS.tab.distance.C.ed <- com.RS.tab.distance.C %>%
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
com.RS.tab.distance.C.ed2 <- com.RS.tab.distance.C.ed %>% 
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Between cont.K and rain.K",
                                   "Between cont.M and rain.M")))
com.RS.tab.distance.C.ed2$Label <- factor(com.RS.tab.distance.C.ed2$Label)

# Check ANOVA assumptions

# check assumption (outliers)
com.RS.tab.distance.C.ed2.out <- com.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  identify_outliers(Distance) # no extreme outliers
# Saphiro-Wilk for normality
com.RS.tab.distance.C.ed2.SW <- com.RS.tab.distance.C.ed2 %>%
  group_by(Label) %>%
  shapiro_test(Distance)
ggqqplot(com.RS.tab.distance.C.ed2, "Distance", ggtheme = theme_bw())#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.RS.tab.distance.C.ed2.Lave <- com.RS.tab.distance.C.ed2 %>%
  levene_test(Distance ~ Label)

### SAPHIRO-WILKS TEST FOR NORMALITY IS VIOLATED BUT ONLY ONE POINT --> Use ANOVA or KW

#_______________________________________________________________________________
# One-way ANOVA
set.seed(13)
com.RS.dist.cap.aov <- aov(Distance ~ Label, data = com.tab.distance.C.rh.ed2)
summary(com.dist.cap.rh.aov) # significant 3.6e-11 ***

# Post-Hoc Test
com.dist.cap.rh.tuk <- TukeyHSD(com.dist.cap.rh.aov)
com.dist.cap.rh.tuk
Tukey.all <- as.data.frame(com.dist.cap.rh.tuk$Label)
Tukey.all <- rownames_to_column(Tukey.all, var = "Comparison")
colnames(Tukey.all)[5] <- "p.adj"

#_____________________________________________________________________________________________________________
# Kruskal-Wallis Test
set.seed(1333)
kruskal.test(Distance ~ Label,
             data = com.RS.tab.distance.C.ed2) #Kruskal-Wallis chi-squared = 116.76, df = 2, p-value < 2.2e-16
# Post Hoc Dunn Test
set.seed(1333)
com.RS.dunn <- dunnTest(Distance ~ Label,
                        data = com.RS.tab.distance.C.ed2, method = "bh")
com.RS.dunn
#______________________________________________________________________________________________________________
# Make the significance letter
com.RS.CLD.all = cldList(P.adj ~ Comparison, data=com.RS.dunn$res)
com.RS.CLD.all
# Plot
com.RS.sumData.all <- ddply(com.RS.tab.distance.C.ed2, "Label", summarise,
                            Max = max(Distance),
                            N    = length(Distance),
                            Mean = mean(Distance),
                            Sd   = sd(Distance),
                            Se   = Sd / sqrt(N))
rownames(com.RS.CLD.all) = com.RS.sumData.all$Label
com.RS.CLD.all.ed = rownames_to_column(com.RS.CLD.all, var="Label")
com.RS.CLD.all.ed

com.RS.CLD.all.ed2 <- com.RS.CLD.all.ed %>%
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

com.RS.tab.distance.C.ed2$Treatment <- as.factor(com.RS.tab.distance.C.ed2$Treatment)
# Plot
COM.RS.CAP.dist<- ggplot(com.RS.tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  theme_classic()+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_x_discrete(labels = c("control vs drought \nBIODYN","control vs drought \nCONFYM","control vs drought \nCONMIN"))+
  scale_x_discrete(labels = c("BIOD","CONF","CONM"))+
  geom_text(data = com.RS.CLD.all.ed2,aes(x = Label, y = com.RS.sumData.all$Max + 1, label = Letter), vjust=0, size=6) +
  ylab("Between treatment distance") +
  labs(subtitle = "F. Comammox")+
  ylim(0,11)+
  #theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.subtitle = element_text(hjust = 0, size = 20, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15), 
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.title.y=element_text(size=16),
        axis.title=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))+
  annotate("text",x=0.5,y=11,label= "Kruskal-Wallis,\nP< 0.0001", hjust = 0, size = 4.5, fontface='italic')
  #annotate("text",x=0.5,y=11,label= "Kruskal-Wallis, P< 0.0001", hjust = 0, size = 4.5, fontface='italic')
COM.RS.CAP.dist
setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("COM_dist_rhizo_CAP2.tiff",
       com.dist.rhizo.all, device = "tiff",
       width = 6.5, height =4, 
       units= "in", dpi = 600)

# Combine All figures

CAP_dist.All <- (AOB.BS.CAP.dist / AOA.BS.CAP.dist / COM.BS.CAP.dist) | (AOB.RS.CAP.dist / AOA.RS.CAP.dist / COM.RS.CAP.dist)
CAP_dist.All
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Supp.Fig.2dpi300.tiff",
       CAP_dist.All, device = "tiff",
       width = 12.6, height = 14.2, 
       units= "in", dpi = 300, compression="lzw")
CAP.BS.Dist <- (AOB.BS.CAP.dist / AOA.BS.CAP.dist / COM.BS.CAP.dist)
ggsave("Supp.Fig.2CAP.BS.Dist.tiff",
       CAP.BS.Dist, device = "tiff",
       width = 2, height = 11, 
       units= "in", dpi = 300, compression="lzw")

CAP.RS.Dist <- (AOB.RS.CAP.dist / AOA.RS.CAP.dist / COM.RS.CAP.dist)
ggsave("Supp.Fig.2CAP.RS.Dist.tiff",
       CAP.RS.Dist, device = "tiff",
       width =1.8, height = 11, 
       units= "in", dpi = 300, compression="lzw")


