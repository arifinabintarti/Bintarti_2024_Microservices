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
library(usedist)
library(PMCMRplus)
library(gdata)
BCstep1 = distance(aoa.physeq_bulk1,  "wunifrac")
BCstep1
d <- as.dist(BCstep1)
#BCstep1 <- as.matrix(BCstep1)
#write.table(BCstep1, file="weighted UniFrac distance matrix.csv")

item_groups <- sample_data(aoa.physeq_bulk1)
item_groups <- item_groups$x
#calculate dist between groups
d.calcul <- dist_groups(d, item_groups)


#Control
tab.distance = as_tibble(d.calcul) 
#tab.distance$Label <- gsub('\\s+', '', tab.distance$Label)


tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D", "Between cont.K and rain.K"))



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

treatOrderC <- factor(c("Between cont.M and rain.M","Between cont.D and rain.D", "Between cont.K and rain.K"))

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
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_text(data=tuckeyGroups,aes(x=Label, y = sumData$Mean + sumData$Sd, label=Groups), vjust=-1) +
  
  ylab("Weighted Unifrac Distances") +
  theme_bw() + 
  
  theme(plot.title = element_text(hjust = 0.5))
