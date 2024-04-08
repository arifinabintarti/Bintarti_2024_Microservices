#install.packages("parallel")
#install.packages("ggforce")
#install.packages("BiodiversityR")
library(parallel)
library(BiodiversityR) # ALWAYS LOAD FROM THE R CONSOLE!!!!
library(ggforce)

### 1. AOA

### 1 A. Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
aoa.bulk_dist_bc
# metadata
str(aoa.meta.bulk.ed)
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
aoa.meta.bulk.ed$Block <- as.factor(aoa.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ x, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 44 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk <- CAPdiscrim(aoa.bulk_dist_bc ~ x, data = aoa.meta.bulk.ed, m = 44, permutations = 9999, add = TRUE) # 94.16667% ,Significance of this percentage was 0.00010001 
aoa.cap.bulk$x

#aoa.cap.bulk.dist <- dist(aoa.cap.bulk$PCoA)
#x: the positions of the sites provided by the discriminant analysis

dist_matrix.aoa <- dist(aoa.cap.bulk$x)
dist_matrix.aoa


success <- cbind(data.frame(aoa.cap.bulk$group), data.frame(aoa.cap.bulk$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.bulk <- paste("CAP1 (", round((100/sum(aoa.cap.bulk$lda.other$svd^2) * aoa.cap.bulk$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2.bulk <- paste("CAP2 (", round((100/sum(aoa.cap.bulk$lda.other$svd^2) * aoa.cap.bulk$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.bulk$x[, 1:2], xlab = aoa.cap1.bulk, ylab = aoa.cap2.bulk, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.bulk.ed$x])
ordiellipse(aoa.cap.bulk$x[, 1:2], groups = aoa.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (90%)", "Drought-induced BIODYN (95%)",
                                 "Control CONFYM (100%)", "Drought-induced CONFYM (100%)",
                                 "Control CONMIN (95%)", "Drought-induced CONMIN (85%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.plot <- ggplot(as.data.frame(aoa.cap.bulk$x), aes(x = aoa.cap.bulk$x[,1], y = aoa.cap.bulk$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk.ed$Treatment, shape = aoa.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.bulk) + ylab(aoa.cap2.bulk) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  #labs(subtitle = "C. AOA")+
  labs(subtitle = "B. AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20)) +
  annotate("text",x=-9.5,y=-1,label= "90%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-7.5,y=3,label= "95%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-1.7,y=-7,label= "100%", hjust = 0, size = 4, color="#FF618C")+
  annotate("text",x=-1.7,y=-1.5,label= "100%", hjust = 0, size = 4,color="#FF618C")+
  annotate("text",x=2,y=6,label= "95%", hjust = 0, size = 4, color="#E69F00")+
  annotate("text",x=6.5,y=0.4,label= "85%", hjust = 0, size = 4,color="#E69F00")+
  annotate("text",x=-14,y=-10,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) +
  annotate("text", x=-14, y=-11.5, label= "Pillai's test=4.1***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
aoa.cap.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_CAP_bulk_bray.tiff",
       aoa.cap.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)

### 1 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aoa.rh_dist_bc <- vegdist(t(aoa.asv.rh1), method = "bray")
aoa.rh_dist_bc
# metadata
aoa.meta.rh
#aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
#aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.rh_dist_bc ~ x, data = aoa.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 38 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.rh <- CAPdiscrim(aoa.rh_dist_bc ~ x, data = aoa.meta.rh.ed, m = 38, permutations = 9999, add = TRUE) # 90.27778 % ,Significance of this percentage was 0.0001
aoa.cap.rh

dist_matrix.aoa.rh <- dist(aoa.cap.rh$x)

success <- cbind(data.frame(aoa.cap.rh$group), data.frame(aoa.cap.rh$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.rh$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.rh <- paste("CAP1 (", round((100/sum(aoa.cap.rh$lda.other$svd^2) * aoa.cap.rh$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2.rh <- paste("CAP2 (", round((100/sum(aoa.cap.rh$lda.other$svd^2) * aoa.cap.rh$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.rh$x[, 1:2], xlab = aoa.cap1.rh, ylab = aoa.cap2.rh, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.rh.ed$x])
ordiellipse(aoa.cap.rh$x[, 1:2], groups = aoa.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (100%)", "Drought-induced BIODYN (83.33%)",
                                "Control CONFYM (91.7%)", "Drought-induced CONFYM (75%)",
                                "Control CONMIN (100%)", "Drought-induced CONMIN (91.7%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.rh.plot <- ggplot(as.data.frame(aoa.cap.rh$x), aes(x = aoa.cap.rh$x[,1], y = aoa.cap.rh$x[,2])) +
  geom_point(aes(color = aoa.meta.rh.ed$Treatment, shape = aoa.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.rh.ed$PlotID)+
  xlab(aoa.cap1.rh) + ylab(aoa.cap2.rh) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Drought",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
 labs(subtitle = "D. AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20)) +
 annotate("text",x=-14,y=-1,label= "100%", hjust = 0, size = 4,color="#009E73") +
 annotate("text",x=-14,y=3,label= "83.3%", hjust = 0, size = 4,color="#009E73") +
 annotate("text",x=-3,y=-7,label= "91.7%", hjust = 0, size = 4, color="#FF618C")+
 annotate("text",x=1,y=-2,label= "75%", hjust = 0, size = 4,color="#FF618C")+
 annotate("text",x=3,y=6,label= "100%", hjust = 0, size = 4, color="#E69F00")+
 annotate("text",x=6,y=1,label= "91.7%", hjust = 0, size = 4,color="#E69F00")+
annotate("text",x=-17,y=-9,label= "Overall reclassification rate: 90.3%", hjust = 0, size = 4) +
annotate("text", x=-17, y=-10.5, label= "Pillai's test=4.4***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
aoa.cap.rh.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_CAP_rhizo_bray.tiff",
       aoa.cap.rh.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)
##############################################################################################################################################

### 2. COMAMMOX

# 2 A. Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
com.bulk_dist_bc <- vegdist(t(com.asv.bulk1), method = "bray")
com.bulk_dist_bc
# metadata
com.meta.bulk
com.meta.bulk.ed <- com.meta.bulk[,-30:-45]
com.meta.bulk.ed$x <- as.factor(com.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(com.bulk_dist_bc ~ x, data = com.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 49 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.bulk <- CAPdiscrim(com.bulk_dist_bc ~ x, data = com.meta.bulk.ed, m = 49, permutations = 9999, add = TRUE) # 78.81356% , Significance of this percentage was 0.0001
com.cap.bulk


dist_matrix.com <- dist(com.cap.bulk$x)


success <- cbind(data.frame(com.cap.bulk$group), data.frame(com.cap.bulk$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.bulk$PCoA)
success <- success[order(success$source), ]
success

com.cap1.bulk <- paste("CAP1 (", round((100/sum(com.cap.bulk$lda.other$svd^2) * com.cap.bulk$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
com.cap2.bulk <- paste("CAP2 (", round((100/sum(com.cap.bulk$lda.other$svd^2) * com.cap.bulk$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.bulk$x[, 1:2], xlab = com.cap1.bulk, ylab = com.cap2.bulk, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.bulk.ed$x])
ordiellipse(com.cap.bulk$x[, 1:2], groups = com.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (90%)", "Drought-induced BIODYN (89.5%)",
                                "Control CONFYM (85%)", "Drought-induced CONFYM (65%)",
                                "Control CONMIN (80%)", "Drought-induced CONMIN (63.2%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.plot <- ggplot(as.data.frame(com.cap.bulk$x), aes(x = com.cap.bulk$x[,1], y = com.cap.bulk$x[,2])) +
  geom_point(aes(color = com.meta.bulk.ed$Treatment, shape = com.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(com.cap1.bulk) + ylab(com.cap2.bulk) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Drought",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
 #labs(subtitle = "E. Comammox")+
  labs(subtitle = "C. Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20)) +
  annotate("text",x=-9,y=2.5,label= "90%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-15,y=3,label= "89.5%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=1,y=-2,label= "85%", hjust = 0, size = 4, color="#FF618C")+
  annotate("text",x=6,y=-7,label= "65%", hjust = 0, size = 4,color="#FF618C")+
  annotate("text",x=6,y=5,label= "80%", hjust = 0, size = 4, color="#E69F00")+
  annotate("text",x=-2.5,y=5,label= "63.2%", hjust = 0, size = 4,color="#E69F00")+
  annotate("text",x=-17,y=-7,label= "Overall reclassification rate: 78.8%", hjust = 0, size = 4) +
  annotate("text", x=-17, y=-8.5, label= "Pillai's test=3.7***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
com.cap.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COM_CAP_bulk_bray.tiff",
       com.cap.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)

### 2 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
com.rh_dist_bc <- vegdist(t(com.asv.rh1), method = "bray")
com.rh_dist_bc
# metadata
str(com.meta.rh)
#com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
#com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.rh_dist_bc ~ x, data = com.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 19 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.rh <- CAPdiscrim(com.rh_dist_bc ~ x, data = com.meta.rh, m = 19, permutations = 9999, add = TRUE) #  83.33%, Significance of this percentage was 0.0001
com.cap.rh

dist_matrix.com.rh <- dist(com.cap.rh$x)

success <- cbind(data.frame(com.cap.rh$group), data.frame(com.cap.rh$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.rh$PCoA)
success <- success[order(success$source), ]
success

com.cap1.rh <- paste("CAP1 (", round((100/sum(com.cap.rh$lda.other$svd^2) * com.cap.rh$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
com.cap2.rh <- paste("CAP2 (", round((100/sum(com.cap.rh$lda.other$svd^2) * com.cap.rh$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.rh$x[, 1:2], xlab = com.cap1.rh, ylab = com.cap2.rh, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.rh$x])
ordiellipse(com.cap.rh$x[, 1:2], groups = com.meta.rh$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (100%)", "Drought-induced BIODYN (91.7%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (83.33)",
                                "Control CONMIN (75%)", "Drought-induced CONMIN (75%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.rh.plot <- ggplot(as.data.frame(com.cap.rh$x), aes(x = com.cap.rh$x[,1], y = com.cap.rh$x[,2])) +
  geom_point(aes(color = com.meta.rh$Treatment, shape = com.meta.rh$Irrigation), size = 2) +
  #geom_text(label=com.meta.rh$PlotID)+
  xlab(com.cap1.rh) + ylab(com.cap2.rh) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Drought",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.rh$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(subtitle = "F. Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20))+
  annotate("text",x=-16,y=-2.5,label= "100%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-16,y=3,label= "91.7%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=0.5,y=-2,label= "75%", hjust = 0, size = 4, color="#FF618C")+
  annotate("text",x=4.4,y=-8.3,label= "83.3%", hjust = 0, size = 4,color="#FF618C")+
  annotate("text",x=0,y=7,label= "75%", hjust = 0, size = 4, color="#E69F00")+
  annotate("text",x=-2,y=3,label= "75%", hjust = 0, size = 4,color="#E69F00")+
  annotate("text",x=-17,y=-7,label= "Overall reclassification rate: 83.3%", hjust = 0, size = 4) +
  annotate("text", x=-17, y=-8.5, label= "Pillai's test=3.4***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
com.cap.rh.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COM_CAP_rhizo_legend.tiff",
       com.cap.rh.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)
#################################################################################################################################################

### 3. AOB

# 3 A. Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aob.bulk_dist_bc <- vegdist(t(aob.asv.bulk1), method = "bray")
aob.bulk_dist_bc
# metadata
aob.meta.bulk
aob.meta.bulk.ed <- aob.meta.bulk[,-30:-45]
aob.meta.bulk.ed$x <- as.factor(aob.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.bulk_dist_bc ~ x, data = aob.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 35 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate

set.seed(333)
aob.cap.bulk <- CAPdiscrim(aob.bulk_dist_bc ~ x, data = aob.meta.bulk, m = 35, permutations = 9999, add = TRUE) # 60.5042 % ,Significance of this percentage was 0.0001
aob.cap.bulk

dist_matrix.aob <- dist(aob.cap.bulk$x)

success <- cbind(data.frame(aob.cap.bulk$group), data.frame(aob.cap.bulk$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.bulk$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.bulk <- paste("CAP1 (", round((100/sum(aob.cap.bulk$lda.other$svd^2) * aob.cap.bulk$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aob.cap2.bulk <- paste("CAP2 (", round((100/sum(aob.cap.bulk$lda.other$svd^2) * aob.cap.bulk$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.bulk$x[, 1:2], xlab = aob.cap1.bulk, ylab = aob.cap2.bulk, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.bulk.ed$x])
ordiellipse(aob.cap.bulk$x[, 1:2], groups = aob.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (57.9%)", "Drought-induced BIODYN (80%)",
                                "Control CONFYM (45%)", "Drought-induced CONFYM (55%)",
                                "Control CONMIN (70%)", "Drought-induced CONMIN (55%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.plot <- ggplot(as.data.frame(aob.cap.bulk$x), aes(x = aob.cap.bulk$x[,1], y = aob.cap.bulk$x[,2])) +
  geom_point(aes(color = aob.meta.bulk$Treatment, shape = aob.meta.bulk$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk$PlotID)+
  xlab(aob.cap1.bulk) + ylab(aob.cap2.bulk) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.bulk$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  #labs(title = "Bulk Soil", subtitle = "A. AOB")+
  labs(subtitle = "A. AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold"),
        plot.title = element_text(size=25, face="bold")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20)) +
  annotate("text",x=-3,y=1.6,label= "57.9%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-5.5,y=3.5,label= "80%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-2.2,y=-2.2,label= "45%", hjust = 0, size = 4, color="#FF618C")+
  annotate("text",x=3.5,y=-2.6,label= "55%", hjust = 0, size = 4,color="#FF618C")+
  annotate("text",x=3.6,y=5,label= "70%", hjust = 0, size = 4, color="#E69F00")+
  annotate("text",x=-1,y=3,label= "55%", hjust = 0, size = 4,color="#E69F00")+
  annotate("text",x=-6,y=-5.5,label= "Overall reclassification rate: 60.5%", hjust = 0, size = 4) +
  annotate("text", x=-6, y=-6.5, label= "Pillai's test=2.7***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
aob.cap.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB_CAP_bulk_bray.tiff",
       aob.cap.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)


### 3 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aob.rh_dist_bc <- vegdist(t(aob.asv.rh1), method = "bray")
aob.rh_dist_bc
# metadata
aob.meta.rh
#aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
#aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_bc ~ x, data = aob.meta.rh, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 20 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh <- CAPdiscrim(aob.rh_dist_bc ~ x, data = aob.meta.rh, m = 20, permutations = 9999, add = TRUE) # 54.16667 % ,Significance of this percentage was 0.0001
aob.cap.rh

dist_matrix.aob.rh <- dist(aob.cap.rh$x)


success <- cbind(data.frame(aob.cap.rh$group), data.frame(aob.cap.rh$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.rh$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.rh <- paste("CAP1 (", round((100/sum(aob.cap.rh$lda.other$svd^2) * aob.cap.rh$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
aob.cap2.rh <- paste("CAP2 (", round((100/sum(aob.cap.rh$lda.other$svd^2) * aob.cap.rh$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.rh$x[, 1:2], xlab = aob.cap1.rh, ylab = aob.cap2.rh, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.rh.ed$x])
ordiellipse(aob.cap.rh$x[, 1:2], groups = aob.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (66.7%)", "Drought-induced BIODYN (58.3%)",
                                "Control CONFYM (66.7%)", "Drought-induced CONFYM (50%)",
                                "Control CONMIN (41.7%)", "Drought-induced CONMIN (41.7%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.rh.plot <- ggplot(as.data.frame(aob.cap.rh$x), aes(x = aob.cap.rh$x[,1], y = aob.cap.rh$x[,2])) +
  geom_point(aes(color = aob.meta.rh$Treatment, shape = aob.meta.rh$Irrigation), size = 2) +
  #geom_text(label=aob.meta.rh.ed$PlotID)+
  xlab(aob.cap1.rh) + ylab(aob.cap2.rh) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.rh$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "Rhizosphere", subtitle = "B. AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.subtitle = element_text(size=20, face="bold"),
        plot.title = element_text(size=25, face="bold")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.text = element_text(size=20))+
  annotate("text",x=-5,y=2,label= "66.7%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-5.5,y=-3,label= "58.3%", hjust = 0, size = 4,color="#009E73") +
  annotate("text",x=-2.7,y=3,label= "66.7%", hjust = 0, size = 4, color="#FF618C")+
  annotate("text",x=2.5,y=5,label= "50%", hjust = 0, size = 4,color="#FF618C")+
  annotate("text",x=3.3,y=-3.5,label= "41.7%", hjust = 0, size = 4, color="#E69F00")+
  annotate("text",x=-1.3,y=-2,label= "41.7%", hjust = 0, size = 4,color="#E69F00")+
  annotate("text",x=-6,y=-5.5,label= "Overall reclassification rate: 54.2%", hjust = 0, size = 4) +
  annotate("text", x=-6, y=-6.5, label= "Pillai's test=2.6***", hjust = 0, size = 4)+
 guides(colour=guide_legend(override.aes = list(size=7)),shape=guide_legend(override.aes = list(size=7)))
aob.cap.rh.plot
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB_CAP_rhizo_bray.tiff",
       aob.cap.rh.plot, device = "tiff",
       width = 4, height =3, 
       units= "in", dpi = 600)
############################################################################################################################################
# Save all the plots 
library(patchwork)
cap.all.BS <- aob.cap.plot / aoa.cap.plot / com.cap.plot
cap.all.BS
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.BS.tiff",
       cap.all.BS, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

cap.all.RZ <- aob.cap.rh.plot / aoa.cap.rh.plot / com.cap.rh.plot
cap.all.RZ
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.RZ.tiff",
       cap.all.RZ, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

library(ggpubr)
cap.all.RZ <- ggarrange(aob.cap.rh.plot,aoa.cap.rh.plot,com.cap.rh.plot,
                        nrow = 1,
                        common.legend = T, 
                        legend = "bottom") 
                        #ncol=2, 
                        #nrow = 2,
                        #common.legend = T,
                        #legend = "right")
cap.all.RZ
#setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.RZ.legend.tiff",
       cap.all.RZ, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)



CAP.All <- ((aob.cap.plot / aoa.cap.plot / com.cap.plot) | (aob.cap.rh.plot / aoa.cap.rh.plot / com.cap.rh.plot))+
 plot_layout(guides = "collect") & theme(legend.position = 'right',legend.title = element_blank(), legend.text = element_text(size = 22))
CAP.All
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.2dpi300.tiff",
       CAP.All, device = "tiff",
       width = 10, height = 12, 
       units= "in", dpi = 300, compression="lzw")

ggsave("Fig.2.2dpi300.tiff",
       CAP.All, device = "tiff",
       width = 10, height = 15, 
       units= "in", dpi = 300, compression="lzw")

CAP.BS <- (aob.cap.plot / aoa.cap.plot / com.cap.plot)
CAP.BS
CAP.RS <- (aob.cap.rh.plot / aoa.cap.rh.plot / com.cap.rh.plot)
ggsave("Fig.2.CAP.BS.tiff",
       CAP.BS, device = "tiff",
       width = 5, height = 11, 
       units= "in", dpi = 300, compression="lzw")
ggsave("Fig.2.CAP.RS.tiff",
       CAP.RS, device = "tiff",
       width = 5, height = 11, 
       units= "in", dpi = 300, compression="lzw")

CAP.BS.horiz <- (aob.cap.plot | aoa.cap.plot| com.cap.plot) +
 plot_layout(guides = "collect") & theme(legend.position = 'bottom',legend.title = element_blank(), legend.text = element_text(size = 22))
CAP.BS.horiz
ggsave("Fig.CAP.BS.horiz.tiff",
       CAP.BS.horiz, device = "tiff",
       width = 15, height = 6, 
       units= "in", dpi = 300, compression="lzw")



##################################################################################################################

# Using Weighted UniFrac Distance

### 1. AOA

### 1 A. Bulk Soil

# run Weighted UniFrac beta diversity on bulk soil
aoa.bulk_dist_wUF 
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aoa.bulk_dist_wUF ~ x, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 31 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.wUF <- CAPdiscrim(aoa.bulk_dist_wUF ~ x, data = aoa.meta.bulk.ed, m = 31, permutations = 9999, add = TRUE) # 85% 

success <- cbind(data.frame(aoa.cap.bulk.wUF$group), data.frame(aoa.cap.bulk.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.wUF$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.bulk.wUF <- paste("CAP1 (", round((100/sum(aoa.cap.bulk.wUF$lda.other$svd^2) * aoa.cap.bulk.wUF$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2.bulk.wuF <- paste("CAP2 (", round((100/sum(aoa.cap.bulk.wUF$lda.other$svd^2) * aoa.cap.bulk.wUF$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.bulk.wUF$x[, 1:2], xlab = aoa.cap1.bulk.wUF, ylab = aoa.cap2.bulk.wuF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.bulk.ed$x])
ordiellipse(aoa.cap.bulk.wUF$x[, 1:2], groups = aoa.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (90%)", "Drought-induced BIODYN (95%)",
                                "Control CONFYM (85%)", "Drought-induced CONFYM (80%)",
                                "Control CONMIN (75%)", "Drought-induced CONMIN (85%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.wUF.plot <- ggplot(as.data.frame(aoa.cap.bulk.wUF$x), aes(x = aoa.cap.bulk.wUF$x[,1], y = aoa.cap.bulk.wUF$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk.ed$Treatment, shape = aoa.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.bulk.wUF) + ylab(aoa.cap2.bulk.wuF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aoa.cap.wUF.plot

### 1 B. AOA Rhizosphere

# run Weighted UniFrac beta diversity on rhizosphere
aoa.rh_dist_wUF
# metadata
aoa.meta.rh
aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.rh_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.rh_dist_wUF ~ x, data = aoa.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 35 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.rh.wUF <- CAPdiscrim(aoa.rh_dist_wUF ~ x, data = aoa.meta.rh.ed, m = 35, permutations = 9999, add = TRUE) #  % 

success <- cbind(data.frame(aoa.cap.rh.wUF$group), data.frame(aoa.cap.rh.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.rh.wUF$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.rh.wUF <- paste("CAP1 (", round((100/sum(aoa.cap.rh.wUF$lda.other$svd^2) * aoa.cap.rh.wUF$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
aoa.cap2.rh.wUF <- paste("CAP2 (", round((100/sum(aoa.cap.rh.wUF$lda.other$svd^2) * aoa.cap.rh.wUF$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.rh.wUF$x[, 1:2], xlab = aoa.cap1.rh.wUF, ylab = aoa.cap2.rh.wUF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.rh.ed$x])
ordiellipse(aoa.cap.rh.wUF$x[, 1:2], groups = aoa.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (91.7%)", "Drought-induced BIODYN (91.7%)",
                                "Control CONFYM (91.7%)", "Drought-induced CONFYM (83.33%)",
                                "Control CONMIN (91.7%)", "Drought-induced CONMIN (100%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.rh.wUF.plot <- ggplot(as.data.frame(aoa.cap.rh.wUF$x), aes(x = aoa.cap.rh.wUF$x[,1], y = aoa.cap.rh.wUF$x[,2])) +
  geom_point(aes(color = aoa.meta.rh.ed$Treatment, shape = aoa.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.rh.ed$PlotID)+
  xlab(aoa.cap1.rh.wUF) + ylab(aoa.cap2.rh.wUF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aoa.cap.rh.wUF.plot

#################################################################################################################################################

### 2. AOB

# 3 A. Bulk Soil

# run Weighted UniFrac beta diversity on bulk soil
aob.bulk_dist_wUF
# metadata
aob.meta.bulk
aob.meta.bulk.ed <- aob.meta.bulk[,-30:-45]
aob.meta.bulk.ed$x <- as.factor(aob.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.bulk_dist_wUF ~ x, data = aob.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 24 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate

set.seed(333)
aob.cap.bulk.wUF <- CAPdiscrim(aob.bulk_dist_wUF ~ x, data = aob.meta.bulk.ed, m = 24, permutations = 9999, add = TRUE) # 57.98 % 

success <- cbind(data.frame(aob.cap.bulk.wUF$group), data.frame(aob.cap.bulk.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.bulk.wUF$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.bulk.wUF <- paste("CAP1 (", round((100/sum(aob.cap.bulk.wUF$lda.other$svd^2) * aob.cap.bulk.wUF$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aob.cap2.bulk.wUF <- paste("CAP2 (", round((100/sum(aob.cap.bulk.wUF$lda.other$svd^2) * aob.cap.bulk.wUF$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.bulk.wUF$x[, 1:2], xlab = aob.cap1.bulk.wUF, ylab = aob.cap2.bulk.wUF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.bulk.ed$x])
ordiellipse(aob.cap.bulk.wUF$x[, 1:2], groups = aob.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (57.9%)", "Drought-induced BIODYN (50%)",
                                "Control CONFYM (65%)", "Drought-induced CONFYM (65%)",
                                "Control CONMIN (55%)", "Drought-induced CONMIN (55%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.wUF.plot <- ggplot(as.data.frame(aob.cap.bulk.wUF$x), aes(x = aob.cap.bulk.wUF$x[,1], y = aob.cap.bulk.wUF$x[,2])) +
  geom_point(aes(color = aob.meta.bulk.ed$Treatment, shape = aob.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aob.cap1.bulk.wUF) + ylab(aob.cap2.bulk.wUF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 60.5%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aob.cap.wUF.plot

### 3 B. Rhizosphere

# run Weighted UniFrac beta diversity on rhizosphere
aob.rh_dist_wUF
# metadata
aob.meta.rh
aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_wUF ~ x, data = aob.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 31 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh.wUF <- CAPdiscrim(aob.rh_dist_wUF ~ x, data = aob.meta.rh.ed, m = 31, permutations = 9999, add = TRUE) # 69.44 % 

success <- cbind(data.frame(aob.cap.rh.wUF$group), data.frame(aob.cap.rh.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.rh.wUF$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.rh.wUF <- paste("CAP1 (", round((100/sum(aob.cap.rh.wUF$lda.other$svd^2) * aob.cap.rh.wUF$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
aob.cap2.rh.wUF <- paste("CAP2 (", round((100/sum(aob.cap.rh.wUF$lda.other$svd^2) * aob.cap.rh.wUF$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.rh.wUF$x[, 1:2], xlab = aob.cap1.rh.wUF, ylab = aob.cap2.rh.wUF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.rh.ed$x])
ordiellipse(aob.cap.rh.wUF$x[, 1:2], groups = aob.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (91.7%)", "Drought-induced BIODYN (83.33%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (66.7%)",
                                "Control CONMIN (33.33%)", "Drought-induced CONMIN (66.7%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.rh.wUF.plot <- ggplot(as.data.frame(aob.cap.rh.wUF$x), aes(x = aob.cap.rh.wUF$x[,1], y = aob.cap.rh.wUF$x[,2])) +
  geom_point(aes(color = aob.meta.rh.ed$Treatment, shape = aob.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aob.meta.rh.ed$PlotID)+
  xlab(aob.cap1.rh.wUF) + ylab(aob.cap2.rh.wUF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))

#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aob.cap.rh.wUF.plot

############################################################################################################################
### 3. COMAMMOX

# 3 A. Bulk Soil

# run Weighted UniFrac beta diversity on bulk soil
com.bulk_dist_wUF
# metadata
com.meta.bulk
com.meta.bulk.ed <- com.meta.bulk[,-30:-45]
com.meta.bulk.ed$x <- as.factor(com.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(com.bulk_dist_wUF ~ x, data = com.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 37 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.bulk.wUF <- CAPdiscrim(com.bulk_dist_wUF ~ x, data = com.meta.bulk.ed, m = 37, permutations = 9999, add = TRUE) # 62.711% 

success <- cbind(data.frame(com.cap.bulk.wUF$group), data.frame(com.cap.bulk.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.bulk.wUF$PCoA)
success <- success[order(success$source), ]
success

com.cap1.bulk.wUF <- paste("CAP1 (", round((100/sum(com.cap.bulk.wUF$lda.other$svd^2) * com.cap.bulk.wUF$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
com.cap2.bulk.wUF <- paste("CAP2 (", round((100/sum(com.cap.bulk.wUF$lda.other$svd^2) * com.cap.bulk.wUF$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.bulk.wUF$x[, 1:2], xlab = com.cap1.bulk.wUF, ylab = com.cap2.bulk.wUF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.bulk.ed$x])
ordiellipse(com.cap.bulk.wUF$x[, 1:2], groups = com.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (65%)", "Drought-induced BIODYN (47.3%)",
                                "Control CONFYM (65%)", "Drought-induced CONFYM (45%)",
                                "Control CONMIN (80%)", "Drought-induced CONMIN (73.7%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.wUF.plot <- ggplot(as.data.frame(com.cap.bulk.wUF$x), aes(x = com.cap.bulk.wUF$x[,1], y = com.cap.bulk.wUF$x[,2])) +
  geom_point(aes(color = com.meta.bulk.ed$Treatment, shape = com.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(com.cap1.bulk.wUF) + ylab(com.cap2.bulk.wUF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 79.7%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
com.cap.wUF.plot

### 3 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
com.rh_dist_wUF
# metadata
com.meta.rh
com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.rh_dist_wUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.rh_dist_wUF ~ x, data = com.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 33 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.rh.wUF <- CAPdiscrim(com.rh_dist_wUF ~ x, data = com.meta.rh.ed, m = 33, permutations = 9999, add = TRUE) #  81.94% 

success <- cbind(data.frame(com.cap.rh.wUF$group), data.frame(com.cap.rh.wUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.rh.wUF$PCoA)
success <- success[order(success$source), ]
success

com.cap1.rh.wUF <- paste("CAP1 (", round((100/sum(com.cap.rh.wUF$lda.other$svd^2) * com.cap.rh.wUF$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
com.cap2.rh.wUF <- paste("CAP2 (", round((100/sum(com.cap.rh.wUF$lda.other$svd^2) * com.cap.rh.wUF$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.rh.wUF$x[, 1:2], xlab = com.cap1.rh.wUF, ylab = com.cap2.rh.wUF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.rh.ed$x])
ordiellipse(com.cap.rh.wUF$x[, 1:2], groups = com.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (83.33%)", "Drought-induced BIODYN (91.7%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (83.33)",
                                "Control CONMIN (83.33%)", "Drought-induced CONMIN (75%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.rh.wUF.plot <- ggplot(as.data.frame(com.cap.rh.wUF$x), aes(x = com.cap.rh.wUF$x[,1], y = com.cap.rh.wUF$x[,2])) +
  geom_point(aes(color = com.meta.rh.ed$Treatment, shape = com.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=com.meta.rh.ed$PlotID)+
  xlab(com.cap1.rh.wUF) + ylab(com.cap2.rh.wUF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
#guides(colour=guide_legend(override.aes = list(size=4)))
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
com.cap.rh.wUF.plot

# Save all the plots 
library(patchwork)
cap.all.wUF.BS <- aob.cap.wUF.plot | aoa.cap.wUF.plot | com.cap.wUF.plot
cap.all.wUF.BS
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.wUF.BS.tiff",
       cap.all.wUF.BS, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

cap.all.wUF.RZ <- aob.cap.rh.wUF.plot | aoa.cap.rh.wUF.plot | com.cap.rh.wUF.plot
cap.all.wUF.RZ
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.wUF.RZ.tiff",
       cap.all.wUF.RZ, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

#########################################################################################################################

# run unWeighted UniFrac beta diversity on bulk soil
aoa.bulk_dist_uwUF 
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_uwUF))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aoa.bulk_dist_uwUF ~ x, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 62 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.uwUF <- CAPdiscrim(aoa.bulk_dist_uwUF ~ x, data = aoa.meta.bulk.ed, m = 62, permutations = 9999, add = TRUE) # 61.67% 

success <- cbind(data.frame(aoa.cap.bulk.uwUF$group), data.frame(aoa.cap.bulk.uwUF$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.uwUF$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.bulk.uwUF <- paste("CAP1 (", round((100/sum(aoa.cap.bulk.uwUF$lda.other$svd^2) * aoa.cap.bulk.uwUF$lda.other$svd^2)[1],
                                           digits = 1), "%)", sep = "")
aoa.cap2.bulk.uwuF <- paste("CAP2 (", round((100/sum(aoa.cap.bulk.uwUF$lda.other$svd^2) * aoa.cap.bulk.uwUF$lda.other$svd^2)[2],
                                           digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.bulk.uwUF$x[, 1:2], xlab = aoa.cap1.bulk.uwUF, ylab = aoa.cap2.bulk.uwuF, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.bulk.ed$x])
ordiellipse(aoa.cap.bulk.uwUF$x[, 1:2], groups = aoa.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (65%)", "Drought-induced BIODYN (50%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (55%)",
                                "Control CONMIN (60%)", "Drought-induced CONMIN (65%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.uwUF.plot <- ggplot(as.data.frame(aoa.cap.bulk.uwUF$x), aes(x = aoa.cap.bulk.uwUF$x[,1], y = aoa.cap.bulk.uwUF$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk.ed$Treatment, shape = aoa.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.bulk.uwUF) + ylab(aoa.cap2.bulk.uwuF) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aoa.cap.uwUF.plot

##########################################################################################################################

# Using Jaccard Distance

### 1. AOA

### 1 A. Bulk Soil

# run Jaccard beta diversity on bulk soil
aoa.bulk_dist_jac 
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aoa.bulk_dist_jac ~ x, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 50 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.jac <- CAPdiscrim(aoa.bulk_dist_jac ~ x, data = aoa.meta.bulk.ed, m = 50, permutations = 9999, add = TRUE) # 77.5% 

success <- cbind(data.frame(aoa.cap.bulk.jac$group), data.frame(aoa.cap.bulk.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.jac$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.bulk.jac <- paste("CAP1 (", round((100/sum(aoa.cap.bulk.jac$lda.other$svd^2) * aoa.cap.bulk.jac$lda.other$svd^2)[1],
                                           digits = 1), "%)", sep = "")
aoa.cap2.bulk.jac <- paste("CAP2 (", round((100/sum(aoa.cap.bulk.jac$lda.other$svd^2) * aoa.cap.bulk.jac$lda.other$svd^2)[2],
                                           digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.bulk.jac$x[, 1:2], xlab = aoa.cap1.bulk.jac, ylab = aoa.cap2.bulk.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.bulk.ed$x])
ordiellipse(aoa.cap.bulk.jac$x[, 1:2], groups = aoa.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (75%)", "Drought-induced BIODYN (70%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (75%)",
                                "Control CONMIN (90%)", "Drought-induced CONMIN (80%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.jac.BS.plot <- ggplot(as.data.frame(aoa.cap.bulk.jac$x), aes(x = aoa.cap.bulk.jac$x[,1], y = aoa.cap.bulk.jac$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk.ed$Treatment, shape = aoa.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.bulk.jac) + ylab(aoa.cap2.bulk.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aoa.cap.jac.BS.plot

### 1 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aoa.rh_dist_jac
# metadata
aoa.meta.rh
aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.rh_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.rh_dist_jac ~ x, data = aoa.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 33 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.rh.jac <- CAPdiscrim(aoa.rh_dist_jac ~ x, data = aoa.meta.rh.ed, m = 33, permutations = 9999, add = TRUE) # 54.17% 

success <- cbind(data.frame(aoa.cap.rh.jac$group), data.frame(aoa.cap.rh.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.rh.jac$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.rh.jac <- paste("CAP1 (", round((100/sum(aoa.cap.rh.jac$lda.other$svd^2) * aoa.cap.rh.jac$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
aoa.cap2.rh.jac <- paste("CAP2 (", round((100/sum(aoa.cap.rh.jac$lda.other$svd^2) * aoa.cap.rh.jac$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(aoa.cap.rh.jac$x[, 1:2], xlab = aoa.cap1.rh.jac, ylab = aoa.cap2.rh.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.meta.rh.ed$x])
ordiellipse(aoa.cap.rh.jac$x[, 1:2], groups = aoa.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (33.33%)", "Drought-induced BIODYN (33.33%)",
                                "Control CONFYM (58.33%)", "Drought-induced CONFYM (66.67%)",
                                "Control CONMIN (58.33%)", "Drought-induced CONMIN (75%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.rh.jac.plot <- ggplot(as.data.frame(aoa.cap.rh.jac$x), aes(x = aoa.cap.rh.jac$x[,1], y = aoa.cap.rh.jac$x[,2])) +
  geom_point(aes(color = aoa.meta.rh.ed$Treatment, shape = aoa.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.rh.ed$PlotID)+
  xlab(aoa.cap1.rh.jac) + ylab(aoa.cap2.rh.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aoa.cap.rh.jac.plot

### 1. AOB

### 1 A. Bulk Soil

# run Jaccard beta diversity on bulk soil
aob.bulk_dist_jac 
# metadata
aob.meta.bulk
aob.meta.bulk.ed <- aob.meta.bulk[,-30:-45]
aob.meta.bulk.ed$x <- as.factor(aob.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aob.bulk_dist_jac ~ x, data = aob.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 76 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.bulk.jac <- CAPdiscrim(aob.bulk_dist_jac ~ x, data = aob.meta.bulk.ed, m = 76, permutations = 9999, add = TRUE) # 52.1% 

success <- cbind(data.frame(aob.cap.bulk.jac$group), data.frame(aob.cap.bulk.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.bulk.jac$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.bulk.jac <- paste("CAP1 (", round((100/sum(aob.cap.bulk.jac$lda.other$svd^2) * aob.cap.bulk.jac$lda.other$svd^2)[1],
                                           digits = 1), "%)", sep = "")
aob.cap2.bulk.jac <- paste("CAP2 (", round((100/sum(aob.cap.bulk.jac$lda.other$svd^2) * aob.cap.bulk.jac$lda.other$svd^2)[2],
                                           digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.bulk.jac$x[, 1:2], xlab = aob.cap1.bulk.jac, ylab = aob.cap2.bulk.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.bulk.ed$x])
ordiellipse(aob.cap.bulk.jac$x[, 1:2], groups = aob.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (57.9%)", "Drought-induced BIODYN (60%)",
                                "Control CONFYM (55%)", "Drought-induced CONFYM (45%)",
                                "Control CONMIN (55%)", "Drought-induced CONMIN (40%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.jac.BS.plot <- ggplot(as.data.frame(aob.cap.bulk.jac$x), aes(x = aob.cap.bulk.jac$x[,1], y = aob.cap.bulk.jac$x[,2])) +
  geom_point(aes(color = aob.meta.bulk.ed$Treatment, shape = aob.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aob.cap1.bulk.jac) + ylab(aob.cap2.bulk.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aob.cap.jac.BS.plot

### 2 B. AOB Rhizosphere

# run Jaccard beta diversity on rhizosphere
aob.rh_dist_jac
# metadata
aob.meta.rh
aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_jac ~ x, data = aob.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 10 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh.jac <- CAPdiscrim(aob.rh_dist_jac ~ x, data = aob.meta.rh.ed, m = 10, permutations = 9999, add = TRUE) #  47.22% 

success <- cbind(data.frame(aob.cap.rh.jac$group), data.frame(aob.cap.rh.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.rh.jac$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.rh.jac <- paste("CAP1 (", round((100/sum(aob.cap.rh.jac$lda.other$svd^2) * aob.cap.rh.jac$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
aob.cap2.rh.jac <- paste("CAP2 (", round((100/sum(aob.cap.rh.jac$lda.other$svd^2) * aob.cap.rh.jac$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(aob.cap.rh.jac$x[, 1:2], xlab = aob.cap1.rh.jac, ylab = aob.cap2.rh.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aob.meta.rh.ed$x])
ordiellipse(aob.cap.rh.jac$x[, 1:2], groups = aob.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (50%)", "Drought-induced BIODYN (25%)",
                                "Control CONFYM (50%)", "Drought-induced CONFYM (50%)",
                                "Control CONMIN (50%)", "Drought-induced CONMIN (58.3%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aob.cap.rh.jac.plot <- ggplot(as.data.frame(aob.cap.rh.jac$x), aes(x = aob.cap.rh.jac$x[,1], y = aob.cap.rh.jac$x[,2])) +
  geom_point(aes(color = aob.meta.rh.ed$Treatment, shape = aob.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aob.meta.rh.ed$PlotID)+
  xlab(aob.cap1.rh.jac) + ylab(aob.cap2.rh.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aob.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))

#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
aob.cap.rh.jac.plot


### 3. COMAMMOX

# 3 A. Bulk Soil

# run Jaccard beta diversity on bulk soil
com.bulk_dist_jac
# metadata
com.meta.bulk
com.meta.bulk.ed <- com.meta.bulk[,-30:-45]
com.meta.bulk.ed$x <- as.factor(com.meta.bulk.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(com.bulk_dist_jac ~ x, data = com.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 35 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.bulk.jac <- CAPdiscrim(com.bulk_dist_jac ~ x, data = com.meta.bulk.ed, m = 35, permutations = 9999, add = TRUE) # 73.73% 

success <- cbind(data.frame(com.cap.bulk.jac$group), data.frame(com.cap.bulk.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.bulk.jac$PCoA)
success <- success[order(success$source), ]
success

com.cap1.bulk.jac <- paste("CAP1 (", round((100/sum(com.cap.bulk.jac$lda.other$svd^2) * com.cap.bulk.jac$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
com.cap2.bulk.jac <- paste("CAP2 (", round((100/sum(com.cap.bulk.jac$lda.other$svd^2) * com.cap.bulk.jac$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.bulk.jac$x[, 1:2], xlab = com.cap1.bulk.jac, ylab = com.cap2.bulk.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.bulk.ed$x])
ordiellipse(com.cap.bulk.jac$x[, 1:2], groups = com.meta.bulk.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (80%)", "Drought-induced BIODYN (78.9%)",
                                "Control CONFYM (55%)", "Drought-induced CONFYM (75%)",
                                "Control CONMIN (75%)", "Drought-induced CONMIN (78.9%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.bulk.jac.plot <- ggplot(as.data.frame(com.cap.bulk.jac$x), aes(x = com.cap.bulk.jac$x[,1], y = com.cap.bulk.jac$x[,2])) +
  geom_point(aes(color = com.meta.bulk.ed$Treatment, shape = com.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(com.cap1.bulk.jac) + ylab(com.cap2.bulk.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.bulk.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) 
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 79.7%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
com.cap.bulk.jac.plot

### 2 B. COMAMMOX Rhizosphere

# run Jaccard beta diversity on rhizosphere
com.rh_dist_jac
# metadata
com.meta.rh
com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.rh_dist_jac))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.rh_dist_jac ~ x, data = com.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 27 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.rh.jac <- CAPdiscrim(com.rh_dist_jac ~ x, data = com.meta.rh.ed, m = 27, permutations = 9999, add = TRUE) # 62.5 % 

success <- cbind(data.frame(com.cap.rh.jac$group), data.frame(com.cap.rh.jac$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.rh.jac$PCoA)
success <- success[order(success$source), ]
success

com.cap1.rh.jac <- paste("CAP1 (", round((100/sum(com.cap.rh.jac$lda.other$svd^2) * com.cap.rh.jac$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
com.cap2.rh.jac <- paste("CAP2 (", round((100/sum(com.cap.rh.jac$lda.other$svd^2) * com.cap.rh.jac$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")

# Plot with basic R

plot(com.cap.rh.jac$x[, 1:2], xlab = com.cap1.rh.jac, ylab = com.cap2.rh.jac, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.rh.ed$x])
ordiellipse(com.cap.rh.jac$x[, 1:2], groups = com.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (50%)", "Drought-induced BIODYN (66.7%)",
                                "Control CONFYM (58.3%)", "Drought-induced CONFYM (50%)",
                                "Control CONMIN (58.3%)", "Drought-induced CONMIN (91.7%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.rh.jac.plot <- ggplot(as.data.frame(com.cap.rh.jac$x), aes(x = com.cap.rh.jac$x[,1], y = com.cap.rh.jac$x[,2])) +
  geom_point(aes(color = com.meta.rh.ed$Treatment, shape = com.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=com.meta.rh.ed$PlotID)+
  xlab(com.cap1.rh.jac) + ylab(com.cap2.rh.jac) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = com.meta.rh.ed$x), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  labs(title = "Comammox")+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
#guides(colour=guide_legend(override.aes = list(size=4)))
#annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) 
#annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)
com.cap.rh.jac.plot

# Save all the plots 
library(patchwork)
cap.all.jac.BS <- aob.cap.jac.BS.plot | aoa.cap.jac.BS.plot | com.cap.bulk.jac.plot
cap.all.jac.BS
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.jac.BS.tiff",
       cap.all.jac.BS, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

cap.all.jac.RZ <- aob.cap.rh.jac.plot | aoa.cap.rh.jac.plot | com.cap.rh.jac.plot
cap.all.jac.RZ
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.jac.RZ.tiff",
       cap.all.jac.RZ, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

############################################################################################
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
data(dataset)
# PCoA
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
# extract the axis scores
tmp <- t1$res_ordination$scores
# differential test with trans_env class
t2 <- trans_env$new(dataset = dataset, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "Group", method = "anova")





















