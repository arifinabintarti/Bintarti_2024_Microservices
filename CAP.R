#install.packages("parallel")
#install.packages("BiodiversityR")
library(parallel)
library(BiodiversityR)
library(ggforce)

### 1. AOA

### 1 A. Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
aoa.bulk_dist_bc
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
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
aoa.cap.bulk <- CAPdiscrim(aoa.bulk_dist_bc ~ x, data = aoa.meta.bulk.ed, m = 44, permutations = 9999, add = TRUE) # 94.16667% 

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
aoa.cap.plot

### 1 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aoa.rh_dist_bc <- vegdist(t(aoa.asv.rh1), method = "bray")
aoa.rh_dist_bc
# metadata
aoa.meta.rh
aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
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
aoa.cap.rh <- CAPdiscrim(aoa.rh_dist_bc ~ x, data = aoa.meta.rh.ed, m = 38, permutations = 9999, add = TRUE) # 90.27778 % 

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
aoa.cap.rh.plot
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
com.cap.bulk <- CAPdiscrim(com.bulk_dist_bc ~ x, data = com.meta.bulk.ed, m = 49, permutations = 9999, add = TRUE) # 78.81356% 

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
com.cap.plot

### 2 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
com.rh_dist_bc <- vegdist(t(com.asv.rh1), method = "bray")
com.rh_dist_bc
# metadata
com.meta.rh
com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
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
com.cap.rh <- CAPdiscrim(com.rh_dist_bc ~ x, data = com.meta.rh.ed, m = 19, permutations = 9999, add = TRUE) #  83.33% 

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
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[com.meta.rh.ed$x])
ordiellipse(com.cap.rh$x[, 1:2], groups = com.meta.rh.ed$x, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (100%)", "Drought-induced BIODYN (91.7%)",
                                "Control CONFYM (75%)", "Drought-induced CONFYM (83.33)",
                                "Control CONMIN (75%)", "Drought-induced CONMIN (75%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

com.cap.rh.plot <- ggplot(as.data.frame(com.cap.rh$x), aes(x = com.cap.rh$x[,1], y = com.cap.rh$x[,2])) +
  geom_point(aes(color = com.meta.rh.ed$Treatment, shape = com.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=com.meta.rh.ed$PlotID)+
  xlab(com.cap1.rh) + ylab(com.cap2.rh) +
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
com.cap.rh.plot
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
aob.cap.bulk <- CAPdiscrim(aob.bulk_dist_bc ~ x, data = aob.meta.bulk.ed, m = 35, permutations = 9999, add = TRUE) # 60.5042 % 

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
  geom_point(aes(color = aob.meta.bulk.ed$Treatment, shape = aob.meta.bulk.ed$Irrigation), size = 2) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aob.cap1.bulk) + ylab(aob.cap2.bulk) +
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
aob.cap.plot

### 3 B. Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aob.rh_dist_bc <- vegdist(t(aob.asv.rh1), method = "bray")
aob.rh_dist_bc
# metadata
aob.meta.rh
aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_bc ~ x, data = aob.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 20 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh <- CAPdiscrim(aob.rh_dist_bc ~ x, data = aob.meta.rh.ed, m = 20, permutations = 9999, add = TRUE) # 54.16667 % 

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
  geom_point(aes(color = aob.meta.rh.ed$Treatment, shape = aob.meta.rh.ed$Irrigation), size = 2) +
  #geom_text(label=aob.meta.rh.ed$PlotID)+
  xlab(aob.cap1.rh) + ylab(aob.cap2.rh) +
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
aob.cap.rh.plot

############################################################################################################################################
# Save all the plots 
library(patchwork)
cap.all.BS <- aob.cap.plot | aoa.cap.plot | com.cap.plot
cap.all.BS
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.BS.tiff",
       cap.all.BS, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

cap.all.RZ <- aob.cap.rh.plot | aoa.cap.rh.plot | com.cap.rh.plot
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
                        legend = "right") 
                        #ncol=2, 
                        #nrow = 2,
                        #common.legend = T,
                        #legend = "right")
cap.all.RZ
setwd('D:/Fina/INRAE_Project/microservices_fig')
ggsave("cap.all.RZ.legend.tiff",
       cap.all.RZ, device = "tiff",
       width = 13, height = 4, 
       units= "in", dpi = 600)

##################################################################################################################

# ConQuR removing batch effect (in this case, treatment effect)









