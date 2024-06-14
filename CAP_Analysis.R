### CAP Analysis ###

library(vegan)
library(BiodiversityR) # ALWAYS LOAD FROM THE R CONSOLE NOT R STUDIO!!!!
library(parallel)
library(ggplot2)
library(ggforce)

# Read the rarefied ASV table
aoa.asv.tab <- read.csv("aoa.asv.tab.csv",row.names = 1)
aoa.asv.tab[] <- lapply(aoa.asv.tab, as.numeric)
# Read the metadata
aoa.metadat <- read.csv("aoa.metadat.csv",row.names = 1)
# change the explanatory variables as factor
aoa.metadat[sapply(aoa.metadat, is.integer)] <- lapply(aoa.metadat[sapply(aoa.metadat, is.integer)], as.factor)
aoa.metadat[sapply(aoa.metadat, is.character)] <- lapply(aoa.metadat[sapply(aoa.metadat, is.character)], as.factor)
str(aoa.metadat)

# Calculate Bray-Curtis dissimilarity matrix on the rarefied ASV table
aoa.dist.bray <- vegdist(t(aoa.asv.tab), method = "bray")

# Run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.dist.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.dist.bray ~ DxF, data = aoa.metadat, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 39 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
aoa.cap <- CAPdiscrim(aoa.dist.bray ~ DxF, data = aoa.metadat, m = 39, permutations = 9999, add = TRUE) # 94.16667% 
aoa.cap

success <- cbind(data.frame(aoa.cap$group), data.frame(aoa.cap$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1 <- paste("CAP1 (", round((100/sum(aoa.cap$lda.other$svd^2) * aoa.cap$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2 <- paste("CAP2 (", round((100/sum(aoa.cap$lda.other$svd^2) * aoa.cap$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")


# Plot with basic R

plot(aoa.cap$x[, 1:2], xlab = aoa.cap1, ylab = aoa.cap2, pch = c(16, 17),
     col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")[aoa.metadat$DxF])
ordiellipse(aoa.cap$x[, 1:2], groups = aoa.metadat$DxF, draw = "polygon", 
            col = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomleft", legend = c("Control BIODYN (95%)", "Drought-induced BIODYN (90%)",
                                "Control CONFYM (100%)", "Drought-induced CONFYM (100%)",
                                "Control CONMIN (95%)", "Drought-induced CONMIN (85%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.5)

# Plot with ggplot2

aoa.cap.Plot <- ggplot(as.data.frame(aoa.cap$x), aes(x = aoa.cap$x[,1], y = aoa.cap$x[,2])) +
  geom_point(aes(color = aoa.metadat$Fertilization, shape = aoa.metadat$Drought), size = 2) +
  xlab(aoa.cap1) + ylab(aoa.cap2) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Drought",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  geom_mark_ellipse(aes(fill = aoa.metadat$DxF), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) +
  annotate("text",x=-14,y=-10,label= "Overall reclassification rate: 94.2%", hjust = 0, size = 4) +
  annotate("text", x=-14, y=-11.5, label= "Pillai's test=4.1***", hjust = 0, size = 4)
aoa.cap.Plot

