# calculate bray-curtis dissimilarity on rarefied data

ISS.iters.bray.B.SP <- mclapply(ISS.iters.B.SP, function(x) vegdist(x, method = "bray"))
ISS.array.bray.B.SP <- laply(ISS.iters.bray.B.SP, as.matrix)
ISS.array.bray.B.SP <- readRDS("1_DATA/ISS.array.bray.bacteria.samplingplot.rds")
ISS.bray.B.SP <- apply(ISS.array.bray.B.SP, 2:3, median)

# separate the bulk soil data

ISS.bray.B.BS <- ISS.bray.B.SP[c(1:72) , c(1:72)]
ISS.bray.B.BS <- as.dist(ISS.bray.B.BS)

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).

nc <- nrow(as.matrix(ISS.bray.B.BS))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:50) {
  cap <- CAPdiscrim(ISS.bray.B.BS ~ TxI, data = design.B.BS, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

 #here for example you would choose 14 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate

ISS.cap.B.BS <- CAPdiscrim(ISS.bray.B.BS ~ TxI, data = design.B.BS, m = 14, permutations = 9999, add = TRUE) # 97%

success <- cbind(data.frame(ISS.cap.B.BS$group), data.frame(ISS.cap.B.BS$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.B.BS$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.B.BS <- paste("CAP1 (", round((100/sum(ISS.cap.B.BS$lda.other$svd^2) * ISS.cap.B.BS$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
ISS.cap2.B.BS <- paste("CAP2 (", round((100/sum(ISS.cap.B.BS$lda.other$svd^2) * ISS.cap.B.BS$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with basic R

plot(ISS.cap.B.BS$x[, 1:2], xlab = ISS.cap1.B.BS, ylab = ISS.cap2.B.BS, pch = c(16, 17),
     col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00")[design.B.BS$TxI])
ordiellipse(ISS.cap.B.BS$x[, 1:2], groups = design.B.BS$TxI, draw = "polygon", 
            col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17),
            kind = "ehull",   border = NA, alpha = 50)
legend("bottomright", legend = c("Control BIODYN (100%)", "Drought-induced BIODYN (100%)",
                                 "Control CONFYM (100%)", "Drought-induced CONFYM (100%)",
                                 "Control CONMIN (92%)", "Drought-induced CONMIN (92%)"),
       col = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00"), pch = c(16, 17, 16, 17, 16, 17), ncol = 1, cex = 0.6)

 

# Plot with ggplot2

ggplot(as.data.frame(ISS.cap.B.BS$x), aes(x = ISS.cap.B.BS$x[,1], y = ISS.cap.B.BS$x[,2])) +
  geom_point(aes(color = design.B.BS$Treatment, shape = design.B.BS$Irrigation), size = 2) +
  xlab(ISS.cap1.B.BS) + ylab(ISS.cap2.B.BS) +
  scale_color_manual(values = colors,
                     name = "Farming system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought-induced")) + theme_classic() +
  scale_fill_manual(values = c("#009E73", "#009E73", "#FF618C", "#FF618C", "#E69F00", "#E69F00")) +
  geom_mark_ellipse(aes(fill = design.B.BS$TxI), 
                    expand = 0, linewidth = NA, show.legend = FALSE)  +
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) +
  annotate("text",x=-12,y=-15,label= "Overall reclassification rate: 97%", hjust = 0, size = 4) +
  annotate("text", x=-12, y=-16.5, label= "Pillai's test=3.1***", hjust = 0, size = 4)

 

#1.	Color code for DOK.

#BIODYN = '#009E73'
#BIOORG = '#56B4E9'
#CONFYM = '#FF618C'
#CONMIN = '#E69F00'
#NOFERT = '#C1C1C1'

#2.	Raw calculated N2O emissions are attached in the file DOK-FLUX.

#F-N2O is the calculated N2O flux. We calculated the flux using 4 timepoints 20 apart. It is the same method as mentioned in the paper of Matti Barthel, which I also cite: https://doi.org/10.1038/s41467-022-27978-6. Let me know if you have questions. 
