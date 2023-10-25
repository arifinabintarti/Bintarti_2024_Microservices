install.packages("parallel")
install.packages("BiodiversityR")
library(parallel)
library(BiodiversityR)
# run Bray-Curtis beta diversity on bulk soil
aoa.bulk_dist_bc
# metadata
aoa.meta.bulk
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:120) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ TxI, data = design.B.BS, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)