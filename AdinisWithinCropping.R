# AOA
aoa.asv.biodyn <- aoa.asv.bulk1[,c(3:4,11:14,21:22,27:28,35:38,45:46,
                                   51:52,59:62,69:70,75:76,83:86,
                                   93:94,99:100,107:110,117:118)]
aoa.asv.biodyn1 <- aoa.asv.biodyn[rowSums(aoa.asv.biodyn)>0,]
sort(rowSums(aoa.asv.biodyn1, na.rm = FALSE, dims = 1), decreasing = FALSE)
aoa.biodyn_dist_bc <- vegdist(t(aoa.asv.biodyn1), method = "bray")

aoa.meta.biodyn <- subset(aoa.meta.bulk, Treatment == 'BIODYN')
str(aoa.meta.biodyn)
set.seed(13)
aoa.biodyn.adonis <- adonis2(aoa.biodyn_dist_bc ~ Irrigation*Date, data=aoa.meta.biodyn, 
                           permutation=9999)
aoa.biodyn.adonis 
