library(dplyr)
library(vegan)
library(ggplot2)
library(devtools)

#devtools::install_github("hannet91/ggcor")
library(ggcor)

data("varechem", package = "vegan")
data("varespec", package = "vegan")

# Mantel.test 
library(vegan)
library(dplyr)
data("varechem")
data("varespec")
set.seed(20191224)
sam_grp <- sample(paste0("sample", 1:3), 24, replace = TRUE)
mantel01 <- fortify_mantel(varespec, varechem, group = sam_grp,
                           spec.select = list(spec01 = 1:5, 
                                              spec02 = 6:12,
                                              spec03 = 7:18,
                                              spec04 = 20:29,
                                              spec05 = 30:44),
                           mantel.fun = "mantel.randtest")
quickcor(mantel01, legend.title = "Mantel's r") + 
  geom_colour() + geom_cross() + facet_grid(rows = vars(.group))

mantel02 <- fortify_mantel(varespec, varechem, 
                           spec.select = list(1:10, 5:14, 7:22, 9:32)) %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))
quickcor(varechem, type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_diag_label() + remove_axis("x")
#> Warning: `add_diag_label()` is deprecated. Use `geom_diag_label()` instead.




aob.data.cor <- aob.meta.df.sub
aob.data.cor.bulk <- aob.data.cor[1:119,]


install.packages("microeco")
library(microeco)
library(magrittr)
data(dataset)
data(env_data_16S)
dataset$sample_table
dataset$sample_table <- data.frame(dataset$sample_table, env_data_16S[rownames(dataset$sample_table), ])
d1 <- clone(dataset)
d1$tax_table
d1$tax_table <- d1$tax_table[d1$tax_table$Phylum == "p__Proteobacteria", ]
d1$tidy_dataset()
d1$cal_betadiv()
d2 <- clone(dataset)
d2$tax_table <- d2$tax_table[d2$tax_table$Phylum == "p__Bacteroidetes", ]
d2$tidy_dataset()
d2$cal_betadiv()

# first perform mantel test
t1 <- trans_env$new(dataset = d1, env_cols = 8:15)
t1$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
t2 <- trans_env$new(dataset = d2, env_cols = 8:15)
t2$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

# extract a part of the results 
x1 <- data.frame(spec = "Proteobacteria", t1$res_mantel) %>% .[, c(1, 3, 6, 8)]
x2 <- data.frame(spec = "Bacteroidetes", t2$res_mantel) %>% .[, c(1, 3, 6, 8)]
# rename columns
colnames(x1) <- colnames(x2) <- c("spec", "env", "r", "p.value")
# generate interval data
x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine two tables
plot_table <- rbind(x1, x2)
# install ggcor following the steps (https://chiliubio.github.io/microeco_tutorial/intro.html#github-packages)
library(ggplot2)
library(ggcor)
set_scale()

g1 <- quickcor(t1$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

g1

quickcor(t1$data_env, type = "upper") + 
  geom_square() + 
  add_link(plot_table, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  geom_diag_label() + remove_axis("x")
















