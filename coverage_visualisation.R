# Coverage study HCov
## Alen Suljič (alen.suljic@mf.uni-lj.si)
### 27. 8. 2024

# libraries
library(tidyverse)

# data loading
setwd("/path/to/coverage/data")

# Coverage
## 229E
e <- read_csv("e229/coverage.csv")

### data transformation and calculations
e1 <- e %>% 
  rename(samp = 1,
         pos = 2,
         cov = 3)

### calculations
e2 <- e1 %>% 
  group_by(pos) %>% 
  mutate(median_cov_pos = median(cov),
         mean_cov_pos = mean(cov)) %>% 
  ungroup() %>% 
  mutate(global_cov = median(cov),
         consensus = if_else(cov > 9, 1, 0)) %>% 
  group_by(samp) %>% 
  mutate(median_samp_cov = median(cov),
         cons_perc = sum(consensus)/27317 * 100) %>%
  ungroup()

### visualisations
e3 <- e2 %>% 
  mutate(gene = case_when(pos <= 292 ~ "5'UTR",
                          pos >= 293 & pos <= 20568 ~ "ORF1ab",
                          pos >= 20569 & pos <= 24091 ~ "S",
                          pos >= 24092 & pos <= 24748 ~ "4ab",
                          pos >= 24749 & pos <= 24983 ~ "E",
                          pos >= 24984 & pos <= 25672 ~ "M",
                          pos >= 25673 & pos <= 26855 ~ "N",
                          pos >= 26856 & pos <= 27317 ~ "3'UTR")) %>% 
  group_by(samp,  gene) %>% 
  mutate(median_cov_gene_samp = median(cov),
         mean_cov_gene_samp = mean(cov)) %>% 
  ungroup()

#### palette
pal_229e <- c(`5'UTR` = "#000000", ORF1ab = "#3cb44b", S = "#4363d8", '4ab' = "#808000",
              E = "#ffe119", M = "#42d4f4", N = "#f032e6", `3'UTR` = "#000000")

#### plotting
mean_cov <- e2$global_cov[1]

e3 %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  ggplot() +
  geom_line(aes(pos, median_cov_pos, color = gene)) +
  geom_hline(aes(yintercept = global_cov), color = "red", lty = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  scale_x_continuous(limits = c(-450, 27600)) +
  scale_y_continuous(n.breaks = 15, limits = c(-2100, 24000)) +
  scale_color_manual(values = pal_229e) +
  annotate("segment", x = 0, xend = 292, y = -500, yend = -500, color = "#000000", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 293, xend = 20568, y = -500, yend = -500, color = "#3cb44b", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 20569, xend = 24091, y = -500, yend = -500, color = "#4363d8", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 24092, xend = 24748, y = -500, yend = -500, color = "#808000", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 24749, xend = 24983, y = -500, yend = -500, color = "#ffe119", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 24984, xend = 25672, y = -500, yend = -500, color = "#42d4f4", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 25673, xend = 26855, y = -500, yend = -500, color = "#f032e6", linewidth = 4, alpha = 0.7) +
  annotate("segment", x = 26856, xend = 27317, y = -500, yend = -500, color = "#000000", linewidth = 4, alpha = 0.7) +
  annotate("text", x = c(-350, 10910, 22500, 25300, 26250, 27600), y = -1200, label = c("5'UTR", "ORF1ab", "S", "M", "N", "3'UTR"), color = "black", size = 4, fontface = 1) +
  annotate("text", x = 24400, y = -1200, label = c("4ab"), color = "black", size = 3, fontface = 2) +
  annotate("text", x = 24850, y = -1700, label = c("E"), color = "black", size = 3.5, fontface = 1) +
  labs(x = "hCoV-229e reference genome (NC_002645.1)",
       y = "Median genome coverage per base")

ggsave("coverage_229e.tiff", device = "tiff",
       path = "/path/to/output/directory",
       width = 24, height = 18, units = "cm", dpi = 300)

## NL63
n <- read_csv("nl63/coverage.csv")

### data transformation and calculations
n1 <- n %>% 
  rename(samp = 1,
         pos = 2,
         cov = 3)

### calculations
n2 <- n1 %>% 
  group_by(pos) %>% 
  mutate(median_cov_pos = median(cov),
         mean_cov_pos = mean(cov)) %>% 
  ungroup() %>% 
  mutate(global_cov = median(cov),
         consensus = if_else(cov > 9, 1, 0)) %>% 
  group_by(samp) %>% 
  mutate(mean_samp_cov = median(cov),
         cons_perc = sum(consensus)/27553 * 100) %>%
  ungroup()

### visualizations
n3 <- n2 %>% 
  mutate(gene = case_when(pos <= 286 ~ "5'UTR",
                          pos >= 287 & pos <= 20475 ~ "ORF1ab",
                          pos >= 20476 & pos <= 24542 ~ "S",
                          pos >= 24543 & pos <= 25219 ~ "ORF3",
                          pos >= 25220 & pos <= 25433 ~ "E",
                          pos >= 25434 & pos <= 26122~ "M",
                          pos >= 26123 & pos <= 27266 ~ "N",
                          pos >= 27267 & pos <= 27553 ~ "3'UTR")) %>% 
  group_by(samp,  gene) %>% 
  mutate(median_cov_gene_samp = median(cov),
         mean_cov_gene_samp = mean(cov)) %>% 
  ungroup()

#### palette
pal_nl63 <- c(`5'UTR` = "#000000", ORF1ab = "#3cb44b", S = "#4363d8", ORF3 = "#f58231",
              E = "#ffe119", M = "#42d4f4", N = "#f032e6", `3'UTR` = "#000000")

#### plotting
n3 %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  ggplot() +
  geom_line(aes(pos, median_cov_pos, color = gene)) +
  geom_hline(aes(yintercept = global_cov), color = "red", lty = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  scale_x_continuous(limits = c(-450, 28100)) +
  scale_y_continuous(n.breaks = 15, limits = c(-2100, 27000)) +
  scale_color_manual(values = pal_nl63) +
  annotate("segment", x = 0, xend = 286, y = -500, yend = -500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 287, xend = 20475, y = -500, yend = -500, color = "#3cb44b", size = 4, alpha = 0.7) +
  annotate("segment", x = 20476, xend = 24542, y = -500, yend = -500, color = "#4363d8", size = 4, alpha = 0.7) +
  annotate("segment", x = 24543, xend = 25219, y = -500, yend = -500, color = "#f58231", size = 4, alpha = 0.7) +
  annotate("segment", x = 25220, xend = 25433, y = -500, yend = -500, color = "#ffe119", size = 4, alpha = 0.7) +
  annotate("segment", x = 25434, xend = 26122, y = -500, yend = -500, color = "#42d4f4", size = 4, alpha = 0.7) +
  annotate("segment", x = 26123, xend = 27266, y = -500, yend = -500, color = "#f032e6", size = 4, alpha = 0.7) +
  annotate("segment", x = 27267, xend = 27553, y = -500, yend = -500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("text", x = c(-350, 10910, 22400, 25750, 26700, 28100), y = -1200, label = c("5'UTR", "ORF1ab", "S", "M", "N", "3'UTR"), color = "black", size = 4, fontface = 1) +
  annotate("text", x = 24700, y = -1200, label = c("ORF3"), color = "black", size = 3, fontface = 1) +
  annotate("text", x = 25300, y = -1800, label = c("E"), color = "black", size = 3.5, fontface = 2) +
  labs(x = "hCoV-NL63 reference genome (NC_005831.2)", y = "Median genome coverage per base")

ggsave("coverage_nl63.tiff", device = "tiff",
       path = "/path/to/output/directory",
       width = 24, height = 18, units = "cm", dpi = 300)

## OC43
o <- read_csv("oc43/coverage.csv")

### data transformation and calculations
o1 <- o %>% 
  rename(samp = 1,
         pos = 2,
         cov = 3)

### calculations
o2 <- o1 %>% 
  group_by(pos) %>% 
  mutate(median_cov_pos = median(cov),
         mean_cov_pos = mean(cov)) %>% 
  ungroup() %>% 
  mutate(global_cov = median(cov),
         consensus = if_else(cov > 9, 1, 0)) %>% 
  group_by(samp) %>% 
  mutate(mean_samp_cov = median(cov),
         cons_perc = sum(consensus)/30741 * 100) %>%
  ungroup()

o3 <- o2 %>% 
  mutate(gene = case_when(pos <= 209 ~ "5'UTR",
                          pos >= 210 & pos <= 21496 ~ "ORF1ab",
                          pos >= 21497 & pos <= 22342 ~ "NS2",
                          pos >= 22343 & pos <= 23628 ~ "HE",
                          pos >= 23629 & pos <= 27704 ~ "S",
                          pos >= 27705 & pos <= 27791 ~ "ITR",
                          pos >= 27792 & pos <= 28121 ~ "ORF5",
                          pos >= 28122 & pos <= 28362 ~ "E",
                          pos >= 28363 & pos <= 29069 ~ "M",
                          pos >= 29070 & pos <= 30423 ~ "N",
                          pos >= 30424 & pos <= 30741 ~ "3'UTR")) %>% 
  group_by(samp,  gene) %>% 
  mutate(median_cov_gene_samp = median(cov),
         mean_cov_gene_samp = mean(cov)) %>% 
  ungroup()

#### palette
pal_oc43 <- c(`5'UTR` = "#000000", ORF1ab = "#3cb44b", NS2 = "#ffd8b1", HE = "#e6194B", S = "#4363d8",
              ITR = "#000000", ORF5 = "#9AFF9A", E = "#ffe119", M = "#42d4f4", N = "#f032e6", `3'UTR` = "#000000")

#### plotting
o3  %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  ggplot() +
  geom_line(aes(pos, median_cov_pos, color = gene)) +
  geom_hline(aes(yintercept = global_cov), color = "red", lty = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  scale_x_continuous(limits = c(-450, 31500)) +
  scale_y_continuous(n.breaks = 15, limits = c(-2100, 22000)) +
  scale_color_manual(values = pal_oc43) +
  annotate("segment", x = 0, xend = 209, y = -500, yend = -500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 210, xend = 21496, y = -500, yend = -500, color = "#3cb44b", size = 4, alpha = 0.7) +
  annotate("segment", x = 21497, xend = 22342, y = -500, yend = -500, color = "#ffd8b1", size = 4, alpha = 0.7) +
  annotate("segment", x = 22343, xend = 23628, y = -500, yend = -500, color = "#e6194B", size = 4, alpha = 0.7) +
  annotate("segment", x = 23629, xend = 27704, y = -500, yend = -500, color = "#4363d8", size = 4, alpha = 0.7) +
  annotate("segment", x = 27705, xend = 27791, y = -500, yend = -500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 27792, xend = 28121, y = -500, yend = -500, color = "#9AFF9A", size = 4, alpha = 0.7) +
  annotate("segment", x = 28122, xend = 28362, y = -500, yend = -500, color = "#ffe119", size = 4, alpha = 0.7) +
  annotate("segment", x = 28363, xend = 29069, y = -500, yend = -500, color = "#42d4f4", size = 4, alpha = 0.7) +
  annotate("segment", x = 29070, xend = 30423, y = -500, yend = -500, color = "#f032e6", size = 4, alpha = 0.7) +
  annotate("segment", x = 30424, xend = 30741, y = -500, yend = -500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("text", x = x = c(-350, 10910, 21800, 23000, 25800, 28700, 29800, 31200), y = -1200, label = c("5'UTR", "ORF1ab", "NS2", "HE", "S", "M", "N", "3'UTR"), color = "black", size = 3.5, fontface = 1) +
  annotate("text", x = 27750, y = -1200, label = label1, color = "black", size = 3, fontface = 1) +
  annotate("text", x = c(27500, 28300), y = -1800, label = c("ITR", "E"), color = "black", size = 3, fontface = 2) +
  labs(x = "hCoV-OC43 reference genome (NC_006213.1)", y = "Median genome coverage per base")

ggsave("coverage_oc43.tiff", device = "tiff",
       path = "/path/to/output/directory",
       width = 24, height = 18, units = "cm", dpi = 300)

## HKU1
h <- read_csv("HKU1/coverage.csv")

### data transformation and calculations
h1 <- h %>% 
  rename(samp = 1,
         pos = 2,
         cov = 3)

### calculations
h2 <- h1 %>% 
  group_by(pos) %>% 
  mutate(median_cov_pos = median(cov),
         mean_cov_pos = mean(cov)) %>% 
  ungroup() %>% 
  mutate(global_cov = median(cov),
         consensus = if_else(cov > 9, 1, 0)) %>% 
  group_by(samp) %>% 
  mutate(mean_samp_cov = median(cov),
         cons_perc = sum(consensus)/29926 * 100) %>%
  ungroup()

h3 <- h2 %>% 
  mutate(gene = case_when(pos <= 205 ~ "5'UTR",
                          pos >= 206 & pos <= 21753 ~ "ORF1ab",
                          pos >= 21754 & pos <= 22933 ~ "HE",
                          pos >= 22934 & pos <= 27012 ~ "S",
                          pos >= 27013 & pos <= 27380 ~ "ORF4",
                          pos >= 27381 & pos <= 27621 ~ "E",
                          pos >= 27622 & pos <= 28304 ~ "M",
                          pos >= 28305 & pos <= 29645 ~ "N",
                          pos >= 29646 & pos <= 29926 ~ "3'UTR")) %>% 
  group_by(samp,  gene) %>% 
  mutate(median_cov_gene_samp = median(cov),
         mean_cov_gene_samp = mean(cov)) %>% 
  ungroup()

#### palette
pal_hku1 <- c(`5'UTR` = "#000000", ORF1ab = "#3cb44b", HE = "#e6194B", S = "#4363d8",
              ORF4 = "#FF8247", E = "#ffe119", M = "#42d4f4", N = "#f032e6", `3'UTR` = "#000000")

#### plotting
h3 %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  ggplot() +
  geom_line(aes(pos, mean_cov_pos, color = gene)) +
  geom_hline(aes(yintercept = global_cov), color = "red", lty = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  scale_x_continuous(limits = c(-450, 30500)) +
  scale_y_continuous(n.breaks = 15, limits = c(-6000, 95000),
                     breaks = c(seq(0, 100000, 10000)),
                     labels = comma_format()) +
  scale_color_manual(values = pal_hku1) +
  annotate("segment", x = 0, xend = 205, y = -1500, yend = -1500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 206, xend = 21753, y = -1500, yend = -1500, color = "#3cb44b", size = 4, alpha = 0.7) +
  annotate("segment", x = 21754, xend = 22933, y = -1500, yend = -1500, color = "#e6194B", size = 4, alpha = 0.7) +
  annotate("segment", x = 22934, xend = 27012, y = -1500, yend = -1500, color = "#4363d8", size = 4, alpha = 0.7) +
  annotate("segment", x = 27013, xend = 27380, y = -1500, yend = -1500, color = "#FF8247", size = 4, alpha = 0.7) +
  annotate("segment", x = 27381, xend = 27621, y = -1500, yend = -1500, color = "#ffe119", size = 4, alpha = 0.7) +
  annotate("segment", x = 27622, xend = 28304, y = -1500, yend = -1500, color = "#42d4f4", size = 4, alpha = 0.7) +
  annotate("segment", x = 28305, xend = 28959, y = -1500, yend = -1500, color = "#f032e6", size = 4, alpha = 0.7) +
  annotate("segment", x = 28960, xend = 29926, y = -1500, yend = -1500, color = "#000000", size = 4, alpha = 0.7) +
  annotate("text", x = x = c(-350, 10910, 22350, 25000, 28650, 30200), y = -3800, label = c("5'UTR", "ORF1ab", "HE", "S", "N", "3'UTR"), color = "black", size = 4, fontface = 1) +
  annotate("text", x = 27250, y = -6000, label = c("ORF4"), color = "black", size = 3, fontface = 2) +
  annotate("text", x = c(27500, 27950), y = -3800, label = c("E", "M"), color = "black", size = 3.5, fontface = 1) +
  labs(x = "hCoV-HKU1 reference genome (NC_006577.2)",
       y = "Median genome coverage per base")

ggsave("coverage_hku1.tiff", device = "tiff",
       path = "/path/to/output/directory",
       width = 24, height = 18, units = "cm", dpi = 300)
