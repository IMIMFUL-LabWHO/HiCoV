# Variant data prep (iVar)
## 19.08.2024
### Alen Suljič (alen.suljic@mf.uni-lj.si)

# libraries
library(tidyverse)

# data loading
setwd("/path/to/sleek_variants.tsv")
d <- read_tsv("sleek_variants.tsv")

d1 <- d %>% 
  rename(sample = 1) %>% 
# delete this line if you want to look at all variants, regardless of Fishers test
  filter(PASS == TRUE) %>% 
  mutate(mut_level = if_else(ALT_FREQ >= 0.5, "major", "minor"),
         gene = str_extract(GFF_FEATURE, ".*(?=:)"),
         alt2 = if_else(grepl("-", ALT), "-", ALT)) %>% 
  unite("ref_POS", REF:POS, sep = "", remove = FALSE) %>% 
  unite("nt_mut", ref_POS, sep = "", alt2, remove = FALSE) %>% 
  select(-ref_POS, GFF_FEATURE) %>% 
# check deletions ALT and alt2!
  mutate(mut_type = case_when(grepl("-", ALT) ~ "Deletion",
                              grepl("\\+", ALT) ~ "Insertion",
                              !is.na(REF_AA) ~ "SNP"),
         mut_type = if_else(is.na(mut_type), "Non-coding", mut_type),
         snp_type = if_else(!REF_AA == ALT_AA, "Non-synonymous", "Synonymous"),
         mut = if_else(is.na(snp_type), mut_type, snp_type),
         across(c(sample, nt_mut, REF, ALT, REF_CODON, REF_AA, ALT_CODON, ALT_AA, mut_level, gene, alt2,
                  mut_type, snp_type, mut), factor)) %>%
  distinct(nt_mut, samp, .keep_all = TRUE)

summary(d1)

write_csv(d1, "path/to/output/variants_enhanced.csv")