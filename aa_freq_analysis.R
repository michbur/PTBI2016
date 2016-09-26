library(dplyr)
library(reshape2)
library(ggplot2)


freq_nondeg <- read.csv2("./dat/freq_nondeg.csv")[, -1]
freq_deg <- read.csv2("./dat/freq_deg.csv")[, -1]

assign_gr <- function(x) {
  gr <- lapply(list(c("D, E, H, K, N, Q, R"), 
                    c("G, P, S, T, Y"),  
                    c("F, I, L, M, V, W"),  
                    c("A, C")), function(i) strsplit(i, split = ", ")) %>% 
    unlist(recursive = FALSE)
  ch_x <- as.character(x)
  sapply(ch_x, function(x_id) which(sapply(gr, function(i) x_id %in% i)))
}

mfreq_deg <- melt(freq_nondeg, variable.name = "residue", value.name = "freq") %>% 
  group_by(taxon, type, residue) %>% 
  summarise(freq = mean(freq)) %>% 
  ungroup %>% 
  mutate(gr = assign_gr(residue))

ggplot(mfreq_deg, aes(x = residue, y = freq, fill = taxon)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(gr~ type)


mfreq <- filter(freq_deg, type == "signal") %>% 
  melt(variable.name = "residue", value.name = "freq") %>% 
  mutate(taxon = factor(taxon, labels = c("other", "Plasmodium")))

plas_freq <- group_by(mfreq, taxon, residue) %>% 
  summarise(med = median(freq)) %>% 
  dcast(residue ~ taxon) %>% 
  mutate(plas = ifelse(Plasmodium - other > 0, "plas", ifelse(Plasmodium - other == 0, "eq", "other"))) %>% 
  select(residue, plas) %>% 
  mutate(plas = factor(plas, labels = c("No difference", "Typical for other\nsignal peptides", 
                                        "Typical for Plasmodium\nsignal peptides")))


ggplot(inner_join(mfreq, plas_freq), aes(x = residue, y = freq, color = taxon)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  facet_wrap(~ plas, ncol = 3, scales = "free_x") +
  scale_x_discrete("Residue") +
  scale_y_continuous("Normalized requency") +
  scale_color_discrete("Taxon")
