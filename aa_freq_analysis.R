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

