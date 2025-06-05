library(tidyverse)
library(here)
#################################################################################
# Load data
################################################################################
dt.euclid <- read_rds(here("res_exp_final","dt_euclid.rds"))
dt.cheb <-  read_rds(here("res_exp_final","dt_chebyshev.rds"))
dt.network <-  read_rds(here("res_exp_final", "dt_network.rds"))
dt.euclid$dist_measure <- "Manhattan-Euclidean"
dt.cheb$dist_measure <- "Chebyshev-Euclidean"
dt.network$dist_measure <- "Manhattan-Network"

dt <- rbind(dt.euclid, dt.cheb, dt.network)

#################################################################################
# Compare all Algorithms on all CVI
################################################################################
viridis_pal <- c(viridis::viridis(n = 4)[2:4])
ggplot(dt, aes(x = algo_lab, y = median_cvi_val,
               fill = algo_lab, alpha = type_lab)) +
  geom_boxplot(outlier.shape = 1
               ,
               outlier.size= 1,
               position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = viridis_pal) +
  scale_alpha_manual(values = c("Original" = 1, "PaCMAP" = 0.5),
                     guide = guide_legend(override.aes = list(fill = c("grey40", "grey70"))))+
  facet_grid(cvi_name_lab ~ dist_measure, scales = "free_y") +
  ylab("CVI value\n") +
  labs(fill = "Algorithm:", alpha = "Data:") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(vjust = 0.8, size = 14L),
    legend.text = element_text(size = 12L),
    axis.text = element_text(size = 12L),
    axis.title = element_text(size = 14L),
    strip.text = element_text(size = 12L),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#################################################################################
# HDBSCAN
################################################################################
dt.hdbscan.ari <- dt %>%
  filter(algo == "hdbscan", cvi_name == "ari", dist_measure == "Manhattan-Euclidean")

viridis_pal <- c(viridis::viridis(n = 3)[2:3])
ggplot(data = dt.hdbscan.ari) +
  geom_rect(data = data.frame(xmin = seq(1.5, 48.5, 2), 
                              xmax = seq(2.5, 49.5, 2), 
                              ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightgrey", alpha = 0.3) +
  geom_boxplot(mapping = aes(x = as.factor(minpts), 
                             y = median_cvi_val, fill = type_lab),
               outlier.size = 0.5) +
  geom_vline(xintercept = seq(1.5, 48.5, 1), color = "darkgrey", linewidth = 0.2) +
  scale_x_discrete(breaks = seq(2, 49, 1), 
                   labels = ifelse(seq(2, 49, 1) %% 6 == 4, seq(2, 49, 1), "")) +
  #scale_fill_viridis_d("Data Type") +
  scale_fill_manual(values = viridis_pal, name = "Data") +
  xlab("\nminpts") +
  ylab("ARI value\n") +
  # facet_grid(ari ~ data_type) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(vjust = 0.8, size = 14L),
    legend.text = element_text(size = 12L),
    axis.text = element_text(size = 12L),
    axis.title = element_text(size = 14L),
    strip.text = element_text(size = 12L)
  )





