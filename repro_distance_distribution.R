#library(MASS)
library(tidyverse)
library(here)
library(lubridate)
library(sf)
library(dbscan)
library(leaflet)
library(RColorBrewer)
library(reticulate)
library(parallelDist)
library(scales)
library(viridis)
library(dbscan)
library(princurve)

reticulate::py_config()

use_virtualenv("r-reticulate", required = TRUE)

hdb   <- import("hdbscan")
np    <- import("numpy")
pacmap <- import("pacmap")

synth.data.path <- here("data",
                   "flows",
                   "synthetic",
                   "synthetic_euclid")
synth.data <- list.dirs(synth.data.path, full.names = FALSE, recursive = FALSE)
path.od.flow.distmat <- here("data",
                             "flows",
                             "synthetic",
                             "synthetic_euclid",
                             "elbe_ea0caa86-ee07-11ef-ac6b-e00af670e1c9",
                             "euclid")
#path.od.flow.distmat <- here("data","flows","realworld","mo_do_9","euclid")
int.id <- 1
dist.org <- np$load(file = here(path.od.flow.distmat,
                                       "dist_mat.npy"))
pacmap.emb <- np$load(here(path.od.flow.distmat, paste0("embedding_4d_", int.id,".npy")))
dist.emb <- dist(pacmap.emb)%>% as.matrix()

d.org.vec <- dist.org[lower.tri(dist.org)]
d.emb.vec <- dist.emb[lower.tri(dist.emb)]

d.org.vec <- rescale(d.org.vec, to = c(0, 1))
d.emb.vec <- rescale(d.emb.vec, to = c(0, 1))

df.distances <- bind_rows(
  data.frame(dist = d.org.vec, type = "Original"),
  data.frame(dist = d.emb.vec, type = "PaCMAP")
)

df.distances.wide <- data.frame(Original = d.org.vec,
                                PaCMAP = d.emb.vec)

rm(d.org.vec)
rm(d.emb.vec)
gc()

viridis.pal = c(viridis::viridis(n = 3)[2:3])
# viridis.pal <- c("grey80", "grey50")






peak.labels.heatmap <- data.frame(
  x = c(
     0.02 # A
    ,0.19   # B
    ,0.34  # C
    ,0.38  # D
    # ,0.30  # E
    # ,0.28  # F
    # ,0.18  # G
    # ,0.33  # H
    # ,0.51  # I
    # ,0.55  # J
  ),
  y = c(
     0.055 # A
    ,0.48  # B
    ,0.62  # C
    ,0.73  # D
    # ,0.71  # E
    # ,0.76  # F
    # ,0.62  # G
    # ,0.78  # H
    # ,0.85  # I
    # ,0.91   #J
  ),
  label = c(
    "A",
    "B",
    "C",
    "D"
    # "E",
    # "F",
    # "G",
    # "H",
    # "I",
    # "J"
  )
)


ggplot(df.distances.wide) +
  stat_density_2d(
    aes(x = Original, y = PaCMAP, fill = after_stat(level)),
    geom    = "polygon",
    n       = 300,                    # höhere Auflösung (300×300 statt 100×100)
    h       = c(0.05, 0.05),          # kleinere Bandbreite in x- und y-Richtung
    contour = TRUE,                   # zeichnet geschlossene Konturen
    alpha   = 0.4
  ) +
  #geom_point(aes(x = Original, y = PaCMAP), size = 0.5, alpha = 0.01) +
  annotate(
    "polygon",
    x     = c(0, 0, 1),
    y     = c(0, 1, 1),
    fill  = "grey90",  # sehr helles Grau
    alpha = 0.2
  ) +
  # 2) Polygon für untere Hälfte in mittlerem Grau
  annotate(
    "polygon",
    x     = c(0, 1, 1),
    y     = c(0, 0, 1),
    fill  = "grey60",  # etwas dunkleres Grau
    alpha = 0.2
  ) +
  annotate(
    "text",
    x     = 0.07,
    y     = 0.88,
    label = "Repulse",
    size  = 6,
    color = "grey55"
  ) +
  annotate(
    "text",
    x     = 0.93,
    y     = 0.05,
    label = "Attract",
    size  = 6,
    color = "grey40"   # etwas dunkler für guten Kontrast
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Original distances", y = "PaCMAP distances") +
  scale_fill_viridis_c(option = "D", 
                       name = "Density",
                       labels = scales::label_number(accuracy = 1)) +  # Viridis-Skala
  geom_text(
    data = peak.labels.heatmap,
    aes(x = x, y = y, label = label),
    parse = TRUE,
    inherit.aes = FALSE,
    size = 8,
    fontface = "bold"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 16L),
    legend.text = element_text(size = 14L),
    axis.text = element_text(size = 18L),
    axis.title.x = element_text(size = 20L, margin = margin(t = 12)),
    axis.title.y = element_text(size = 20L, margin = margin(r = 12)),
    strip.text = element_text(size = 20L),
    plot.margin = margin(t = 10, r = 12, l = 10)
  )

ggsave("./output/distance_heatmap.pdf", device = "pdf", width = 8, height = 8)
df.vlines <- data.frame(
  pos   = c(peak.labels.heatmap$x, peak.labels.heatmap$y),
  label = rep(peak.labels.heatmap$label, 2),
  which = rep(c("x", "y"), each = nrow(peak.labels.heatmap))
)





df.vlines <- data.frame(
  pos   = c(peak.labels.heatmap$x, peak.labels.heatmap$y),
  label = rep(peak.labels.heatmap$label, 2),
  which = rep(c("x", "y"), each = nrow(peak.labels.heatmap))
)

labels <- unique(peak.labels.heatmap$label)
pal <- colorFactor(
  palette = brewer.pal(length(labels), "Set1"),
  domain  = labels
)
label.colors <- setNames(pal(labels), labels)
label.colors

type.colors  <- viridis.pal

all.colors <- c(
  Original = type.colors[1],
  PaCMAP   = type.colors[2],
  label.colors
)

peak.labels.density <- data.frame(
  x = c(0.01, 0.05,   # A[P], A[O]
        0.21, 0.48,   # B[O], B[P]
        0.3, 0.62, # C[O], C[P]
        0.36, 0.74  # D[O], D[P]
  ),
  y = c(0.26, 0.19,    # A[P], A[O]
        0.88, 0.8,  # B[O], B[P]
        0.96, 0.92,  # C[O], C[P]
        1.03, 1.03   # D[O], D[P]
  ),
  label = c(
    "A[P]", "A[O]",
    "B[O]", "B[P]",
    "C[O]", "C[P]",
    "D[O]", "D[P]"
  ),
  stringsAsFactors = FALSE
)

df.arrows <- peak.labels.density %>%
  mutate(letter = sub("\\[.*$", "", label),
         type   = sub("^.*\\[(.*)\\]$", "\\1", label)) %>%
  select(letter, type, x, y) %>%
  pivot_wider(
    names_from  = type,
    values_from = c(x, y)
  ) %>%
  rename(
    x_start = x_O,  y_start = y_O,
    x_end   = x_P,  y_end   = y_P
  )

# df.arrows[1, c("x_start", "y_start", "x_end", "y_end")] <-
#   df.arrows[1, c("x_end",   "y_end",   "x_start", "y_start")]

df.arrows.adj <- df.arrows %>%
  rowwise() %>%
  mutate(
    dx   = x_end - x_start,
    dy   = y_end - y_start,
    dist = sqrt(dx^2 + dy^2),          
    ux = dx / dist,
    uy = dy / dist,
    x_start_adj = x_start + ux * 0.03,
    y_start_adj = y_start + uy * 0.03,
    x_end_adj   = x_end   - ux * 0.03,
    y_end_adj   = y_end   - uy * 0.03
  ) %>%
  ungroup()

# df.arrows.adj[1, "x_end"] <- 0.0025
# df.arrows.adj[1, "x_start"] <- 0.025
# df.arrows.adj[1, "y_end"] <- 
# df.arrows.adj[1, "y_start"] <- 0.21

# delta_x <- 0.00005   # um 0.05 nach links
# delta_y <- 0.00005   # um 0.05 nach oben
# df.arrows.adj$x_start_adj[1] <- df.arrows.adj$x_start_adj[1] - delta_x
# df.arrows.adj$y_start_adj[1] <- df.arrows.adj$y_start_adj[1] + delta_y
# 4) Jetzt rufen wir ggplot auf und setzen genau diese eine scale_color_manual:
ggplot(df.distances, aes(x = dist, fill = type, color = type)) +
  geom_density(aes(y = after_stat(scaled)),
               alpha = 0.7, color = "black", linewidth = .5) +
  labs(x = "Distance", y = "Density") +
  scale_fill_manual(values = viridis.pal) +   # fill für densities
  scale_color_manual(values = all.colors) +   # eine gemeinsame Farbskala
  # für beide: 'type' und 'label'
  geom_text(
    data      = peak.labels.density,
    aes(x = x, y = y, label = label),
    parse     = TRUE,
    inherit.aes= FALSE,
    size      = 7,
    fontface  = "bold"
  ) +
  # geom_vline(
  #   data      = df.vlines,
  #   aes(xintercept = pos, color = label, linetype = which),
  #   linewidth = 1
  # ) +
  geom_segment(
    data        = df.arrows.adj,
    inherit.aes = FALSE,
    aes(
      x    = x_start_adj,
      y    = y_start_adj,
      xend = x_end_adj,
      yend = y_end_adj
    ),
    arrow      = arrow(length = unit(0.01, "npc"), type = "closed"),
    color      = "black",
    linewidth  = 0.4
  )+
  scale_linetype_manual(values = c(x = "solid", y = "dashed"), guide = "none") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position     = "bottom",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 20L),
    axis.text           = element_text(size = 18L),
    axis.title.x        = element_text(size = 20L, margin = margin(t = 12)),
    axis.title.y        = element_text(size = 20L, margin = margin(r = 12)),
    strip.text          = element_text(size = 20L)
  )


ggsave("./output/distance_distribution.pdf", device = "pdf", width = 8, height = 8)

















