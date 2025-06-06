library(tidyverse)
library(shiny)
library(leaflet)
library(sf)
library(RColorBrewer)
library(dbscan)
library(data.table)
library(reticulate)
library(here)
library(aricode)
library(cluster)
set.seed(123)
reticulate::use_virtualenv("r-reticulate")
np <- import("numpy")

#################################################################################
# Load data
################################################################################
char.dirs <- list.dirs(here("data", "synthetic_data_generated"), full.names = FALSE, recursive = FALSE)
char.dirs
char.data <- char.dirs[1]
path.data <- here("data", "synthetic_data_generated", char.data)
char.dist.measure <- "chebyshev_euclid"
int.id <- 1
matrix.org.dist.mat <- np$load(here(path.data,  char.dist.measure, "dist_mat.npy")) %>%
  as.dist()
matrix.embedding <- np$load(here(path.data, 
                                 char.dist.measure,
                                 paste0("embedding_4d_", int.id,".npy")))
sf.trips <- read_rds(here(path.data, "synthetic_flows.rds"))
sf.trips <- st_transform(sf.trips, crs = 4326)
sf.trips$flow_id <- 1:nrow(sf.trips)


#################################################################################
# Shiny
################################################################################
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .leaflet-container {
        border: 2px solid #ddd; /* Rahmen */
        border-radius: 5px; /* Abgerundete Ecken */
        margin: 10px; /* Abstand zwischen den Karten */
        padding: 5px; /* Innenabstand */
      }
      .map-container {
        display: flex;
        flex-wrap: wrap;
        justify-content: space-between;
      }
      .map-box {
        flex: 0 0 48%; /* Jede Karte nimmt 48% der Breite ein (kleiner Abstand) */
      }
    "))
  ),
  
  titlePanel("SNN Clustering mit PaCMAP"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("minpts_emb", "minPts (Emb):", min = 2, max = 30, value = 10),
      sliderInput("minpts_org", "minPts (Org):", min = 2, max = 30, value = 10),
      actionButton("go", "Cluster berechnen"),
      hr(),
      textInput("cluster_filter_org", "Cluster auswählen (org, kommagetrennt):", ""),
      textInput("cluster_filter_emb", "Cluster auswählen (emb, kommagetrennt):", ""),
      actionButton("show_cluster", "Show Cluster"),
      actionButton("reset_cluster", "Reset View"),
      hr(),
      h4("CVI & Cluster Details"),
      textOutput("ari_org"),
      textOutput("ari_emb"),
      textOutput("num_clusters_org"),
      textOutput("num_clusters_emb"),
      textOutput("num_noise_org"),
      textOutput("num_noise_emb")
    ),
    mainPanel(
      div(class = "map-container",
          div(class = "map-box",
              h4("Original-Clustering"), 
              leafletOutput("map_org", height = "400px")
          ),
          div(class = "map-box",
              h4("PaCMAP-Clustering"), 
              leafletOutput("map_emb", height = "400px")
          ),
          div(class = "map-box",
              h4("Noise with Original, but Cluster wtih PaCMAP"), 
              leafletOutput("map_noise_org", height = "400px")
          ),
          div(class = "map-box",
              h4("Noise wtih PaCMAP, but Cluster with Original"), 
              leafletOutput("map_noise_emb", height = "400px")
          )
      )
    )
    
  )
)



server <- function(input, output, session) {
  cluster_results <- eventReactive(input$go, {

    minpts_emb <- input$minpts_emb
    minpts_org <- input$minpts_org
    
    hdbscan.res.org <- hdbscan(matrix.org.dist.mat, minPts = minpts_org)
    hdbscan.res.emb <- hdbscan(matrix.embedding, minPts = minpts_emb)
    sf.trips$cluster_org <- hdbscan.res.org$cluster
    sf.trips$cluster_emb <- hdbscan.res.emb$cluster

    list(
      sf.trips = sf.trips,
      ari_org = aricode::ARI(sf.trips$cluster_id, hdbscan.res.org$cluster),
      ari_emb = aricode::ARI(sf.trips$cluster_id, hdbscan.res.emb$cluster),
      num_clusters_org = length(unique(hdbscan.res.org$cluster[hdbscan.res.org$cluster > 0])),
      num_clusters_emb = length(unique(hdbscan.res.emb$cluster[hdbscan.res.emb$cluster > 0])),
      num_noise_org = sum(hdbscan.res.org$cluster == 0),
      num_noise_emb = sum(hdbscan.res.emb$cluster == 0)
    )
  })
  
  observeEvent(input$go, {
    results <- cluster_results()
    sf.trips <- results$sf.trips
    
    output$ari_org <- renderText({ paste("Original ARI:", round(results$ari_org, 3)) })
    output$ari_emb <- renderText({ paste("Embedded ARI:", round(results$ari_emb, 3)) })
    output$num_clusters_org <- renderText({ paste("Anzahl Cluster (org):", results$num_clusters_org) })
    output$num_clusters_emb <- renderText({ paste("Anzahl Cluster (emb):", results$num_clusters_emb) })
    output$num_noise_org <- renderText({ paste("Anzahl Noise-Flows (org):", results$num_noise_org) })
    output$num_noise_emb <- renderText({ paste("Anzahl Noise-Flows (emb):", results$num_noise_emb) })
    
    output$map_org <- renderLeaflet({
      cluster_colors <- colorFactor(brewer.pal(8, "Dark2"), domain = unique(sf.trips$cluster_org[sf.trips$cluster_org > 0]))
      
      leaflet(sf.trips[sf.trips$cluster_org > 0, ]) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_org),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_org)
        )
    })
    
    output$map_emb <- renderLeaflet({
      cluster_colors <- colorFactor(brewer.pal(8, "Dark2"), domain = unique(sf.trips$cluster_emb[sf.trips$cluster_emb > 0]))
      
      leaflet(sf.trips[sf.trips$cluster_emb > 0, ]) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_emb),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_emb)
        )
    })
    
    # Map for cluster with PaCMAP but noise with original data
    output$map_noise_org <- renderLeaflet({
      sf_diff_emb <- sf.trips %>% filter(cluster_emb != 0 & cluster_org == 0)
      print(sf_diff_emb)
      cluster_colors <- colorFactor(brewer.pal(min(length(unique(sf_diff_emb$cluster_emb)), 8), "Set1"), domain = sf_diff_emb$cluster_emb)
      
      leaflet(sf_diff_emb) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_emb),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_emb)
        )
    })
    
    # Map for cluster with original data but noise with PaCMAP data
    output$map_noise_emb <- renderLeaflet({
      sf_diff_org <- sf.trips %>% filter(cluster_org != 0 & cluster_emb == 0)
      print(sf_diff_org)
      cluster_colors <- colorFactor(brewer.pal(min(length(unique(sf_diff_org$cluster_org)), 8), "Set2"), domain = sf_diff_org$cluster_org)
      
      leaflet(sf_diff_org) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_org),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_org)
        )
    })
    
  })
  
  observeEvent(input$reset_cluster, {
    results <- cluster_results()
    sf.trips <- results$sf.trips
    
    output$map_org <- renderLeaflet({
      cluster_colors <- colorFactor(brewer.pal(8, "Dark2"), domain = unique(sf.trips$cluster_org[sf.trips$cluster_org > 0]))
      
      leaflet(sf.trips[sf.trips$cluster_org > 0, ]) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_org),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_org)
        )
    })
    
    output$map_emb <- renderLeaflet({
      cluster_colors <- colorFactor(brewer.pal(8, "Dark2"), domain = unique(sf.trips$cluster_emb[sf.trips$cluster_emb > 0]))
      
      leaflet(sf.trips[sf.trips$cluster_emb > 0, ]) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_emb),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_emb)
        )
    })
    
    # Karte für Cluster in `emb`, die in `org` Noise waren
    output$map_noise_org <- renderLeaflet({
      sf_diff_emb <- sf.trips %>% filter(cluster_emb != 0 & cluster_org == 0)
      cluster_colors <- colorFactor(brewer.pal(min(length(unique(sf_diff_emb$cluster_emb)), 8), "Set1"), domain = sf_diff_emb$cluster_emb)
      
      leaflet(sf_diff_emb) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_emb),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_emb)
        )
    })
    
    # Karte für Cluster in `org`, die in `emb` Noise waren
    output$map_noise_emb <- renderLeaflet({
      sf_diff_org <- sf.trips %>% filter(cluster_org != 0 & cluster_emb == 0)
      cluster_colors <- colorFactor(brewer.pal(min(length(unique(sf_diff_org$cluster_org)), 8), "Set2"), domain = sf_diff_org$cluster_org)
      
      leaflet(sf_diff_org) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~cluster_colors(cluster_org),
          weight = 1,
          opacity = 0.8,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_org)
        )
    })
    
  })
  
  
  observeEvent(input$show_cluster, {
    sf.trips <- cluster_results()$sf.trips
    selected_clusters_org <- as.numeric(unlist(strsplit(input$cluster_filter_org, ",")))
    selected_clusters_emb <- as.numeric(unlist(strsplit(input$cluster_filter_emb, ",")))
    # cat("selected_clusters_org: ", selected_clusters_org, "\n")
    # cat("selected_clusters_emb: ", selected_clusters_emb, "\n")
    selected_flow_ids_org <- sf.trips$flow_id[sf.trips$cluster_org %in% selected_clusters_org]
    selected_flow_ids_emb <- sf.trips$flow_id[sf.trips$cluster_emb %in% selected_clusters_emb]
    # cat("selected_flow_ids_org: ", selected_flow_ids_org, "\n")
    # cat("selected_flow_ids_emb: ", selected_flow_ids_emb, "\n")
    corresponding_clusters_in_emb <- if (length(selected_clusters_org) > 0) {
      unique(sf.trips$cluster_emb[sf.trips$flow_id %in% selected_flow_ids_org])
    } else {
      numeric(0) 
    }
    
    corresponding_clusters_in_org <- if (length(selected_clusters_emb) > 0) {
      unique(sf.trips$cluster_org[sf.trips$flow_id %in% selected_flow_ids_emb])
    } else {
      numeric(0)  
    }
    # cat("corresponding_clusters_in_org: ", corresponding_clusters_in_org, "\n")
    # cat("corresponding_clusters_in_emb: ", corresponding_clusters_in_emb, "\n")
    extra_flow_ids_emb <- if (length(corresponding_clusters_in_emb) > 0) {
      sf.trips$flow_id[(sf.trips$cluster_emb %in% corresponding_clusters_in_emb & sf.trips$cluster_emb != 0) |
                         (sf.trips$cluster_emb == 0 & sf.trips$flow_id %in% selected_flow_ids_org)]
    } else {
      numeric(0)
    }
    
    extra_flow_ids_org <- if (length(corresponding_clusters_in_org) > 0) {
      sf.trips$flow_id[(sf.trips$cluster_org %in% corresponding_clusters_in_org & sf.trips$cluster_org != 0) |
                         (sf.trips$cluster_org == 0 & sf.trips$flow_id %in% selected_flow_ids_emb)]
    } else {
      numeric(0)
    }
    
    # cat("extra_flow_ids_org: ", extra_flow_ids_org, "\n")
    # cat("extra_flow_ids_emb: ", extra_flow_ids_emb, "\n")
    output$map_org <- renderLeaflet({
      sf_filtered <- sf.trips %>%
        filter(cluster_org %in% selected_clusters_org | flow_id %in% extra_flow_ids_org)
      
      if (nrow(sf_filtered) == 0) return(NULL)
      
      n_colors_org <- max(3, length(unique(sf_filtered$cluster_org)))  
      cluster_palette_org <- colorFactor(
        brewer.pal(n_colors_org, "Set2"), 
        domain = unique(sf_filtered$cluster_org)
      )
      
      leaflet(sf_filtered) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~ifelse(cluster_org == 0, "gray", cluster_palette_org(cluster_org)),
          weight = 1,
          opacity = 1,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_org)
        )
    })
    
    output$map_emb <- renderLeaflet({
      sf_filtered <- sf.trips %>%
        filter(cluster_emb %in% selected_clusters_emb | flow_id %in% extra_flow_ids_emb)
      
      if (nrow(sf_filtered) == 0) return(NULL)
      
      n_colors_emb <- max(3, length(unique(sf_filtered$cluster_emb)))  
      cluster_palette_emb <- colorFactor(
        brewer.pal(n_colors_emb, "Set2"), 
        domain = unique(sf_filtered$cluster_emb)
      )
      
      leaflet(sf_filtered) %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addPolylines(
          color = ~ifelse(cluster_emb == 0, "gray", cluster_palette_emb(cluster_emb)),
          weight = 1,
          opacity = 1,
          popup = ~paste0("Flow ID: ", flow_id, "<br>Cluster: ", cluster_emb)
        )
    })
  })
  
  
  
  
  
  
  
  
  
  
}

shinyApp(ui, server)