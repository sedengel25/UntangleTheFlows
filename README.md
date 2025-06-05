# Untangle the Flows 
This repository provides supplementary code and data for the paper **"Untangle the Flows: Enhanced Clustering of Origin-Destination Flows using PaCMAP"**, 
supporting the implementation and analysis presented in the study.  

The project consists of two independent parts:

1. **Instance Generation & Clustering** (`01_*` to `03_*` files)  
   Enables users to generate OD flow datasets, apply clustering methods, and visualize results interactively.

2. **Reproducibility** (`repro_*` files)  
   Provides scripts to reproduce the results from the research paper.


---

## 1. Instance Generation & Clustering

### Instance Generation (`01_instance_generator.R`)
This script enables users to generate synthetic **origin-destination (OD) flow datasets** for a specific section of Dresden. It relies on the following input files:

- `sf_network.rds`: An `sf` object containing the street network.
- `local_node_dist_mat.rds`: A  `data.table` containing network distances between all nodes in the street network.

The main parameters include `num_min_length` & `num_max_length` for flow length boundaries, `int_min_cl_size` & `int_max_cl_size` for the number of flows per cluster, and `int_min_n_cl` & `int_max_n_cl` for the possible number of clusters. Generated OD flow datasets are stored under unique UUID-based directories as subfolders of `/data/synthetic_data_generated`, each containing an `rds` file with the generated **OD flows** as an `sf` object with the associated `cluster_id`, representing the ground truth.


### Clustering (`02_hdbscan.R`)
For a selected synthetic OD dataset, this script computes OD flow distance matrices using `manhattan_euclid`, `chebyshev_euclid`, or `manhattan_network`. The chosen distance measure is then processed with **PaCMAP**, generating the following outputs:


- A `dist()` object containing the original OD flow distance matrix.
- A `matrix()` object storing the PaCMAP 4D embedding result.

The `dbscan::hdbscan()` function is applied with varying `minpts` values (ranging from 2 to 50) to evaluate clustering performance. The results are visualized using `ggplot2`, displaying a line chart that compares the ARI (Adjusted Rand Index) values across different `minpts` settings.

### Shiny-App (`03_shiny_app.R`)
This `shiny` dashboard enables users to interactively set the `minpts` parameter for `hdbscan`, dynamically updating four `leaflet` maps visualizing the corresponding cluster results for one specific synthetic dataset, based on the original OD flow distance matrix and the PaCMAP 4D embedding. Users can explore individual clusters and compare their assignments across both representations.

---

## Reproducibility
Unfortunately, due to data protection and copyright reasons, the data used for the qualitative real-world analysis cannot be provided. To enable the reproducability of the synthetic experiments conducted as well as the distance-distribution plots, the 100 synthetic OD flow datasets used can be downloaded from [Zenodo](https://zenodo.org/records/15025399?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjIwYzY5NzgwLTRlN2EtNGI3Yy04OGJkLWU2NzNmMTdlMWEyYiIsImRhdGEiOnt9LCJyYW5kb20iOiIzMDY0NGZhN2ViNzkwZmExMDE1NzQxZDA0NmRmN2E2MiJ9.FZxejH2itkRpWa-T4MRH62ZDfkAu8QtiQksBd0pzQ7rlRiXnLTkzo3kZbblljaiGSUUBA7leOD-utgRJ_PrUXw). It is recommended to unzip `experiment_synthetic.zip` in the folder `/data/synthetic_data_experiment`.


### Run the experiment (`repro_exp_main.R`)
This script applies `dbscan::dbscan()`, `dbscan::hdbscan()`, and `dbscan::snnClust()` to both the OD flow distance matrix and the PaCMAP embedding of the selected distance measure across all 100 datasets. It tracks clustering performance using ARI (`aricode::ARI()`), FMI (`clusterCrit::extCriteria()`), and SWC (`cluster::silhouette()`) for the specified parameter setting.

The experiments were executed on an HPC-CPU cluster with jobs running in parallel, each corresponding to a parameter configuration from `/input/exp_param.csv`. The provided source code is structured accordingly, meaning that running `repro_exp_main.R` once represents a single job. Since each job takes ~15 minutes and `/input/exp_param.csv` contains 6,628 configurations, parallel execution is recommended. The parameter `int_n_samples` controls how often PaCMAP is applied to a single dataset. Since `experiment_synthetic.zip` includes 10 precomputed embeddings, `int_n_samples` must be set to a value within this limit. The script generates an `rds` file storing cluster validation indices (CVI) for each dataset and distance measure. The file is saved to `/res_exp_intermediate/{distance_measure}/{UUID dataset}/cvi_h_dbscan_snn_{k}_{eps}_{minpts}.rds`.


### Process experiment results (`repro_exp_process_01.R`)
Running this script takes the median value for each CVI, algorithm and dataset over the `int_n_samples` PaCMAP-runs. The processed results are written to `/res_exp_processed/{distance_measure}/{UUID dataset}/cvi_h_dbscan_snn_grouped_{k}_{eps}_{minpts}.rds`.

### Prepare experiment results for visualization (`repro_exp_process_02.R`)
This script simply further processes the results to prepare for visualization. Results are written to `/res_exp_final/dt_{distance_measure}.rds`.

### Visualize experiment results (`repro_exp_vis.R`)
This script creates the boxplot from the paper, comparing the CVI results for all three algorithms and all three different distance measures, as well as the boxplot for all `hdbscan`-configurations. **Note:** If the experiment cannot be run, the second dataset from [Zenodo](https://zenodo.org/records/15025399?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjIwYzY5NzgwLTRlN2EtNGI3Yy04OGJkLWU2NzNmMTdlMWEyYiIsImRhdGEiOnt9LCJyYW5kb20iOiIzMDY0NGZhN2ViNzkwZmExMDE1NzQxZDA0NmRmN2E2MiJ9.FZxejH2itkRpWa-T4MRH62ZDfkAu8QtiQksBd0pzQ7rlRiXnLTkzo3kZbblljaiGSUUBA7leOD-utgRJ_PrUXw) (`experiment_results.zip`) can be extracted to `res_exp_final` to enable the visualizations.

### Reproduce distance-distribution plots (`repro_distance_distribution.R`)
Running `repro_distance_distribution.R` loads the original OD flow distance matrix and one specific PaCMAP embedding of the provided synthetic OD flow dataset and creats the distance heatmap (Figure 5) as well as the distance distribution  (Figure 6) shown in the paper. To inspect these plots for additional synthetic datasets, one of the downloaded datasets in `/data/synthetic_data_experiment` can be copied to `data/distance_distribution`. Vertical lines, as discussed in the paper to locate the peaks A–D, can optionally be included. This can be done by uncommenting the corresponding `geom_vline(...)` commands in the script. The positions of these lines (and the associated peaks) of course vary depending on the dataset used.

## Requirements
The files `R_requirements.txt` and `py_requirements.txt` list the `R` and `Python` packages, along with their respective versions, used in this project.  

**Note:** To save files as `.npy` arrays in `R`, the `reticulate` package is required. Therefore, the user has to set up a Python virtual environment (`venv`).  

