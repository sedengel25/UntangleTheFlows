#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <map>
#include <mutex>
#include <chrono>
#include <fstream>

using namespace Rcpp;
using namespace RcppParallel;

// Hilfsfunktion zur Erstellung der Netzwerk-Map
std::map<int, std::pair<int, int>> create_network_map(DataFrame dt_network) {
  IntegerVector ids = dt_network["id"];
  IntegerVector starts = dt_network["source"];
  IntegerVector ends = dt_network["target"];
  std::map<int, std::pair<int, int>> network_map;
  for (int i = 0; i < ids.size(); i++) {
    network_map[ids[i]] = std::make_pair(starts[i], ends[i]);
  }
  return network_map;
}

// Funktion zum Laden der dist_map
void initialize_map(std::map<std::pair<int, int>, int>& dist_map, DataFrame dt_dist_mat) {
  IntegerVector sources = dt_dist_mat["source"];
  IntegerVector targets = dt_dist_mat["target"];
  IntegerVector distances = dt_dist_mat["m"];
  for (int i = 0; i < sources.size(); ++i) {
    int source = sources[i];
    int target = targets[i];
    int distance = distances[i];
    dist_map[std::make_pair(std::min(source, target), std::max(source, target))] = distance;
  }
}

// Funktion zum Abrufen der Distanz
int get_distance(const std::map<std::pair<int, int>, int>& dist_map, int source, int target) {
  auto it = dist_map.find(std::make_pair(std::min(source, target), std::max(source, target)));
  return (it != dist_map.end()) ? it->second : -1;
}

// Struktur für parallele Verarbeitung
struct NetworkProcessor : public Worker {
  const IntegerVector& od_pts_full_id;
  const IntegerVector& od_pts_full_line_id;
  const IntegerVector& od_pts_full_start;
  const IntegerVector& od_pts_full_end;
  const std::map<int, std::pair<int, int>>& network_map;
  const std::map<std::pair<int, int>, int>& dist_map;
  
  RMatrix<int> distance_matrix;
  
  NetworkProcessor(const IntegerVector& od_pts_full_id,
                   const IntegerVector& od_pts_full_line_id,
                   const IntegerVector& od_pts_full_start,
                   const IntegerVector& od_pts_full_end,
                   const std::map<int, std::pair<int, int>>& network_map,
                   const std::map<std::pair<int, int>, int>& dist_map,
                   RMatrix<int> distance_matrix)
  : od_pts_full_id(od_pts_full_id),
  od_pts_full_line_id(od_pts_full_line_id),
  od_pts_full_start(od_pts_full_start),
  od_pts_full_end(od_pts_full_end),
  network_map(network_map),
  dist_map(dist_map),
  distance_matrix(distance_matrix) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int point_ij = od_pts_full_id[i];
      int line_ij = od_pts_full_line_id[i];
      int start_ij = od_pts_full_start[i];
      int end_ij = od_pts_full_end[i];
      
      auto it_ij = network_map.find(line_ij);
      if (it_ij == network_map.end()) continue;
      
      int source_ij = it_ij->second.first;
      int target_ij = it_ij->second.second;
      
      
      
      for (std::size_t j = i; j < od_pts_full_id.size(); ++j) {
        int point_kl = od_pts_full_id[j];
        int line_kl = od_pts_full_line_id[j];
        int start_kl = od_pts_full_start[j];
        int end_kl = od_pts_full_end[j];
        if (point_ij == point_kl) continue;
        
        auto it_kl = network_map.find(line_kl);
        if (it_kl == network_map.end()) continue;
        
        int source_kl = it_kl->second.first;
        int target_kl = it_kl->second.second;
        
        int nd_ik = get_distance(dist_map, source_ij, source_kl);
        int nd_jl = get_distance(dist_map, target_ij, target_kl);
        int nd_il = get_distance(dist_map, source_ij, target_kl);
        int nd_jk = get_distance(dist_map, source_kl, target_ij);
        
        if (nd_ik == -1 || nd_jl == -1 || nd_il == -1 || nd_jk == -1) continue;
        
        int om_on_distance = std::min({
          start_ij + start_kl + nd_ik,
          end_ij + end_kl + nd_jl,
          start_ij + end_kl + nd_il,
          end_ij + start_kl + nd_jk
        });
        
        // Atomar gespeicherte Werte für parallelen Zugriff
        distance_matrix(i, j) = om_on_distance;
        distance_matrix(j, i) = om_on_distance;
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix parallel_process_networks_matrix(DataFrame dt_od_pts_full,
                                               DataFrame dt_network,
                                               DataFrame dt_dist_mat,
                                               int num_cores) {
  auto start = std::chrono::system_clock::now();
  
  static std::map<std::pair<int, int>, int> dist_map;
  initialize_map(dist_map, dt_dist_mat);
  
  IntegerVector od_pts_full_id = dt_od_pts_full["flow_id"];
  IntegerVector od_pts_full_line_id = dt_od_pts_full["id_edge"];
  IntegerVector od_pts_full_start = dt_od_pts_full["dist_to_start"];
  IntegerVector od_pts_full_end = dt_od_pts_full["dist_to_end"];
  
  auto network_map = create_network_map(dt_network);
  
  // Maximale ID für Matrixgröße bestimmen
  int max_id = max(od_pts_full_id);
  IntegerMatrix distance_matrix_R(max_id, max_id);
  std::fill(distance_matrix_R.begin(), distance_matrix_R.end(), 0);
  
  // Umwandlung in RMatrix<int>
    RMatrix<int> distance_matrix(distance_matrix_R);
  
  std::size_t chunk_size = (od_pts_full_id.size() + num_cores - 1) / num_cores;
  
  // Paralleler Worker mit RMatrix<int>
    NetworkProcessor worker(od_pts_full_id, od_pts_full_line_id, od_pts_full_start, 
                            od_pts_full_end, network_map, dist_map, distance_matrix);
  
  parallelFor(0, od_pts_full_id.size(), worker, chunk_size);
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  Rcpp::Rcout << "Time needed for network distance calculation: " << elapsed_seconds.count() << " seconds." << std::endl;
  
  return distance_matrix_R;  // `RMatrix<int>` zurück in `IntegerMatrix` konvertieren
}