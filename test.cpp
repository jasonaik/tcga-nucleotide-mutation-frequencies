#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector optimized_setdiff(IntegerVector x, IntegerVector y) {
  std::unordered_set<int> y_set(y.begin(), y.end());
  IntegerVector result;
  for (int value : x) {
    if (y_set.find(value) == y_set.end()) {
      result.push_back(value);
    }
  }
  return result;
}

// [[Rcpp::export]]
List greedy_panel_curator(IntegerMatrix patient_data, int top_number) {
  std::unordered_map<int, std::unordered_set<int>> mutation_patient_map;
  for (int i = 0; i < patient_data.nrow(); ++i) {
    int patient = patient_data(i, 0);
    int mutation = patient_data(i, 1);
    mutation_patient_map[mutation].insert(patient);
  }
  
  std::vector<std::pair<int, int>> mutation_coverage;
  for (const auto& entry : mutation_patient_map) {
    mutation_coverage.push_back({entry.first, (int)entry.second.size()});
  }
  
  std::sort(mutation_coverage.begin(), mutation_coverage.end(), 
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
              return b.second < a.second;
            });
  
  IntegerVector panel_number;
  IntegerVector panel_mutation;
  
  for (int i = 0; i < top_number; ++i) {
    if (mutation_coverage.empty()) break;
    
    int best_mutation = mutation_coverage.front().first;
    panel_number.push_back(i + 1);
    panel_mutation.push_back(best_mutation);
    
    std::unordered_set<int> patients_with_best_mutation = mutation_patient_map[best_mutation];
    
    mutation_patient_map.erase(best_mutation);
    
    for (auto it = mutation_patient_map.begin(); it != mutation_patient_map.end();) {
      IntegerVector updated_patients = optimized_setdiff(IntegerVector(it->second.begin(), it->second.end()), 
                                                         IntegerVector(patients_with_best_mutation.begin(), patients_with_best_mutation.end()));
      if (updated_patients.size() == 0) {
        it = mutation_patient_map.erase(it);
      } else {
        it->second = std::unordered_set<int>(updated_patients.begin(), updated_patients.end());
        ++it;
      }
    }
    
    mutation_coverage.clear();
    for (const auto& entry : mutation_patient_map) {
      mutation_coverage.push_back({entry.first, (int)entry.second.size()});
    }
    
    std::sort(mutation_coverage.begin(), mutation_coverage.end(), 
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                return b.second < a.second;
              });
  }
  
  return List::create(Named("number") = panel_number,
                      Named("panel_mutation") = panel_mutation);
}
