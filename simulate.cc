#include <vector>
#include <string>
#include <unordered_map>

using std::string;
using std::unordered_map;
using std::vector;
using adjacency_matrix = std::vector<std::vector<size_t>>;

void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
              const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
              const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
              size_t total_processes) {
    // TODO: Implement this with MPI using total_processes
}