#include "cayley_graph.h"

#include <H5Cpp.h>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace multiply_gates {
namespace {
struct CayleyGraphRaw {
  std::vector<uint64_t> my_hash;
  std::vector<uint8_t> my_level;
  std::vector<uint8_t> other_gate_id;
  std::vector<uint64_t> other_hash;
  std::vector<uint64_t> same_gate_ptr;
  std::vector<uint64_t> in_gate_ptr;
  std::vector<uint8_t> id_to_invid;
  std::vector<std::string> id_to_name;
};

/**
 * Read a 1D dataset from an HDF5 file into a vector.
 * @param file The HDF5 file.
 * @param dataset_name The name of the dataset.
 * @param datatype The datatype of the dataset.
 * @param vec The vector to read the dataset into.
 * @param err The error stream. Used iff the return value is false.
 * @return True if the dataset was read successfully, false otherwise.
 */
template <typename T>
bool readHdf5Vec(H5::H5File &file, const char *dataset_name,
                 H5::DataType datatype, std::vector<T> &vec,
                 std::ostream &err) {
  H5::DataSet dataset = file.openDataSet(dataset_name);
  H5::DataSpace dataspace = dataset.getSpace();
  const int ndims = dataspace.getSimpleExtentNdims();
  if (ndims != 1) {
    err << "Dataset " << dataset_name << " is not 1D." << std::endl;
    return false;
  }
  hsize_t dims_out[1];
  dataspace.getSimpleExtentDims(dims_out, NULL);
  vec.resize(dims_out[0]);
  dataset.read(vec.data(), datatype);
  return true;
}

bool readHdf5Vec(H5::H5File &file, const char *dataset_name,
                 std::vector<std::string> &vec, std::ostream &err) {
  H5::DataSet dataset =
      file.openDataSet("id_to_name"); // TODO: avoid hardcoding
  H5::DataSpace dataspace = dataset.getSpace();
  const int ndims = dataspace.getSimpleExtentNdims();
  if (ndims != 1) {
    err << "Dataset " << dataset_name << " is not 1D." << std::endl;
    return false;
  }
  hsize_t dims[1];
  dataspace.getSimpleExtentDims(dims, NULL);
  hsize_t size = dims[0];
  std::vector<char *> rdata(size);
  H5::StrType strdatatype(0, H5T_VARIABLE);
  strdatatype.setCset(dataset.getStrType().getCset());
  try {
    dataset.read(rdata.data(), strdatatype);
  } catch (const H5::Exception &e) {
    err << "Failed to read dataset \"" << dataset_name
        << "\": " << e.getDetailMsg() << std::endl;
    return false;
  }
  vec.resize(size);
  for (hsize_t i = 0; i < size; ++i) {
    vec[i] = std::string(rdata[i]);
    free(rdata[i]);
  }

  return true;
}

bool readCayleyGraphRaw(const std::string &filename,
                        CayleyGraphRaw &cayley_graph, std::ostream &err) {
  H5::H5File file(filename, H5F_ACC_RDONLY);
  return readHdf5Vec(file, "my_hash", H5::PredType::NATIVE_UINT64,
                     cayley_graph.my_hash, err) &&
         readHdf5Vec(file, "my_level", H5::PredType::NATIVE_UINT8,
                     cayley_graph.my_level, err) &&
         readHdf5Vec(file, "other_gate_id", H5::PredType::NATIVE_UINT8,
                     cayley_graph.other_gate_id, err) &&
         readHdf5Vec(file, "other_hash", H5::PredType::NATIVE_UINT64,
                     cayley_graph.other_hash, err) &&
         readHdf5Vec(file, "same_gate_ptr", H5::PredType::NATIVE_UINT64,
                     cayley_graph.same_gate_ptr, err) &&
         readHdf5Vec(file, "in_gate_ptr", H5::PredType::NATIVE_UINT64,
                     cayley_graph.in_gate_ptr, err) &&
         readHdf5Vec(file, "id_to_invid", H5::PredType::NATIVE_UINT8,
                     cayley_graph.id_to_invid, err) &&
         readHdf5Vec(file, "id_to_name", cayley_graph.id_to_name, err);
}

bool checkGateIds(std::size_t num_gates, const std::vector<uint8_t> &gate_ids,
                  std::ostream &err) {
  for (uint8_t gate_id : gate_ids) {
    if (gate_id >= num_gates) {
      err << "Invalid gate ID " << static_cast<int>(gate_id) << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace

bool CayleyGraph::read(const std::string &filename, std::ostream &err) {
  CayleyGraphRaw graph_raw;
  if (!readCayleyGraphRaw(filename, graph_raw, err)) {
    return false;
  }

  constexpr uint32_t kInvalidId = static_cast<uint32_t>(-1);
  const std::size_t num_vertices = graph_raw.my_hash.size();
  if (num_vertices * 2 > kInvalidId) {
    err << "Too many vertices: " << num_vertices * 2 << " > " << kInvalidId
        << std::endl;
    return false;
  }
  if (graph_raw.my_level.size() != num_vertices ||
      graph_raw.same_gate_ptr.size() != num_vertices * 2 ||
      graph_raw.in_gate_ptr.size() != num_vertices * 2) {
    err << "Inconsistent numbers of vertices." << std::endl;
    return false;
  }

  const std::size_t num_edges = graph_raw.other_hash.size();
  if (num_edges > kInvalidId) {
    err << "Too many edges: " << num_edges << " > " << kInvalidId << std::endl;
    return false;
  }
  if (graph_raw.other_gate_id.size() != num_edges) {
    err << "Inconsistent numbers of edges." << std::endl;
    return false;
  }

  const std::size_t num_gates = graph_raw.id_to_name.size();
  if (num_gates >
      static_cast<std::size_t>(std::numeric_limits<uint8_t>::max()) + 1) {
    err << "Too many gates: " << num_gates << " > 256" << std::endl;
    return false;
  }
  if (graph_raw.id_to_invid.size() != num_gates) {
    err << "Inconsistent numbers of gates." << std::endl;
    return false;
  }

  my_hash = std::move(graph_raw.my_hash);
  my_level = std::move(graph_raw.my_level);
  other_gate_id = std::move(graph_raw.other_gate_id);
  id_to_invid = std::move(graph_raw.id_to_invid);
  id_to_name = std::move(graph_raw.id_to_name);

  if (!(checkGateIds(num_gates, other_gate_id, err) &&
        checkGateIds(num_gates, id_to_invid, err))) {
    return false;
  }
  // Inverse of inverse is the original gate:
  for (uint8_t i = 0; i < num_gates; ++i) {
    if (id_to_invid[id_to_invid[i]] != i) {
      err << "Invalid gate inversion." << std::endl;
      return false;
    }
  }
  const std::size_t num_ptrs = 2 * num_vertices;
  for (std::size_t i = 0; i < num_ptrs; ++i) {
    if (graph_raw.same_gate_ptr[i] > num_edges ||
        graph_raw.in_gate_ptr[i] > num_edges) {
      err << "Invalid edge pointer." << std::endl;
      return false;
    }
  }
  same_gate_ptr.resize(num_ptrs);
  in_gate_ptr.resize(num_ptrs);
  for (std::size_t i = 0; i < num_ptrs; ++i) {
    same_gate_ptr[i] = static_cast<uint32_t>(graph_raw.same_gate_ptr[i]);
    in_gate_ptr[i] = static_cast<uint32_t>(graph_raw.in_gate_ptr[i]);
  }

  // Compute vertex_hash_to_id and root_id.
  root_id = kInvalidId;
  max_level = 0;
  for (uint32_t i = 0; i < my_hash.size(); ++i) {
    // We check for duplicates in my_hash.
    auto it = vertex_hash_to_id.find(my_hash[i]);
    if (it != vertex_hash_to_id.end()) {
      err << "Duplicate hash " << my_hash[i] << " at vertices " << it->second
          << " and " << i << std::endl;
      return false;
    }
    vertex_hash_to_id[my_hash[i]] = i;
    uint8_t cur_level = my_level[i];
    max_level = std::max(max_level, cur_level);
    if (cur_level == 0) {
      if (root_id != kInvalidId) {
        err << "Multiple roots at vertices " << root_id << " and " << i
            << std::endl;
        return false;
      }
      root_id = i;
    }
  }
  if (root_id == kInvalidId) {
    err << "No root found." << std::endl;
    return false;
  }

  other_id.resize(graph_raw.other_hash.size());
  for (uint32_t i = 0; i < graph_raw.other_hash.size(); ++i) {
    auto it = vertex_hash_to_id.find(graph_raw.other_hash[i]);
    if (it == vertex_hash_to_id.end()) {
      err << "Hash " << graph_raw.other_hash[i] << " not found." << std::endl;
      return false;
    }
    other_id[i] = it->second;
  }

  // my_level should decrease by 1 when following the 'in' edges and stay the
  // same when following the 'same' edges.
  for (uint32_t i = 0; i < my_hash.size(); ++i) {
    auto start = in_gate_ptr[2 * i];
    auto end = in_gate_ptr[2 * i + 1];
    if (end < start || (end == start && i != root_id)) {
      err << "At least one incoming edge expected for vertex " << i
          << std::endl;
      return false;
    }
    auto cur_level = my_level[i];
    for (uint32_t j = start; j < end; ++j) {
      if (my_level[other_id[j]] != cur_level - 1) {
        err << "Invalid level for in edge " << j << " of vertex " << i
            << std::endl;
        return false;
      }
    }
    start = same_gate_ptr[2 * i];
    end = same_gate_ptr[2 * i + 1];
    if (end < start) {
      err << "The number of horizontal edges for vertex " << i
          << " is negative." << std::endl;
      return false;
    }
    for (uint32_t j = start; j < end; ++j) {
      if (my_level[other_id[j]] != cur_level) {
        err << "Invalid level for same edge " << j << " of vertex " << i
            << std::endl;
        return false;
      }
    }
  }
  return true;
}

void CayleyGraph::fill_downstream_ids() {
  std::size_t num_vertices = my_hash.size();
  downstream_id_ptr.resize(num_vertices + 1);
  std::size_t num_ds_ids = 0;
  for (std::size_t i = 0; i < num_vertices; ++i) {
    downstream_id_ptr[i] = num_ds_ids;
    uint8_t cur_level = my_level[i];
    if (cur_level <= 1) {
      continue;
    }
    num_ds_ids += cur_level - 1;
  }
  downstream_id_ptr[num_vertices] = num_ds_ids;
  if (num_ds_ids > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("Too many downstream IDs.");
  }
  downstream_ids.resize(num_ds_ids);

  // v_stack contains pairs (vertex_id, needs_expanding)
  // describing the currently enqueued vertices. The vertices are in the order
  // of non-increasing level and the main purpose of the stack is to ensure that
  // the downstream vertices are processed before the current vertex.
  std::vector<std::pair<uint32_t, bool>> v_stack;
  std::unordered_set<uint32_t> visited;
  constexpr uint32_t kInvalidId = static_cast<uint32_t>(-1);
  for (std::size_t i = 0; i < num_vertices; ++i) {
    uint8_t cur_level = my_level[i];
    if (cur_level <= 1 || visited.count(i) > 0) {
      continue;
    }
    v_stack.emplace_back(i, cur_level > 2);
    visited.insert(i);
    while (v_stack.size() > 0) {
      auto [v, needs_expanding] = v_stack.back();
      uint8_t cur_level = my_level[v];
      auto igp_start = in_gate_ptr[2 * v];
      auto igp_end = in_gate_ptr[2 * v + 1];
      if (needs_expanding) {
        // We need to ensure that downstream vertices are processed before v.
        v_stack.back().second = false;
        bool next_needs_expanding = cur_level > 3;
        for (uint32_t j = igp_start; j < igp_end; ++j) {
          uint32_t u = other_id[j];
          if (visited.count(u) == 0) {
            v_stack.emplace_back(u, next_needs_expanding);
            visited.insert(u);
          }
        }
      } else {
        // Downstream vertices have been processed, we can now process v.
        v_stack.pop_back();
        uint32_t ds_start = downstream_id_ptr[v];
        if (!(igp_end > igp_start)) {
          throw std::runtime_error("At least one incoming edge expected.");
        }
        uint32_t u = other_id[igp_start];
        uint32_t u_ds_start = downstream_id_ptr[u];
        for (uint8_t j = 0; j < cur_level - 2; ++j) {
          downstream_ids[ds_start + j] = u_ds_start + j;
        }
        if (igp_end == igp_start + 1) {
          downstream_ids[ds_start + cur_level - 2] = u;
        } else {
          downstream_ids[ds_start + cur_level - 2] = kInvalidId;
        }
        for (uint32_t j = igp_start + 1; j < igp_end; ++j) {
          uint32_t u = other_id[j];
          uint32_t u_ds_start = downstream_id_ptr[u];
          for (uint8_t k = 0; k < cur_level - 2; ++k) {
            // If currently identified bottleneck for v is different from the
            // bottleneck for u, then there is no bottleneck for v.
            if (downstream_ids[ds_start + k] !=
                downstream_ids[u_ds_start + k]) {
              downstream_ids[ds_start + k] = kInvalidId;
            }
          }
        }
      }
    }
  }
}

CycleTops CayleyGraph::find_cycle_tops() const {
  // 1. List all vertices of level >= 1 without bottlenecks.
  // 2. List all horizontal edges without common bottlenecks.
  CycleTops cycle_tops;
  // Most vertices are cycle tops (I guess); the root is not:
  cycle_tops.top_vertex_ids.reserve(my_hash.size() - 1);
  std::size_t num_h_edges = 0;
  constexpr uint32_t kInvalidId = static_cast<uint32_t>(-1);
  for (std::size_t i = 0; i < my_hash.size(); ++i) {
    uint8_t cur_level = my_level[i];
    if (cur_level == 0) {
      continue;
    }
    num_h_edges += same_gate_ptr[2 * i + 1] - same_gate_ptr[2 * i];
    bool has_bottleneck = false;
    const uint32_t *cur_ds_ids = downstream_ids.data() + downstream_id_ptr[i];
    // Note: unsigned comparison won't work for cur_level <= 0, but
    // cur_level >= 0 and we checked for 0 above:
    for (std::size_t j = 0; j < std::size_t(cur_level - 1); ++j) {
      if (cur_ds_ids[j] != kInvalidId) {
        has_bottleneck = true;
        break;
      }
    }
    if (!has_bottleneck) {
      cycle_tops.top_vertex_ids.push_back(i);
    }
  }
  cycle_tops.top_edges.reserve(num_h_edges);
  for (std::size_t v = 0; v < my_hash.size(); ++v) {
    uint8_t cur_level = my_level[v];
    if (cur_level == 0) {
      continue;
    }
    const uint32_t *v_ds_ids = downstream_ids.data() + downstream_id_ptr[v];
    std::span<const uint32_t> v_same_ids(other_id.data() + same_gate_ptr[2 * v],
                                         same_gate_ptr[2 * v + 1] -
                                             same_gate_ptr[2 * v]);
    std::span<const uint8_t> v_same_gates(
        other_gate_id.data() + same_gate_ptr[2 * v],
        same_gate_ptr[2 * v + 1] - same_gate_ptr[2 * v]);
    for (std::size_t j = 0; j < v_same_ids.size(); ++j) {
      uint32_t u = v_same_ids[j];
      const uint32_t *u_ds_ids = downstream_ids.data() + downstream_id_ptr[u];
      bool has_common_bottleneck = false;
      // cur_level == 0 was checked above.
      for (std::size_t k = 0; k < std::size_t(cur_level - 1); ++k) {
        if (v_ds_ids[k] != kInvalidId && v_ds_ids[k] == u_ds_ids[k]) {
          has_common_bottleneck = true;
          break;
        }
      }
      if (!has_common_bottleneck) {
        cycle_tops.top_edges.emplace_back(v, u, v_same_gates[j]);
      }
    }
  }
  return cycle_tops;
}

void CycleBuilder::_fill_0() {
  uint32_t v0 = vertex_ids[0][pos[0]];
  std::span<const uint32_t> v0_in_ids(
      cayley_graph.other_id.data() + cayley_graph.in_gate_ptr[2 * v0],
      cayley_graph.in_gate_ptr[2 * v0 + 1] - cayley_graph.in_gate_ptr[2 * v0]);
  uint32_t edge_id = 0;
  // If there are multiple edge candidates, pick one:
  if (v0_in_ids.size() > 1) {
    uint32_t v1 = vertex_ids[1][pos[1]];
    std::span<const uint32_t> v1_ds_ids(
        cayley_graph.downstream_ids.data() + cayley_graph.downstream_id_ptr[v1],
        cayley_graph.downstream_ids.data() +
            cayley_graph.downstream_id_ptr[v1 + 1]);

    edge_buf.resize(v0_in_ids.size());
    std::iota(edge_buf.begin(), edge_buf.end(), 0);
    // Try all downstream edges in a random order until a valid one is found.
    while (edge_buf.size() > 1) {
      uint32_t j0 =
          std::uniform_int_distribution<int>(0, edge_buf.size() - 1)(rng);
      uint32_t j = edge_buf[j0];
      uint32_t u = v0_in_ids[j];
      bool is_disjoint =
          (u != v1_ds_ids.back()) && _check_disjoint(u, v1_ds_ids);
      if (is_disjoint) {
        edge_id = j;
        break;
      }
      edge_buf[j0] = edge_buf.back();
      edge_buf.pop_back();
    }
    if (edge_buf.size() == 1) {
      edge_id = edge_buf[0];
    }
  }
  pos[0] -= 1;
  vertex_ids[0][pos[0]] = v0_in_ids[edge_id];
  gate_ids[0][pos[0]] =
      cayley_graph.other_gate_id[cayley_graph.in_gate_ptr[2 * v0] + edge_id];
}

void CycleBuilder::_fill_1() {
  uint32_t v1 = vertex_ids[1][pos[1]];
  std::span<const uint32_t> v1_in_ids(
      cayley_graph.other_id.data() + cayley_graph.in_gate_ptr[2 * v1],
      cayley_graph.in_gate_ptr[2 * v1 + 1] - cayley_graph.in_gate_ptr[2 * v1]);
  uint32_t edge_id = 0;
  // If there are multiple edge candidates, pick one:
  if (v1_in_ids.size() > 1) {
    uint32_t v0 = vertex_ids[0][pos[0]];
    std::span<const uint32_t> v0_ds_ids(
        cayley_graph.downstream_ids.data() + cayley_graph.downstream_id_ptr[v0],
        cayley_graph.downstream_ids.data() +
            cayley_graph.downstream_id_ptr[v0 + 1]);

    edge_buf.resize(v1_in_ids.size());
    std::iota(edge_buf.begin(), edge_buf.end(), 0);
    // Try all downstream edges in a random order until a valid one is found.
    while (edge_buf.size() > 1) {
      uint32_t j0 =
          std::uniform_int_distribution<int>(0, edge_buf.size() - 1)(rng);
      uint32_t j = edge_buf[j0];
      uint32_t u = v1_in_ids[j];
      bool is_disjoint = (u != v0) && _check_disjoint(u, v0_ds_ids);
      if (is_disjoint) {
        edge_id = j;
        break;
      }
      edge_buf[j0] = edge_buf.back();
      edge_buf.pop_back();
    }
    if (edge_buf.size() == 1) {
      edge_id = edge_buf[0];
    }
  }
  pos[1] -= 1;
  vertex_ids[1][pos[1]] = v1_in_ids[edge_id];
  gate_ids[1][pos[1]] =
      cayley_graph.other_gate_id[cayley_graph.in_gate_ptr[2 * v1] + edge_id];
}

bool CycleBuilder::_check_disjoint(uint32_t u,
                                   std::span<const uint32_t> downstream_ids) {
  std::span<const uint32_t> u_ds_ids(cayley_graph.downstream_ids.data() +
                                         cayley_graph.downstream_id_ptr[u],
                                     cayley_graph.downstream_ids.data() +
                                         cayley_graph.downstream_id_ptr[u + 1]);
  constexpr uint32_t kInvalidId = static_cast<uint32_t>(-1);
  for (std::size_t i = 0; i < u_ds_ids.size(); ++i) {
    if (u_ds_ids[i] != kInvalidId && u_ds_ids[i] == downstream_ids[i]) {
      return false;
    }
  }
  return true;
}

namespace {
constexpr char kOpenQasm3QubitHeader[] =
    "OPENQASM 3;\ninclude \"stdgates.inc\";\nqubit[3] q;";
constexpr char kOpenQasm3GateSxdg[] = "\ngate sxdg q { u2(3*pi/2,pi/2) q; }";
} // namespace

/**
 * Initialize class members corresponding to a given gate ID.
 *
 * This is an auxiliary function for OpenQasmGenerator constructor.
 * It accepts a partially initialized OpenQasmGenerator object and,
 * thus, its implementation should be understood in the context of the
 * constructor.
 *
 * @param generator A partially constructed OpenQasmGenerator object to
 * initialize.
 * @param gate_str_to_qasm_lookup A map from gate names to their OpenQASM names.
 * @param i The gate ID.
 * @param gate_str A string describing the gate (e.g., "cs_01").
 *
 * Precondition:
 * * `id_to_qasm` and `is_sxdg` are resized to the number of gates
 *   but have an unspecified value at index `i`.
 * * `id_to_tokens_ptr` has length `i` and is filled with indices in `tokens`
 *   corresponding to the previous values of `i`.
 * * `id_to_tokens_data` contains the lists of tokens corresponding to the
 *   previous values of `i`.
 * * `gate_name_to_token` maps gate names (e.g. "cs", "x") to their token ids
 *   for previously seen gates.
 * * `qubit_j_to_token` is fully initialized.
 *
 * Postcondition:
 * * `id_to_qasm[i]` is set to the OpenQASM representation of the gate.
 * * `is_sxdg[i]` is set to true if the gate is an sxdg gate.
 * * One element is added to `id_to_tokens_ptr` corresponding to the previous
 *   size of `id_to_tokens_data`.
 * * List of tokens encoding `gate_str` are added to `id_to_tokens_data`.
 * * Newly created tokens (if any) are added to `tokens` and
 *   `gate_name_to_token`.
 */
void OpenQasmGenerator::initialize_gate_info(
    const std::unordered_map<std::string, std::string> &gate_str_to_qasm_lookup,
    std::size_t i, const std::string &gate_str) {
  // Comments annotate parsing of 2 examples: gate_str = "cs_01" and "x2".
  std::string gate;        // We want to set `gate` to "cs", "x".
  std::string qubits_part; // We want to set `qubits_part` to "01", "2".
  size_t pos = gate_str.find('_');
  if (pos != std::string::npos) {
    gate = gate_str.substr(0, pos);         // "cs"
    qubits_part = gate_str.substr(pos + 1); // "01"
  } else {
    gate = gate_str.substr(0, gate_str.size() - 1); // "x"
    qubits_part = gate_str.back();                  // "2"
  }
  for (char c : qubits_part) {
    if (!isdigit(c)) {
      std::ostringstream err;
      err << "Unknown gate name format: \"" << gate_str << "\".";
      throw std::runtime_error(err.str());
    }
  }

  // gate_str is now split into `gate` and `qubits_part`.
  uint32_t num_qubits = qubits_part.size();

  std::ostringstream qasm_oss;
  // We add '\n' at the beginning to simplify concatenation into OpenQASM
  // script:
  qasm_oss << "\n";
  {
    auto it = gate_str_to_qasm_lookup.find(gate);
    if (it != gate_str_to_qasm_lookup.end()) {
      gate = it->second;
    }
  }
  qasm_oss << gate << " ";
  TokenId gate_token_id;
  {
    auto it = gate_name_to_token.find(gate);
    if (it == gate_name_to_token.end()) {
      gate_token_id = tokens.size();
      gate_name_to_token[gate] = gate_token_id;
      tokens.emplace_back(TokenData::Type::kGate, gate, num_qubits);
    } else {
      gate_token_id = it->second;
      TokenData &token_data = tokens[gate_token_id];
      if (token_data.type != TokenData::Type::kGate) {
        std::ostringstream err;
        err << "Token \"" << gate << "\" is not a gate.";
        throw std::runtime_error(err.str());
      }
      if (token_data.text != gate) {
        std::ostringstream err;
        err << "Token \"" << gate << "\" maps to " << gate_token_id
            << " with text \"" << token_data.text << "\".";
        throw std::runtime_error(err.str());
      }
      if (token_data.data != num_qubits) {
        std::ostringstream err;
        err << "Token \"" << gate << "\" maps to " << gate_token_id << ". ";
        err << "There is a mismatch in the number of qubits: " << num_qubits
            << " vs " << token_data.data << ".";
        throw std::runtime_error(err.str());
      }
    }
  }
  id_to_tokens_ptr.push_back(id_to_tokens_data.size());
  id_to_tokens_data.push_back(gate_token_id);
  for (size_t i = 0; i < qubits_part.size(); ++i) {
    if (i > 0) {
      qasm_oss << ", ";
    }
    qasm_oss << "q[" << qubits_part[i] << "]";
    id_to_tokens_data.push_back(qubit_j_to_token[qubits_part[i] - '0']);
  }
  qasm_oss << ";";
  id_to_qasm[i] = qasm_oss.str();
  is_sxdg[i] = gate_str.starts_with("sxdg");
}

using namespace std::literals::string_literals;

OpenQasmGenerator::OpenQasmGenerator(
    const std::vector<std::string> &id_to_name) {
  tokens.clear();
  // Note: the order should match kTokenId* constants:
  tokens.emplace_back(TokenData::Type::kSpecial, "\0"s, 0);
  tokens.emplace_back(TokenData::Type::kSpecial, "MULTIPLY_GATES\n", 0);
  tokens.emplace_back(TokenData::Type::kSpecial, "ANSWER:\n", 0);

  constexpr uint32_t kMaxNumQubits = 20;
  qubit_j_to_token.resize(kMaxNumQubits);
  for (uint32_t j = 0; j < kMaxNumQubits; ++j) {
    qubit_j_to_token[j] = tokens.size();
    tokens.emplace_back(TokenData::Type::kQubit, "q[" + std::to_string(j) + "]",
                        j);
  }

  gate_name_to_token.clear();
  id_to_tokens_ptr.clear();
  id_to_tokens_data.clear();

  std::unordered_map<std::string, std::string> gate_str_to_qasm_lookup{
      {"cs", "ctrl @ s"}, {"csdg", "ctrl @ sdg"}};
  std::size_t num_gate_ids = id_to_name.size();

  id_to_qasm.resize(num_gate_ids);
  is_sxdg.resize(num_gate_ids);

  for (size_t i = 0; i < id_to_name.size(); ++i) {
    const std::string &gate_str = id_to_name[i];
    initialize_gate_info(gate_str_to_qasm_lookup, i, gate_str);
  }
  id_to_tokens_ptr.push_back(id_to_tokens_data.size());
}

void OpenQasmGenerator::cycle_to_openqasm(
    const std::vector<uint8_t> &cycle_gates, std::ostream &out) const {
  out << kOpenQasm3QubitHeader;

  bool has_sxdg = false;
  for (auto c : cycle_gates) {
    if (is_sxdg[c]) {
      has_sxdg = true;
      break;
    }
  }

  if (has_sxdg) {
    out << kOpenQasm3GateSxdg;
  }

  for (auto c : cycle_gates) {
    out << id_to_qasm[c];
  }
}

void OpenQasmGenerator::safe_cycle_to_tokens(
    std::span<const uint8_t> cycle_gates, std::vector<TokenId> &out) const {
  out.clear();
  for (std::size_t gate_id : cycle_gates) {
    if (gate_id + 1 >= id_to_tokens_ptr.size()) {
      throw std::runtime_error("Invalid gate ID.");
    }
    const TokenId *start = id_to_tokens_data.data() + id_to_tokens_ptr[gate_id];
    const TokenId *end =
        id_to_tokens_data.data() + id_to_tokens_ptr[gate_id + 1];
    if (end > id_to_tokens_data.data() + id_to_tokens_data.size()) {
      throw std::runtime_error("Invalid token pointer.");
    }
    out.insert(out.end(), start, end);
  }
}

void OpenQasmGenerator::cycle_to_tokens(std::span<const uint8_t> cycle_gates,
                                        std::vector<TokenId> &out) const {
  out.clear();
  for (std::size_t gate_id : cycle_gates) {
    auto start = id_to_tokens_data.begin() + id_to_tokens_ptr[gate_id];
    auto end = id_to_tokens_data.begin() + id_to_tokens_ptr[gate_id + 1];
    out.insert(out.end(), start, end);
  }
}

void OpenQasmGenerator::detokenize(const std::vector<TokenId> &in,
                                   std::ostream &out) const {
  if (in.size() == 0) {
    return;
  }
  out << tokens[in[0]].text;
  TokenData::Type state = tokens[in[0]].type;
  for (std::size_t pos = 1; pos < in.size(); ++pos) {
    const TokenId token_id = in[pos];
    const TokenData &token_data = tokens[token_id];
    const TokenData::Type new_state = token_data.type;
    switch (2 * (state == TokenData::Type::kQubit) +
            (new_state == TokenData::Type::kQubit)) {
    case 0:
      // Non-qubit to non-qubit transition.
      if (state == TokenData::Type::kGate) {
        out << ";\n";
      }
      break;
    case 1:
      // Non-qubit to qubit transition.
      if (state == TokenData::Type::kGate) {
        out << " ";
      }
      break;
    case 2:
      // Qubit to non-qubit transition.
      out << ";\n";
      break;
    case 3:
      // Qubit to qubit transition.
      out << ", ";
      break;
    }
    out << token_data.text;
    state = new_state;
  }
  if (state == TokenData::Type::kGate || state == TokenData::Type::kQubit) {
    out << ";";
  }
}

void OpenQasmGenerator::shuffle_qubits(std::mt19937_64 &rng,
                                       std::span<TokenId> token_ids) const {
  std::vector<TokenId> available_qubit_tids(qubit_j_to_token);
  std::vector<TokenId> qubit_map(qubit_j_to_token.size(), kTokenIdEos);
  for (auto &token_id : token_ids) {
    const TokenData &token_data = tokens[token_id];
    if (token_data.type != TokenData::Type::kQubit) {
      continue;
    }
    const uint32_t token_j = token_data.data;
    TokenId new_token_id = qubit_map[token_j];
    if (new_token_id == kTokenIdEos) {
      std::size_t l = std::uniform_int_distribution<int>(
          0, qubit_j_to_token.size() - 1)(rng);
      new_token_id = qubit_map[token_j] = available_qubit_tids[l];
      available_qubit_tids[l] = available_qubit_tids.back();
      available_qubit_tids.pop_back();
    }
    token_id = new_token_id;
  }
}

int32_t gen_multiply_gates_tasks(const OpenQasmGenerator &generator,
                                 CycleBuilder &cycle_builder,
                                 const CycleTops &cycle_tops, SplitType split,
                                 int32_t num_tasks, int32_t max_size,
                                 int32_t *part_lengths, TokenId *task_tokens) {
  if (num_tasks == 0) {
    return 0;
  }
  if (num_tasks < 0) {
    return -1;
  }
  if (max_size < 5) {
    return -2;
  }
  std::vector<uint8_t> cycle_gates(
      static_cast<std::size_t>((max_size + 3) / 4));
  std::vector<TokenId> tokens_buf(static_cast<std::size_t>(max_size));
  for (int i = 0; i < num_tasks; ++i) {
    int32_t *cur_part_lengths = part_lengths + i * 2;
    TokenId *task_tokens_start = task_tokens + i * max_size;
    TokenId *cur_task_tokens = task_tokens_start;
    cycle_builder.from_random_cycle_top(cycle_tops, split);
    int32_t cg_size = cycle_builder.get_cycle_gates(true, cycle_gates.data(),
                                                    cycle_gates.size());
    if (cg_size <= 0) {
      return -3;
    }
    *(cur_task_tokens++) = kTokenIdMultiplyGates;
    generator.cycle_to_tokens(
        std::span<const uint8_t>(cycle_gates).first(cg_size - 1), tokens_buf);
    if (tokens_buf.size() > static_cast<std::size_t>(max_size - 3)) {
      return -4;
    }
    for (TokenId token_id : tokens_buf) {
      *(cur_task_tokens++) = token_id;
    }
    *(cur_task_tokens++) = kTokenIdAnswer;
    *(cur_part_lengths++) = 2 + tokens_buf.size();
    int32_t remaining_size = max_size - 2 - tokens_buf.size();
    generator.cycle_to_tokens(std::span(&cycle_gates[cg_size - 1], 1),
                              tokens_buf);
    // 1 token is needed for '\0' at the end.
    if (tokens_buf.size() >= static_cast<std::size_t>(remaining_size)) {
      return -5;
    }
    for (TokenId token_id : tokens_buf) {
      *(cur_task_tokens++) = token_id;
    }
    *(cur_part_lengths++) = 1 + tokens_buf.size();
    generator.shuffle_qubits(cycle_builder.rng,
                             std::span(task_tokens_start, cur_task_tokens));
    // Fill the remaining space with '\0'.
    remaining_size -= tokens_buf.size();
    while (remaining_size-- > 0) {
      *(cur_task_tokens++) = kTokenIdEos;
    }
  }
  return 0;
}

} // namespace multiply_gates
