#ifndef CAYLEY_GRAPH_H
#define CAYLEY_GRAPH_H

#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace multiply_gates {

struct CycleTops {
  // Vertex can be the top of a cycle if it is not the root and does not have
  // any bottlenecks on lower levels.
  std::vector<uint32_t> top_vertex_ids;

  // The edge is described by tuple (i, j, g), where U_i = gate_g @ U_j and i
  // and j are vertex IDs on the same level.
  std::vector<std::tuple<uint32_t, uint32_t, uint8_t>> top_edges;
};

struct CayleyGraph {
  // For vertex i, we have the following information:
  // * my_hash[i] is the hash of the unitary corresponding to the vertex.
  // * my_level[i] is the level of the vertex in the Cayley graph (i.e. the
  //   distance from the root).
  // * other_gate_id[start:end] for
  //   start = in_gate_ptr[2 * i]
  //   end = in_gate_ptr[2 * i + 1]
  //   are the gate IDs incoming to the vertex from the previous level.
  // * other_gate_id[start:end] for
  //   start = same_gate_ptr[2 * i]
  //   end = same_gate_ptr[2 * i + 1]
  //   are the gate IDs at the same level as the vertex.
  // * other_id[start:end] for start and end as above are the corresponding
  //   vertex IDs.
  std::vector<uint64_t> my_hash;
  std::vector<uint8_t> my_level;
  std::vector<uint8_t> other_gate_id;
  std::vector<uint32_t> other_id;
  std::vector<uint32_t> same_gate_ptr;
  std::vector<uint32_t> in_gate_ptr;

  uint32_t root_id;  // The vertex ID of the root of the Cayley graph.
  uint8_t max_level; // max(my_level)

  // One can look up a vertex index i given its hash using vertex_hash_to_id:
  std::unordered_map<uint64_t, uint32_t> vertex_hash_to_id;

  // If i is a vertex index and j in other_id[start:end] for
  // start = in_gate_ptr[2 * i], end = in_gate_ptr[2 * i + 1]
  // (or same_gate_ptr), then U_i = gate_g @ U_j, where g is the gate with ID
  // in other_gate_id[start:end]. To obtain U_j from U_i, i.e. find g' s.t.
  // U_j = gate_{g'} @ U_i, one can use the gate g' = id_to_invid[g].
  std::vector<uint8_t> id_to_invid;
  std::vector<std::string> id_to_name; // Look up the gate name given its ID.

  // Ids of downstream vertices. For each vertex i, downstream_ids[start + j] is
  // the ID of the vertex on level 1 + j that is the bottleneck, i.e. a vertex
  // through which all shortest paths to the root must pass, and uint32_t(-1) if
  // there is no such vertex. Here
  // * 0 <= j < my_level[i] - 1
  // * start = downstream_id_ptr[i].
  // The chunks corresponding to different vertices are stored consecutively.
  std::vector<uint32_t> downstream_ids;
  // size is my_level.size() + 1: the last element is downstream_ids.size().
  std::vector<uint32_t> downstream_id_ptr;

public:
  CayleyGraph() : root_id(0) {}
  bool read(const std::string &filename, std::ostream &err);
  void fill_downstream_ids();
  // Requeires downstream_ids to be filled.
  CycleTops find_cycle_tops() const;
};

enum class SplitType : int32_t { kAll = 0, kTrain = 1, kTest = 2 };

// Auxiliary class for building a cycle. The class is used to combine
// multiple variations of the cycle building algorithm. The basic idea
// is to start at a certain level and build two paths down to the root
// not intersecting (other than at the beginning and the end).
//
// The class is designed to be reusable, so that the same instance can be
// re-initialized and used for building multiple cycles. In this way, one
// can build arbitrary amount of cycles using only O(1) memory allocations.
class CycleBuilder {
public:
  CycleBuilder(const CayleyGraph &cayley_graph_, std::mt19937_64 &rng_)
      : rng(rng_), cayley_graph(cayley_graph_) {
    std::size_t max_path_size =
        static_cast<std::size_t>(cayley_graph.max_level) + 1;
    vertex_ids[0].reserve(max_path_size);
    vertex_ids[1].reserve(max_path_size + 1);
    gate_ids[0].reserve(max_path_size - 1);
    gate_ids[1].reserve(max_path_size);
    // We do not reserve edge_buf since identifying its max size would
    // require iterating over all vertices.
  }
  // Maximal number of edges expected to be in a cycle.
  // Add 1 to get the maximal number of vertices (the start and end vertices
  // of the cycle are the same).
  int32_t max_cycle_size() const {
    return static_cast<int32_t>(cayley_graph.max_level) * 2 + 1;
  }

  // We have two implementations of get_cycle_gates:
  // * C-style implementation operates on raw buffers, it is designed to be
  // useful for ctypes interface.
  // * C++-style implementation operates on std::vector.
  //
  // @param invert_last If true, the last gate is not inverted. This is useful
  // for Multiply Gates task: this function generates a sequence of gates
  // which multiply to identity. When the last gate is inverted, the sequence
  // of gates before it multiplies to the last gate. In Multiply Gates task,
  // that last gate is the answer ML needs to find.
  int32_t get_cycle_gates(bool invert_last, uint8_t *cycle_gates,
                          int32_t max_size) const {
    if (max_size < 0) {
      return -1;
    }
    std::size_t num_gates_0 = gate_ids[0].size();
    std::size_t num_gates_1 = gate_ids[1].size();
    if (num_gates_0 + num_gates_1 > static_cast<std::size_t>(max_size)) {
      return -2;
    }
    std::size_t pos;
    for (pos = 0; pos < num_gates_0; ++pos) {
      cycle_gates[pos] = gate_ids[0][pos];
    }
    for (std::size_t j = num_gates_1 - 1; j > 0; --j, ++pos) {
      cycle_gates[pos] = cayley_graph.id_to_invid[gate_ids[1][j]];
    }
    if (invert_last) {
      cycle_gates[pos] = gate_ids[1][0];
    } else {
      cycle_gates[pos] = cayley_graph.id_to_invid[gate_ids[1][0]];
    }
    ++pos;
    return pos;
  }

  void get_cycle_gates(std::vector<uint8_t> &cycle_gates) const {
    std::size_t num_gates_0 = gate_ids[0].size();
    cycle_gates.resize(0);
    cycle_gates.reserve(num_gates_0 + gate_ids[1].size());
    for (std::size_t i = 0; i < num_gates_0; ++i) {
      cycle_gates.push_back(gate_ids[0][i]);
    }
    for (std::size_t j = gate_ids[1].size(); j > 0; --j) {
      cycle_gates.push_back(cayley_graph.id_to_invid[gate_ids[1][j - 1]]);
    }
  }

  int32_t get_cycle_vertices(uint32_t *cycle_vertices, int32_t max_size) const {
    if (max_size < 0) {
      return -1;
    }
    std::size_t num_vertices_0 = vertex_ids[0].size();
    std::size_t num_vertices_1 = vertex_ids[1].size();
    if (num_vertices_0 + num_vertices_1 - 1 >
        static_cast<std::size_t>(max_size)) {
      return -2;
    }
    std::size_t pos;
    for (pos = 0; pos < num_vertices_0; ++pos) {
      cycle_vertices[pos] = vertex_ids[0][pos];
    }
    for (std::size_t j = num_vertices_1 - 1; j > 0; --j, ++pos) {
      cycle_vertices[pos] = vertex_ids[1][j - 1];
    }
    return pos;
  }
  void get_cycle_vertices(std::vector<uint32_t> &cycle_vertices) const {
    std::size_t num_vertices_0 = vertex_ids[0].size();
    cycle_vertices.resize(0);
    cycle_vertices.reserve(num_vertices_0 + vertex_ids[1].size() - 1);
    for (std::size_t i = 0; i < num_vertices_0; ++i) {
      cycle_vertices.push_back(vertex_ids[0][i]);
    }
    // Note that vertex_ids[1].back() == vertex_ids[0].back(),
    // we avoid adding it twice here:
    for (std::size_t j = vertex_ids[1].size() - 1; j > 0; --j) {
      cycle_vertices.push_back(vertex_ids[1][j - 1]);
    }
  }
  void from_vertex(uint32_t vertex_id) {
    uint8_t level = cayley_graph.my_level[vertex_id];
    vertex_ids[0].resize(level + 1);
    vertex_ids[1].resize(level + 1);
    gate_ids[0].resize(level);
    gate_ids[1].resize(level);
    vertex_ids[0].back() = vertex_id;
    vertex_ids[1].back() = vertex_id;
    pos[0] = level;
    pos[1] = level;
    _fill_all();
  }
  // This assumes U(u) = gate(g) @ U(v).
  void from_edge(uint32_t u, uint32_t v, uint8_t g) {
    uint8_t level = cayley_graph.my_level[u];
    if (cayley_graph.my_level[v] != level) {
      throw std::runtime_error("from_edge: u and v are not on the same level.");
    }
    vertex_ids[0].resize(level + 1);
    vertex_ids[1].resize(level + 2);
    gate_ids[0].resize(level);
    gate_ids[1].resize(level + 1);
    vertex_ids[0].back() = u;
    vertex_ids[1].back() = u;
    vertex_ids[1][level] = v;
    gate_ids[1][level] = g;
    pos[0] = level;
    pos[1] = level;
    _fill_all();
  }
  void from_random_cycle_top(const CycleTops &cycle_tops,
                             SplitType split = SplitType::kTrain) {
    std::size_t num_options =
        cycle_tops.top_vertex_ids.size() + cycle_tops.top_edges.size();
    std::size_t choice;
    // kAll: choose any option
    // kTrain: choose even option
    // kTest: choose odd option
    switch (split) {
    case SplitType::kAll:
      choice =
          std::uniform_int_distribution<std::size_t>(0, num_options - 1)(rng);
      break;
    case SplitType::kTrain:
      choice = 2 * std::uniform_int_distribution<std::size_t>(
                       0, (num_options - 1) / 2)(rng);
      break;
    case SplitType::kTest:
      choice = 1 + 2 * std::uniform_int_distribution<std::size_t>(
                           0, num_options / 2 - 1)(rng);
      break;
    default:
      throw std::runtime_error("from_random_cycle_top: unknown split.");
    }
    if (choice < cycle_tops.top_vertex_ids.size()) {
      from_vertex(cycle_tops.top_vertex_ids[choice]);
    } else {
      auto [u, v, g] =
          cycle_tops.top_edges[choice - cycle_tops.top_vertex_ids.size()];
      from_edge(u, v, g);
    }
  }

public:
  std::mt19937_64 &rng;

private:
  const CayleyGraph &cayley_graph;
  // Two paths:
  // vertex_ids[j][l] is usually at level l (in particular, vertex_ids[j][0] is
  // the root). The only exception is that the last vertex
  // vertex_ids[j][l + 1] may be present with level l in the case the top of the
  // cycle is a horizontal edge.
  // In a constructed state we have:
  // U(vertex_ids[j][l + 1]) = gate_ids[j][l] @ U(vertex_ids[j][l]).
  // vertex_ids[0].back() == vertex_ids[1].back().
  std::vector<uint32_t> vertex_ids[2];
  std::vector<uint8_t> gate_ids[2];
  // During the construction,
  // vertex_ids[j][pos[j]:end] and gate_ids[j][pos[j]:end]
  // are already filled.
  std::size_t pos[2];

  // Internal buffer for the downstream edge candidates:
  std::vector<uint32_t> edge_buf;

private:
  // _fill_j decreases pos[j] by 1 by filling one more edge.
  // These have to be called in a very specific order:
  // Initially, pos[0] == pos[1] > 0. We call _fill_0, _fill_1, _fill_0, ...
  // In this way before _fill_0 is called, pos[0] == pos[1] > 0,
  // and before _fill_1 is called, pos[1] == pos[0] + 1.
  void _fill_0();
  void _fill_1();
  void _fill_all() {
    while (pos[0] > 0) {
      _fill_0();
      _fill_1();
    }
  }
  // Check if downstream vertices of u are disjoint from downstream_ids.
  bool _check_disjoint(uint32_t u, std::span<const uint32_t> downstream_ids);
};

struct TokenData {
  enum class Type : uint32_t { kSpecial = 0, kGate = 1, kQubit = 2 };
  TokenData(Type type, const std::string &text, uint32_t data)
      : type(type), text(text), data(data) {}
  Type type;
  std::string text;
  // data has the following meaning depending on type:
  // * SPECIAL: (unused)
  // * GATE: the number of qubits the gate acts on.
  // * QUBIT: the qubit index.
  uint32_t data;
};

using TokenId = uint32_t;
// Note: below constants should match OpenQasmGenerator::OpenQasmGenerator
constexpr TokenId kTokenIdEos = 0; // TokenId for '\0' is 0, which is nice.
constexpr TokenId kTokenIdMultiplyGates = 1;
constexpr TokenId kTokenIdAnswer = 2;
class OpenQasmGenerator {
public:
  explicit OpenQasmGenerator(const std::vector<std::string> &id_to_name);
  void cycle_to_openqasm(const std::vector<uint8_t> &cycle_gates,
                         std::ostream &out) const;
  void cycle_to_tokens(std::span<const uint8_t> cycle_gates,
                       std::vector<TokenId> &out) const;
  void safe_cycle_to_tokens(std::span<const uint8_t> cycle_gates,
                            std::vector<TokenId> &out) const;
  void detokenize(const std::vector<TokenId> &in, std::ostream &out) const;
  const std::vector<TokenId> &get_qubit_j_to_token() const {
    return qubit_j_to_token;
  }
  void shuffle_qubits(std::mt19937_64 &rng, std::span<TokenId> token_ids) const;

private:
  // OpenQASM text generation:
  std::vector<std::string> id_to_qasm;
  std::vector<bool> is_sxdg;
  // Token generation:
  std::vector<uint32_t> id_to_tokens_ptr;
  std::vector<TokenId> id_to_tokens_data;
  std::unordered_map<std::string, TokenId> gate_name_to_token;
  std::vector<TokenId> qubit_j_to_token;
  std::vector<TokenData> tokens;
  void initialize_gate_info(const std::unordered_map<std::string, std::string>
                                &gate_str_to_qasm_lookup,
                            std::size_t i, const std::string &gate_str);
};

int32_t gen_multiply_gates_tasks(const OpenQasmGenerator &generator,
                                 CycleBuilder &cycle_builder,
                                 const CycleTops &cycle_tops, SplitType split,
                                 int32_t num_tasks, int32_t max_size,
                                 int32_t *part_lengths, TokenId *task_tokens);

} // namespace multiply_gates
#endif // CAYLEY_GRAPH_H
