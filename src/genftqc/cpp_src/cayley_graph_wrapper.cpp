#include "cayley_graph_wrapper.h"
#include "cayley_graph.h"
#include <exception>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace multiply_gates;
namespace {
struct CycleBuilderWrapper {
  CycleBuilderWrapper(CayleyGraph &graph, uint64_t seed)
      : rng(seed), builder(graph, rng) {}
  std::mt19937_64 rng;
  CycleBuilder builder;
};
} // namespace

extern "C" {

void *create_cayley_graph() { return new CayleyGraph(); }

void destroy_cayley_graph(void *cayley_graph) {
  delete static_cast<CayleyGraph *>(cayley_graph);
}

bool cayley_graph_read(void *cayley_graph, const char *filename) {
  CayleyGraph *graph = static_cast<CayleyGraph *>(cayley_graph);
  try {
    return graph->read(filename, std::cerr);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
}

void cayley_graph_fill_downstream_ids(void *cayley_graph) {
  static_cast<CayleyGraph *>(cayley_graph)->fill_downstream_ids();
}

void *create_cycle_tops(void *cayley_graph) {
  return new CycleTops(
      static_cast<CayleyGraph *>(cayley_graph)->find_cycle_tops());
}

void destroy_cycle_tops(void *cycle_tops) {
  delete static_cast<CycleTops *>(cycle_tops);
}

int32_t cycle_tops_get_num_vertices(void *cycle_tops) {
  return static_cast<CycleTops *>(cycle_tops)->top_vertex_ids.size();
}

void cycle_tops_get_vertices(void *cycle_tops, uint32_t *vertices) {
  const std::vector<uint32_t> &v =
      static_cast<CycleTops *>(cycle_tops)->top_vertex_ids;
  std::copy(v.begin(), v.end(), vertices);
}

int32_t cycle_tops_get_num_edges(void *cycle_tops) {
  return static_cast<CycleTops *>(cycle_tops)->top_edges.size();
}

void cycle_tops_get_edges(void *cycle_tops, uint32_t *edge_u, uint32_t *edge_v,
                          uint8_t *edge_g) {
  const std::vector<std::tuple<uint32_t, uint32_t, uint8_t>> &edges =
      static_cast<CycleTops *>(cycle_tops)->top_edges;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edge_u[i] = std::get<0>(edges[i]);
    edge_v[i] = std::get<1>(edges[i]);
    edge_g[i] = std::get<2>(edges[i]);
  }
}

void *create_cycle_builder(void *cayley_graph, uint64_t seed) {
  return new CycleBuilderWrapper(*static_cast<CayleyGraph *>(cayley_graph),
                                 seed);
}

void destroy_cycle_builder(void *cycle_builder) {
  delete static_cast<CycleBuilderWrapper *>(cycle_builder);
}

int32_t cycle_builder_get_max_cycle_size(void *cycle_builder) {
  return static_cast<CycleBuilderWrapper *>(cycle_builder)
      ->builder.max_cycle_size();
}

int32_t cycle_builder_get_cycle_gates(void *cycle_builder, bool invert_last,
                                      uint8_t *cycle_gates, int32_t max_size) {
  return static_cast<CycleBuilderWrapper *>(cycle_builder)
      ->builder.get_cycle_gates(invert_last, cycle_gates, max_size);
}

int32_t cycle_builder_get_cycle_vertices(void *cycle_builder,
                                         uint32_t *cycle_vertices,
                                         int32_t max_size) {
  return static_cast<CycleBuilderWrapper *>(cycle_builder)
      ->builder.get_cycle_vertices(cycle_vertices, max_size);
}

void cycle_builder_from_vertex(void *cycle_builder, uint32_t vertex_id) {
  static_cast<CycleBuilderWrapper *>(cycle_builder)
      ->builder.from_vertex(vertex_id);
}

void cycle_builder_from_edge(void *cycle_builder, uint32_t u, uint32_t v,
                             uint8_t g) {
  static_cast<CycleBuilderWrapper *>(cycle_builder)->builder.from_edge(u, v, g);
}

void cycle_builder_from_random_cycle_top(void *cycle_builder,
                                         void *cycle_tops) {
  static_cast<CycleBuilderWrapper *>(cycle_builder)
      ->builder.from_random_cycle_top(*static_cast<CycleTops *>(cycle_tops));
}

void *create_openqasm_generator(void *cayley_graph) {
  return new OpenQasmGenerator(
      static_cast<CayleyGraph *>(cayley_graph)->id_to_name);
}

void destroy_openqasm_generator(void *generator) {
  delete static_cast<OpenQasmGenerator *>(generator);
}

int32_t cycle_to_openqasm(void *generator, const uint8_t *cycle_gates,
                          int32_t cycle_gates_size, char *out,
                          int32_t out_size) {
  if (out_size <= 0) {
    return -1;
  }
  std::ostringstream oss;
  std::vector<uint8_t> gates(cycle_gates, cycle_gates + cycle_gates_size);
  static_cast<OpenQasmGenerator *>(generator)->cycle_to_openqasm(gates, oss);
  std::string s = oss.str();
  int32_t res = s.length() + 1;
  std::size_t s_out_size = res;
  if (res > out_size) {
    res = -1;
    s_out_size = out_size;
  }
  // Note: c_str() returns a pointer to a null-terminated character array,
  // i.e. its length is s.length() + 1. Thus, below code works in both cases:
  // * if s.length() >= out_size, the first out_size characters are copied.
  // * otherwise, s.length() characters and '\0' are copied.
  std::memcpy(out, s.c_str(), s_out_size);
  return res;
}

int32_t cycle_to_tokens(void *generator, const uint8_t *cycle_gates,
                        int32_t cycle_gates_size, uint32_t *out,
                        int32_t out_size) {
  if (!generator || !cycle_gates || !out) {
    return -2;
  }

  multiply_gates::OpenQasmGenerator *gen =
      static_cast<multiply_gates::OpenQasmGenerator *>(generator);
  // TODO: remove unnecessary copy
  std::vector<uint8_t> cycle_gates_vector(cycle_gates,
                                          cycle_gates + cycle_gates_size);
  std::vector<TokenId> tokens;

  gen->cycle_to_tokens(cycle_gates_vector, tokens);

  if (tokens.size() > static_cast<size_t>(out_size)) {
    std::copy(tokens.begin(), tokens.begin() + out_size, out);
    return -1;
  } else {
    std::copy(tokens.begin(), tokens.end(), out);
    return static_cast<int32_t>(tokens.size());
  }
}

int32_t detokenize(void *generator, const uint32_t *tokens, int32_t tokens_size,
                   char *out, int32_t out_size) {
  if (!generator || !tokens || !out) {
    return -2;
  }

  multiply_gates::OpenQasmGenerator *gen =
      static_cast<multiply_gates::OpenQasmGenerator *>(generator);
  std::vector<TokenId> tokens_vector(tokens, tokens + tokens_size);
  std::ostringstream oss;

  gen->detokenize(tokens_vector, oss);
  std::string result = oss.str();

  if (result.size() >= static_cast<size_t>(out_size)) {
    std::memcpy(out, result.c_str(), out_size);
    return -1;
  } else {
    std::memcpy(out, result.c_str(), result.size() + 1);
    return static_cast<int32_t>(result.size());
  }
}

int32_t gen_multiply_gates_tasks(void *generator, void *cycle_builder,
                                 void *cycle_tops, int32_t split,
                                 int32_t num_tasks, int32_t max_size,
                                 int32_t *part_lengths, uint32_t *task_tokens) {
  if (!generator || !cycle_builder || !cycle_tops || !part_lengths ||
      !task_tokens) {
    return -6;
  }
  const multiply_gates::OpenQasmGenerator &gen =
      *static_cast<multiply_gates::OpenQasmGenerator *>(generator);
  multiply_gates::CycleBuilder &builder =
      static_cast<CycleBuilderWrapper *>(cycle_builder)->builder;
  const multiply_gates::CycleTops &tops = *static_cast<CycleTops *>(cycle_tops);
  return multiply_gates::gen_multiply_gates_tasks(
      gen, builder, tops, multiply_gates::SplitType(split), num_tasks, max_size,
      part_lengths, task_tokens);
}

} // extern "C"
