#ifndef CAYLEY_GRAPH_WRAPPER_H
#define CAYLEY_GRAPH_WRAPPER_H

#include <cstdint>

extern "C" {

// Create and destroy CayleyGraph object
void *create_cayley_graph();
void destroy_cayley_graph(void *cayley_graph);

// Read method for CayleyGraph
bool cayley_graph_read(void *cayley_graph, const char *filename);

// Fill downstream ids
void cayley_graph_fill_downstream_ids(void *cayley_graph);

// Find cycle tops
void *create_cycle_tops(void *cayley_graph);
void destroy_cycle_tops(void *cycle_tops);

// Get cycle top vertices. Caller should first call
// `cycle_tops_get_num_vertices` to determine the size of the buffer needed to
// store the vertices. `cycle_tops_get_vertices` fills the buffer of that size.
int32_t cycle_tops_get_num_vertices(void *cycle_tops);
void cycle_tops_get_vertices(void *cycle_tops, uint32_t *vertices);

int32_t cycle_tops_get_num_edges(void *cycle_tops);
void cycle_tops_get_edges(void *cycle_tops, uint32_t *edge_u, uint32_t *edge_v,
                          uint8_t *edge_g);

// Create and destroy CycleBuilderWrapper object
void *create_cycle_builder(void *cayley_graph, uint64_t seed);
void destroy_cycle_builder(void *cycle_builder);

// Get cycle gates and vertices. `cycle_gates` or `cycle_vertices` should point
// to a buffer of size at least `max_size`. The actual size filled is returned.
// Negative return values indicate an error.
int32_t cycle_builder_get_max_cycle_size(void *cycle_builder);
int32_t cycle_builder_get_cycle_gates(void *cycle_builder, bool invert_last,
                                      uint8_t *cycle_gates, int32_t max_size);
int32_t cycle_builder_get_cycle_vertices(void *cycle_builder,
                                         uint32_t *cycle_vertices,
                                         int32_t max_size);

// Start from vertex or edge
void cycle_builder_from_vertex(void *cycle_builder, uint32_t vertex_id);
void cycle_builder_from_edge(void *cycle_builder, uint32_t u, uint32_t v,
                             uint8_t g);

// Start from random cycle top
void cycle_builder_from_random_cycle_top(void *cycle_builder, void *cycle_tops);

// multiply_gates::OpenQasmGenerator generator{graph.id_to_name};
// generator.cycle_to_openqasm(cycle_gates, std::cout);
void *create_openqasm_generator(void *cayley_graph);
void destroy_openqasm_generator(void *generator);

// Populate `out` with the OpenQASM representation of the cycle gates.
// `cycle_gates` should point to a buffer of size `out_size`.
//
// Returns the number of characters written to `out` or -1 if out_size was too
// small to contain the OpenQASM representation (in which case the output is
// truncated to `out_size` characters).
int32_t cycle_to_openqasm(void *generator, const uint8_t *cycle_gates,
                          int32_t cycle_gates_size, char *out,
                          int32_t out_size);

// If out_size is sufficient for the output, returns the number of tokens
// written. Otherwise, returns -1 and writes the first out_size tokens to out.
int32_t cycle_to_tokens(void *generator, const uint8_t *cycle_gates,
                        int32_t cycle_gates_size, uint32_t *out,
                        int32_t out_size);
int32_t detokenize(void *generator, const uint32_t *tokens, int32_t tokens_size,
                   char *out, int32_t out_size);

int32_t gen_multiply_gates_tasks(void *generator, void *cycle_builder,
                                 void *cycle_tops, int32_t split,
                                 int32_t num_tasks, int32_t max_size,
                                 int32_t *part_lengths, uint32_t *task_tokens);
}

#endif // CAYLEY_GRAPH_WRAPPER_H
