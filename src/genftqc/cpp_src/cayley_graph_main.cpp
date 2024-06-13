#include "cayley_graph.h"
#include "cayley_graph_wrapper.h"
#include <H5Cpp.h>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

int main_gmgt() {
  multiply_gates::CayleyGraph graph;
  std::string filename =
      std::getenv("HOME") +
      std::string("/d/work/qc4ml/nosync/genftqc_data/multiply_gates/"
                  "multiply_gates_cayley_graph.h5py");
  if (!graph.read(filename, std::cerr)) {
    return 1;
  }
  std::cout << "Read cayley graph" << std::endl;
  graph.fill_downstream_ids();
  auto cycle_tops = graph.find_cycle_tops();
  multiply_gates::OpenQasmGenerator generator{graph.id_to_name};
  std::mt19937_64 rng(0);
  multiply_gates::CycleBuilder builder(graph, rng);
  // builder.from_random_cycle_top(cycle_tops);
  // std::vector<uint8_t> cycle_gates;
  // std::vector<uint32_t> cycle_vertices;
  // builder.get_cycle_gates(cycle_gates);
  // builder.get_cycle_vertices(cycle_vertices);
  // std::cout << "Cycle gates: ";
  // for (auto gate : cycle_gates) {
  //   std::cout << static_cast<int>(gate) << " ";
  // }
  // std::cout << "\nOpenQASM:\n";
  // generator.cycle_to_openqasm(cycle_gates, std::cout);
  // std::cout << "\nCycle vertices: ";
  // for (auto vertex : cycle_vertices) {
  //   std::cout << vertex << " ";
  // }
  // std::cout << std::endl;
  int32_t num_tasks = 1024;
  int32_t max_size = 43;
  std::vector<int32_t> part_lengths(2 * num_tasks);
  std::vector<multiply_gates::TokenId> task_tokens(num_tasks * max_size);
  int32_t ret = multiply_gates::gen_multiply_gates_tasks(
      generator, builder, cycle_tops, multiply_gates::SplitType::kTrain,
      num_tasks, max_size, part_lengths.data(), task_tokens.data());
  std::cout << "Res: " << ret << std::endl;
  return 0;
}

int main_fail() {
  // This main function uses ctypes wrappers to simulate the Python code below.
  /*
  import sys
  if "." not in sys.path:
      sys.path.append(".")
  import multiply_gates
  import os

  cayley_graph = multiply_gates.CayleyGraph()
  filename = os.path.expanduser(
      "~/d/work/qc4ml/nosync/genftqc_data/multiply_gates/"
      "multiply_gates_cayley_graph.h5py")
  cayley_graph.read(filename)
  cayley_graph.fill_downstream_ids()
  cycle_tops = multiply_gates.CycleTops(cayley_graph)
  builder = multiply_gates.CycleBuilder(cayley_graph, 1)
  part_lengths, task_tokens = cayley_graph.gen_multiply_gates_tasks(
      builder, cycle_tops)
  */
  std::unique_ptr<multiply_gates::CayleyGraph> cayley_graph(
      static_cast<multiply_gates::CayleyGraph *>(create_cayley_graph()));
  std::string filename =
      std::getenv("HOME") +
      std::string("/d/work/qc4ml/nosync/genftqc_data/multiply_gates/"
                  "multiply_gates_cayley_graph.h5py");
  if (!cayley_graph_read(cayley_graph.get(), filename.c_str())) {
    return 1;
  }
  cayley_graph_fill_downstream_ids(cayley_graph.get());
  std::unique_ptr<multiply_gates::CycleTops> cycle_tops(
      static_cast<multiply_gates::CycleTops *>(
          create_cycle_tops(cayley_graph.get())));
  void *cycle_builder = create_cycle_builder(cayley_graph.get(), 1);
  int32_t num_tasks = 1024;
  int32_t max_size = 43;
  std::vector<int32_t> part_lengths(2 * num_tasks);
  std::vector<uint32_t> task_tokens(num_tasks * max_size);
  std::unique_ptr<multiply_gates::OpenQasmGenerator> gen(
      static_cast<multiply_gates::OpenQasmGenerator *>(
          create_openqasm_generator(cayley_graph.get())));
  int32_t ret = gen_multiply_gates_tasks(
      gen.get(), cycle_builder, cycle_tops.get(),
      int32_t(multiply_gates::SplitType::kTrain), num_tasks, max_size,
      part_lengths.data(), task_tokens.data());
  destroy_cycle_builder(cycle_builder);
  return ret;
}

int main() { return main_fail(); }
