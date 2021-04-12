#include "FM.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    std::cout << "Usage: ./partioner [path to *.hgr file] [optional max difference]" << std::endl;
    return 1;
  }
  std::filesystem::path input_path = argv[1],
                        output_path = std::filesystem::current_path() /
                                      input_path.filename().replace_extension(".part.2");
  unsigned max_diff = 1;
  if (argc == 3)
    if (unsigned n = std::stoul(argv[2]))
      max_diff = n;
  // Load hypergraph
  std::ifstream input(input_path.c_str());
  if (!input.is_open()) {
    std::cout << "Can not open " << input_path << std::endl;
    return 1;
  }
  unsigned num_nets, num_cells;
  input >> num_nets >> num_cells;
  Hypergraph hypergraph(num_cells);
  std::string cur_str;
  std::getline(input, cur_str);
  for (unsigned i = 0; i < num_nets; ++i) {
    std::getline(input, cur_str);
    std::istringstream cur_str_buf(cur_str);
    Net cur_net{std::istream_iterator<unsigned>(cur_str_buf), std::istream_iterator<unsigned>()};
    for (unsigned &cell : cur_net)
      --cell;
    hypergraph.addNet(cur_net);
  }
  input.close();
  // Do partioning
  Partitionment partitionment(num_cells);
  unsigned cost = partitionment.getPartitionCost(hypergraph), pass_num = 0;
  std::cout << "Initial cut cost: " << cost << std::endl;
  while (auto cost_reduction = FM_pass(hypergraph, partitionment, max_diff)) {
    cost -= cost_reduction;
    ++pass_num;
    std::cout << "Cut cost after pass #" << std::setw(2) << pass_num << ": " << std::setw(5) << cost
              << std::endl;
    assert(cost == partitionment.getPartitionCost(hypergraph));
  }
  std::cout << "Final cut cost: " << partitionment.getPartitionCost(hypergraph) << std::endl;
  // Save results
  std::ofstream output(output_path.c_str());
  if (!output.is_open()) {
    std::cout << "Can not open " << output_path << std::endl;
    return 1;
  }
  for (unsigned i = 0; i < num_cells; ++i)
    output << +partitionment[i] << std::endl;
  output.close();
  return 0;
}
