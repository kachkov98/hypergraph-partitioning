#include "FM.hpp"
#include <algorithm>
#include <iostream>
#include <set>

void Hypergraph::addNet(const Net &net) {
  assert(net.size() >= 2);
  unsigned net_index = nets_.size();
  for (auto cell_index : net) {
    assert(cell_index < cells_.size());
    cells_[cell_index].push_back(net_index);
  }
  nets_.push_back(net);
}

unsigned Hypergraph::getMaxAdjacency() const {
  unsigned max_adjacency = 0;
  for (const auto &cell : getCells())
    max_adjacency = std::max(max_adjacency, (unsigned)cell.size());
  return max_adjacency;
}

Partitionment::Partitionment(unsigned num_cells) {
  partitionment_.reserve(num_cells);
  for (unsigned i = 0; i < num_cells; ++i)
    partitionment_.push_back(i % 2);
}

unsigned Partitionment::getPartitionCost(const Hypergraph &hypergraph) const {
  assert(partitionment_.size() == hypergraph.getCells().size());
  unsigned cost = 0;
  for (const auto &net : hypergraph.getNets()) {
    std::set<PartId> net_parts;
    for (unsigned cell : net)
      net_parts.insert(partitionment_[cell]);
    cost += net_parts.size() - 1;
  }
  return cost;
}

size_t GainBuckets::addCell(int gain, unsigned cell) {
  ++num_elems_;
  if (gain > best_gain_)
    best_gain_ = gain;
  (*this)[gain].push_back(cell);
  return (*this)[gain].size() - 1;
}

unsigned GainBuckets::delCell(int gain, size_t iter) {
  --num_elems_;
  unsigned cell = UINT_MAX;
  if (iter != (*this)[gain].size() - 1)
    cell = (*this)[gain][iter] = (*this)[gain].back();
  (*this)[gain].pop_back();
  if (best_gain_ == gain)
    while ((*this)[best_gain_].empty()) {
      if (best_gain_ == -(int)max_gain_) {
        best_gain_ = INT_MIN;
        break;
      }
      --best_gain_;
    }
  return cell;
}

GainBuckets::GainedCell GainBuckets::getMaxGainedCell() const {
  assert(best_gain_ != INT_MIN);
  return GainedCell{(*this)[best_gain_].back(), best_gain_};
}

GainContainer::GainContainer(const Hypergraph &hypergraph, const Partitionment &partitionment) {
  // create PartitionGains for both parts
  auto max_gain = hypergraph.getMaxAdjacency();
  for (unsigned i = 0; i < 2; ++i)
    buckets_.emplace_back(max_gain);
  // calculate number of cells in each partition for all nets
  num_cells_in_partitions_.reserve(hypergraph.getNets().size());
  for (const auto &net : hypergraph.getNets()) {
    CellsInPartitions num_cells{0, 0};
    for (unsigned cell : net)
      ++num_cells[partitionment[cell]];
    num_cells_in_partitions_.push_back(num_cells);
  }
  // calculate gains for all cells
  cell_states_.reserve(hypergraph.getCells().size());
  for (unsigned cell = 0; cell < hypergraph.getCells().size(); ++cell) {
    auto cell_part = partitionment[cell];
    int gain = 0;
    for (unsigned net : hypergraph.getCells()[cell])
      if (num_cells_in_partitions_[net][cell_part] == 1)
        ++gain;
      else if (num_cells_in_partitions_[net][OppositePart(cell_part)] == 0)
        --gain;
    cell_states_.push_back(CellState{buckets_[cell_part].addCell(gain, cell), gain, false});
  }
};

GainContainer::GainedMove GainContainer::getBestMove() const {
  PartId biggest_part = buckets_[1].getNumElems() > buckets_[0].getNumElems() ? 1 : 0;
  auto [cell, gain] = buckets_[biggest_part].getMaxGainedCell();
  return GainedMove{{cell, biggest_part, OppositePart(biggest_part)}, gain};
}

void GainContainer::update(const Hypergraph &hypergraph, const Partitionment &partitionment,
                           Move move) {
  auto &cell_state = cell_states_[move.cell];
  if (auto cell = buckets_[move.src].delCell(cell_state.gain, cell_state.iter); cell != UINT_MAX)
    cell_states_[cell].iter = cell_state.iter;
  cell_state.locked = true;
  for (unsigned net : hypergraph.getCells()[move.cell]) {
    // updates before move
    if (num_cells_in_partitions_[net][move.dst] == 0) {
      for (unsigned cell : hypergraph.getNets()[net])
        if (cell != move.cell)
          updateGain(cell, move.src, +1);
    }
    if (num_cells_in_partitions_[net][move.dst] == 1) {
      for (unsigned cell : hypergraph.getNets()[net])
        if (cell != move.cell && partitionment[cell] == move.dst)
          updateGain(cell, move.dst, -1);
    }
    // move
    --num_cells_in_partitions_[net][move.src];
    ++num_cells_in_partitions_[net][move.dst];
    // updates after move
    if (num_cells_in_partitions_[net][move.src] == 0) {
      for (unsigned cell : hypergraph.getNets()[net])
        if (cell != move.cell)
          updateGain(cell, move.dst, -1);
    }
    if (num_cells_in_partitions_[net][move.src] == 1) {
      for (unsigned cell : hypergraph.getNets()[net])
        if (cell != move.cell && partitionment[cell] == move.src)
          updateGain(cell, move.src, +1);
    }
  }
}

void GainContainer::updateGain(unsigned cell, PartId part, int diff) {
  auto &cell_state = cell_states_[cell];
  if (cell_state.locked)
    return;
  if (auto moved_cell = buckets_[part].delCell(cell_state.gain, cell_state.iter);
      moved_cell != UINT_MAX)
    cell_states_[moved_cell].iter = cell_state.iter;
  cell_state.gain += diff;
  cell_state.iter = buckets_[part].addCell(cell_state.gain, cell);
}

unsigned FM_pass(const Hypergraph &hypergraph, Partitionment &partitionment) {
  unsigned best_cost_reduction = 0;
  Partitionment best_partitionment(partitionment);
  GainContainer gain_container(hypergraph, partitionment);
  int cur_cost_reduction = 0;
  for (unsigned i = 0; i < hypergraph.getCells().size(); ++i) {
    auto [move, gain] = gain_container.getBestMove();
    gain_container.update(hypergraph, partitionment, move);
    partitionment.doMove(move);
    cur_cost_reduction += gain;
    if (cur_cost_reduction > (int)best_cost_reduction) {
      best_cost_reduction = (unsigned)cur_cost_reduction;
      best_partitionment = partitionment;
    }
  }
  partitionment = std::move(best_partitionment);
  return best_cost_reduction;
}
