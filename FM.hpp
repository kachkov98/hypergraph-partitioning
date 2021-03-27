#ifndef FM_HPP
#define FM_HPP

#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <list>
#include <vector>

using Net = std::vector<unsigned>;  // list of all cells this net connect
using Cell = std::vector<unsigned>; // list of all nets that are connected to this cell

class Hypergraph {
public:
  Hypergraph(unsigned num_cells) : cells_(num_cells){};
  void addNet(const Net &net);
  const std::vector<Cell> &getCells() const { return cells_; }
  const std::vector<Net> &getNets() const { return nets_; }
  unsigned getMaxAdjacency() const;

private:
  std::vector<Cell> cells_;
  std::vector<Net> nets_;
};

using PartId = uint8_t;

struct Move {
  unsigned cell;
  PartId src, dst;
};

class Partitionment {
public:
  Partitionment(unsigned num_cells);
  PartId operator[](unsigned cell) const { return partitionment_[cell]; }
  void doMove(const Move &move) {
    assert(partitionment_[move.cell] == move.src);
    partitionment_[move.cell] = move.dst;
  }
  unsigned getPartitionCost(const Hypergraph &hypergraph) const;

private:
  std::vector<PartId> partitionment_;
};

static PartId OppositePart(PartId id) { return !id; }

class GainBuckets {
public:
  using Bucket = std::list<unsigned>;
  using GainedCell = std::pair<unsigned, int>;
  GainBuckets(unsigned max_gain)
      : container_(max_gain * 2 + 1), max_gain_(max_gain), num_elems_(0), best_gain_(INT_MIN) {}

  Bucket::const_iterator addCell(int gain, unsigned cell);
  void delCell(int gain, Bucket::const_iterator iter);
  Bucket::const_iterator moveCell(int old_gain, Bucket::const_iterator old_iter, int new_gain);
  GainedCell getMaxGainedCell() const;
  unsigned getNumElems() const { return num_elems_; }

private:
  std::vector<Bucket> container_;
  unsigned max_gain_, num_elems_;
  int best_gain_;

  const Bucket &operator[](int index) const {
    assert(-(int)max_gain_ <= index && index <= (int)max_gain_);
    return container_[index + max_gain_];
  }

  Bucket &operator[](int index) {
    assert(-(int)max_gain_ <= index && index <= (int)max_gain_);
    return container_[index + max_gain_];
  }
};

class GainContainer {
public:
  using GainedMove = std::pair<Move, int>; // move and its gain
  GainContainer(const Hypergraph &hypergraph, const Partitionment &partitionment);
  GainedMove getBestMove() const;
  void update(const Hypergraph &hypergraph, const Partitionment &partitionment, Move move);

private:
  struct CellState {
    GainBuckets::Bucket::const_iterator iter;
    int gain;
    bool locked;
  };
  std::vector<CellState> cell_states_;
  std::vector<GainBuckets> buckets_; // buckets for each partition
  using CellsInPartitions = std::array<unsigned, 2>;
  std::vector<CellsInPartitions>
      num_cells_in_partitions_; // number of cells in each partition for all nets
  void updateGain(unsigned cell, PartId part, int diff);
};

unsigned FM_pass(const Hypergraph &hypergraph, Partitionment &partitionment);

#endif
