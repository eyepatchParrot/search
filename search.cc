#include "oracle.h"
#include "interpolate.h"
#include "benchmark.h"
#include "bin.h"
#include "lin.h"
#include "util.h"
#include "fib.h"

#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>
#include <x86intrin.h>
#include <unordered_map>

std::vector<Key> randNums(int n, long seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
  std::vector<Key> v(n);
  for (auto& x : v) x = dist(rng);
  return v;
}

struct Input {
  private:
    struct KeyRandomSortedIndex {
      Key key;
      int random_index;
      int sorted_index;

      KeyRandomSortedIndex(Key k) { key = k; }
    };
  public:

  std::vector<Key> keys;
  std::vector<int> indexes;

  Input(std::vector<Key> nums, int n_gets) : keys(nums.size()), indexes(n_gets) {
#ifdef DUMP_INPUT
    return;
#endif

    std::vector<KeyRandomSortedIndex> key_indexes(nums.begin(), nums.end());
    for (int i = 0; i < key_indexes.size(); i++)
      key_indexes[i].random_index = i;
    std::sort(key_indexes.begin(), key_indexes.end(), [](auto l, auto r) {
        return l.key < r.key; });
    for (int i = 0; i < key_indexes.size(); i++)
      key_indexes[i].sorted_index = i;

    // save the sorted input while it's here
    for (auto i = 0; i < nums.size(); i++) keys[i] = key_indexes[i].key;

    std::sort(key_indexes.begin(), key_indexes.end(), [](auto l, auto r) {
        return l.random_index < r.random_index; });
    for (int i = 0; i < n_gets; i++) indexes[i] = key_indexes[i].sorted_index;
  }
};

// ./x datasetSz seed nThds [benchmarks ...]
int main(int argc, char *argv[]) {
  if (3 >= argc) return -1;

  using std::istream_iterator;

  int argi = 1;
  int datasetSz = atoi(argv[argi++]);
  long seed = atol(argv[argi++]);
  int n_threads = atoi(argv[argi++]);
  std::vector<std::string> benchmarks;
  while (argi < argc) benchmarks.push_back(argv[argi++]);

  std::cerr << datasetSz << ' ' << seed << '\n';

  int nGets = SUBSET_SIZE < 1 ? datasetSz : SUBSET_SIZE;

  Input input(randNums(datasetSz, seed), nGets);

#ifdef DUMP_INPUT
  for (auto x : input) std::cout << x << '\n';
  return 0;
#endif
  std::unordered_map<std::string, Benchmark*> benchmarkFns{
    {"b0", benchmark<B0>},
      {"b1", benchmark<B1>},
      {"i", benchmark<InterpolationNaive>},
      {"i-precompute", benchmark<InterpolationPrecompute>},
      {"i-lut", benchmark<InterpolationRecurseLut>},
      {"i-seq-fp", benchmark<InterpolationLinearFp> },
      {"i-seq", benchmark<InterpolationLinear>},
      {"i-guard", benchmark<InterpolationRecurseGuard>},
      {"i-seq-simd", benchmark<i_simd>},

      {"b-lr", benchmark<BinaryLR<>>},
      {"b-lr-cond", benchmark<BinaryLRCond>},
      {"b-lr-noeq", benchmark<BinaryLRNoEq>},
      {"b-lr-for", benchmark<BinaryLRFor>},
      {"b-lr-noeq-for", benchmark<BinaryLRNoEqFor>},
      {"b-lr-over", benchmark<BinaryLROver>},
      {"b-lr-lin", benchmark<BinaryLRLin>},

      {"b-sz", benchmark<BinarySz<>>},
      {"b-sz-cond", benchmark<BinarySzCond>},
      {"b-sz-noeq", benchmark<BinarySzNoEq>},
      {"b-sz-for", benchmark<BinarySzFor>},
      {"b-sz-noeq-for", benchmark<BinarySzNoEqFor>},
      {"b-sz-pow", benchmark<BinarySzPow>},
      {"b-sz-lin", benchmark<BinarySzLin>},
      {"b2", benchmark<B2>},
      {"fib", benchmark<Fib>},

      {"oracle", benchmark<Oracle> }
  };

  std::vector<TestStats> tests;
  for (auto& benchmark : benchmarks) {
    TestStats ts = benchmarkFns[benchmark](input.keys, input.indexes, n_threads); 
    if (!ts.ok)
      std::cerr << "mess up " << seed << ' ' << benchmark << '\n';
    ts.name = benchmark;
    tests.push_back(ts);
  }

  // print report
  for (int i = 0; i < benchmarks.size(); i++)
    std::cout << (0 != i ? "," : "") << benchmarks[i];
  std::cout << '\n';

  constexpr int n_samples = N_SAMPLES < 0 ? N_RUNS : std::min(N_RUNS, N_SAMPLES);
#ifndef NSORT
  for (auto& t : tests)
    std::partial_sort(t.ns.begin(), t.ns.begin() + n_samples, t.ns.end());
#endif
  for (int runIx = 0; runIx < n_samples; runIx++, std::cout << '\n')
    for (int testIx = 0; testIx < tests.size(); testIx++)
      printf("%s%.3f", 0 != testIx ? "," : "", tests[testIx].ns[runIx]);
}
