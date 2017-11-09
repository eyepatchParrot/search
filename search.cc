#include "oracle.h"
#include "interpolate.h"
#include "benchmark.h"
#include "bin.h"
#include "lin.h"
#include "util.h"

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

std::vector<Key> randInts(int n, long seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
  std::vector<Key> v(n);
  for (auto& x : v) x = dist(rng);
  return v;
}

// ./x datasetSz seed nThds [benchmarks ...]
int main(int argc, char *argv[]) {
  if (3 >= argc) return -1;

  using std::istream_iterator;

  int argi = 1;
  int datasetSz = atoi(argv[argi++]);
  long seed = atol(argv[argi++]);
  int nThreads = atoi(argv[argi++]);
  std::vector<std::string> benchmarks;
  while (argi < argc) benchmarks.push_back(argv[argi++]);

  std::cerr << datasetSz << ' ' << seed << '\n';

  int nGets = SUBSET_SIZE < 1 ? datasetSz : SUBSET_SIZE;

  // permute the items
  //std::vector<Key> input;
  auto input = randInts(datasetSz, seed);
  std::vector<int> testIndexes(nGets);
  {
    //auto rndInput = randInts(datasetSz, seed);
    // It should be faster to do three sorts than a sort and a search all
    // even with the copies
    std::vector<std::tuple<Key, int, int>> input_ix(input.size());
    for (int i = 0; i < input.size(); i++)
      input_ix[i] = std::make_tuple(input[i], i, 0);
    std::sort(input_ix.begin(), input_ix.end());
    assert(std::is_sorted(input_ix.begin(), input_ix.end()));

    for (int i = 0; i < input.size(); i++)
      std::get<2>(input_ix[i]) = i;
    for (int i = 0; i < input.size(); i++) input[i] = std::get<0>(input_ix[i]);

    auto cmp = [](auto _i_, auto _j_) { return std::get<1>(_i_) < std::get<1>(_j_); };
    std::sort(input_ix.begin(), input_ix.end(), cmp);
    assert(std::is_sorted(input_ix.begin(), input_ix.end(), cmp));

    for (int i = 0; i < nGets; i++) testIndexes[i] = std::get<2>(input_ix[i]);
  }
  assert(std::is_sorted(input.begin(), input.end()));

  std::unordered_map<std::string, Benchmark*> benchmarkFns{
    { "interpolation-naive", benchmark<InterpolationNaive>},
      { "interpolation-recurse", benchmark<InterpolationRecurse>},
      { "binary-naive", benchmark<BinaryNaive> },
      { "binary-size", benchmark<BinarySize> },
      { "binary-linear", benchmark<BinaryLinear> },
      { "interpolation-linear-fp", benchmark<InterpolationLinearFp> },
      { "isIDiv", benchmark<InterpolationIDiv> },
      { "isLin_2", benchmark<InterpolationLin_2> },
      { "isSub", benchmark<InterpolationSub> },
      { "oracle", benchmark<Oracle> },
      { "interpolation-linear", benchmark<InterpolationLinear>}
  };

  std::vector<TestStats> tests;
  for (auto& benchmark : benchmarks) {
    TestStats ts = benchmarkFns[benchmark](input, testIndexes, nThreads); 
    if (!ts.ok)
      std::cerr << "mess up " << seed << ' ' << benchmark << '\n';
    ts.name = benchmark;
    tests.push_back(ts);
  }

  // print report
  for (int i = 0; i < benchmarks.size(); i++)
    std::cout << (0 != i ? "," : "") << benchmarks[i];
  std::cout << '\n';

  constexpr int nSamples = N_SAMPLES < 0 ? N_RUNS : std::min(N_RUNS, N_SAMPLES);
#ifndef NSORT
  for (auto& t : tests)
    std::partial_sort(t.ns.begin(), t.ns.begin() + nSamples, t.ns.end());
#endif
  for (int runIx = 0; runIx < nSamples; runIx++, std::cout << '\n')
    for (int testIx = 0; testIx < tests.size(); testIx++)
      printf("%s%.3f", 0 != testIx ? "," : "", tests[testIx].ns[runIx]);
}
