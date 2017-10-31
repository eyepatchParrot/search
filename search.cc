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

// ./x fileName nThds [benchmarks ...]
int main(int argc, char *argv[]) {
  using std::istream_iterator;
  constexpr int seed = 42;

  if (2 >= argc) return -1;

  std::string fileName = argv[1];
  int nThreads = atoi(argv[2]);
  std::vector<std::string> benchmarks;
  for (int i = 3; i < argc; i++) benchmarks.push_back(argv[i]);

  std::cerr << fileName << '\n';

  std::ifstream f(fileName);
  
  int nNums;
  f >> nNums;
  assert(nNums != 0);
  auto input = std::vector<Key>(istream_iterator<Key>(f), istream_iterator<Key>());
  int nGets = SUBSET_SIZE < 1 ? input.size() : SUBSET_SIZE;

  // permute the items
  std::vector<int> testIndexes(nGets);
  {
    std::vector<int> allIndexes(input.size());
    std::iota(allIndexes.begin(), allIndexes.end(), 0);
    std::shuffle(allIndexes.begin(), allIndexes.end(), std::mt19937{seed});
    std::copy_n(allIndexes.begin(), nGets, testIndexes.begin());
  }

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
      std::cerr << "mess up " << fileName << ' ' << benchmark << '\n';
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
