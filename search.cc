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

// ./x fileName [benchmarks ...]
int main(int argc, char *argv[]) {
  using std::istream_iterator;
  constexpr int seed = 42;

  std::cerr << argv[1] << '\n';

  std::ifstream f(argv[1]);
  
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
  Oracle o(input, testIndexes);

  std::vector<TestStats> tests;
  for (int i = 2; i < argc; i++) {
    TestStats ts; 
    std::string s = argv[i];
    if (s == "binary-naive") ts = benchmark<Binary<BS_EQ>>(input,testIndexes);
    else if (s == "binary-size") ts = benchmark<Binary<BS_LIN,1>>(input,testIndexes);
    else if (s == "binary-linear") ts = benchmark<Binary<BS_LIN>>(input,testIndexes);
    else if (s == "interpolation-naive") ts = benchmark<Interpolation<IS_FP,-1,false,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-recurse") ts = benchmark<Interpolation<IS_FP,-1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-linear-fp") ts = benchmark<Interpolation<IS_FP,1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-linear") ts = benchmark<Interpolation<>>(input, testIndexes);
    else if (s == "isRecurse") ts = benchmark<Interpolation<IS_LIN,-1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "isIDiv") ts = benchmark<Interpolation<IS_IDIV>>(input, testIndexes);
    else if (s == "isLin_1_slow") ts = benchmark<Interpolation<IS_LIN,1,false>>(input, testIndexes);
    else if (s == "isLin_2") ts = benchmark<Interpolation<IS_LIN,2>>(input, testIndexes);
    else if (s == "isSub") ts = benchmark<Interpolation<IS_LIN,1,true>>(input, testIndexes);
    else if (s == "oracle") ts = benchmark(input, testIndexes, o);
    if (!ts.ok)
      std::cerr << "mess up " << argv[1] << ' ' << s << '\n';
    tests.push_back(ts);
  }

  // Set-up the headers and organize data
  std::vector<std::vector<size_t>> testsIxs(tests.size(),
      std::vector<size_t>(N_RUNS));
  for (int i = 0; i < tests.size(); i++) {
    auto& t = tests[i];
    assert(t.ns.size() == N_RUNS);
    printf("%s%s", 0 != i ? "," : "", t.name.c_str());
    std::iota(testsIxs[i].begin(), testsIxs[i].end(), 0);
#ifndef NSORT
    std::sort(testsIxs[i].begin(), testsIxs[i].end(), [t](size_t a, size_t b) {
        return t.ns[a] < t.ns[b]; });
#endif
  }
  for (int runIx = 0; runIx < N_RUNS && (N_SAMPLES < 0 ? true : runIx < N_SAMPLES); runIx++) {
    printf("\n");
    for (int testIx = 0; testIx < tests.size(); testIx++) {
      auto& t = tests[testIx];
      auto nsIx = testsIxs[testIx][runIx];
      printf("%s%.3f", 0 != testIx ? "," : "",
          t.ns[nsIx] / (double)testIndexes.size());
    }
  }
  printf("\n");
}
