#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "omp.h"
#include <string>
#include <vector>
#include "util.h"
#include <numeric>
#include <iostream>

struct TestStats {
  std::string name;
  std::vector<double> ns;
  std::vector<double> thdAvg;
  bool ok;
};

using Benchmark = TestStats(const std::vector<Key>&, const std::vector<int>&);

template <class S>
TestStats benchmark(
    const std::vector<Key>& input,
    const std::vector<int>& indexes) {
  S search(input, indexes);
  constexpr int nRuns = N_RUNS;
  TestStats ts;

  // vals is shuffled, so can't use it. Maybe shuffle indices and use that
  // next time
  std::vector<Key> vals(indexes.size());
  for (int i = 0; i < vals.size(); i++)
    vals[i] = input[indexes[i]];
  
  // get verification info
  auto expSum = 0UL;
  for (auto j=0;j<nRuns;j++) for (auto i : indexes) expSum += vals[i];
  // https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
  int nThreads = 1;
  ts.thdAvg.resize(nThreads);
#pragma omp parallel num_threads(nThreads)
  {
    std::vector<double> ns(nRuns);
    auto valSum = 0UL;
    for (int runIx = 0; runIx < nRuns; runIx++) {
      double t0 = omp_get_wtime();

      for (int i = 0; i < vals.size(); i++) {
        //ERR("%d,%d,", i,indexes[i]);
        auto val = search(vals[i]);
        valSum += val;
        assert(val == vals[i]);
      }
      double t1 = omp_get_wtime();
      ns[runIx] = 1e9*(t1-t0);
    }
#pragma omp barrier
    ts.thdAvg[omp_get_thread_num()] = std::accumulate(ns.begin(), ns.end(), 0.0) / ns.size();
#pragma omp single
    {
      ts.ns = std::move(ns);
      ts.ok = valSum == expSum;
    }
  }
  return ts;
}

#endif
