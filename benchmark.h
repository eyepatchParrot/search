#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "omp.h"
#include <string>
#include <vector>
#include "util.h"
#include <numeric>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>

struct TestStats {
  std::string name;
  std::vector<double> ns;
  std::vector<double> thdAvg;
  bool ok;
};

using Benchmark = TestStats(const std::vector<Key>&, const std::vector<int>&, const int);

double mode(std::vector<double> X) {
  std::unordered_map<int, int> cnt;
  //for (auto x : X) cnt[round(x*100.0)]++;
  for (auto x : X) {
    int ix = round(x*100.0);
    if (0 == cnt.count(ix)) cnt[ix] = 0;
    cnt[ix]++;
  }
  auto it = std::max_element(cnt.begin(), cnt.end(),
      [](std::pair<int, int> l, std::pair<int, int> r) {
      return l.second < r.second;
      });
  std::cerr << it->second << '\n';
  return it->first / 100.0;
}

template <class S>
TestStats benchmark(
    const std::vector<Key>& input,
    const std::vector<int>& indexes,
    const int nThreads = 1) {
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
  ts.thdAvg.resize(nThreads);
#pragma omp parallel num_threads(nThreads)
  {
    std::vector<double> ns(nRuns, 0);
    auto valSum = 0UL;
    for (int runIx = 0; runIx < nRuns; runIx++) {
      auto t0 = std::chrono::steady_clock::now();
      //double t0 = omp_get_wtime();

      for (int i = 0; i < vals.size(); i++) {
        auto val = search(vals[i]);
        valSum += val;
        assert(val == vals[i]);
      }
      auto t1 = std::chrono::steady_clock::now();
      //double t1 = omp_get_wtime();
      //ns[runIx] = (t1-t0) * 1e9/ (double)vals.size();
      ns[runIx] = (std::chrono::nanoseconds(t1-t0)).count() / (double)vals.size();
    }
#pragma omp barrier
    ts.thdAvg[omp_get_thread_num()] = std::accumulate(ns.begin(), ns.end(), 0.0) / ns.size();
#pragma omp single
    {
      ts.ns = std::move(ns);
      ts.ok = valSum == expSum;
    }
  }
  //std::cerr << "Mode " << mode(ts.ns) << " min " << *std::min_element(ts.ns.begin(), ts.ns.end()) << " mean " << ts.thdAvg[0] << '\n';
  return ts;
}

#endif
