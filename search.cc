#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <immintrin.h>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <time.h>
#include <tuple>
#include <vector>
#include <x86intrin.h>
#include <limits>

#include "util.h"
#include "div.h"

typedef uint64_t ll;

DivLut<8> dL;

struct RunStats {
  int sumSteps1, sumSteps2;
  int maxSteps1, maxSteps2;
  //  ll maxV;
};

struct TestStats {
  std::string name;
  std::vector<RunStats> runStats;
  std::vector<std::tuple<ll, int> > cyclesByIx;

  void datum(ll dt) {
      cyclesByIx.emplace_back(dt, cyclesByIx.size());
  }
  void datum(ll st, ll et) { datum(et-st); }
  void data(TestStats& ts) {
      for (auto t_i : ts.cyclesByIx)
          datum(std::get<0>(t_i));
  }
};

template <typename S, int (f) (const ll, S& s)>
TestStats benchmark(
    std::string&& name,
    const std::vector<ll>& vals,
    const std::vector<int>& indexes,
    S& s) {
  constexpr int nRuns = N_RUNS;
  const auto N_NUMS = vals.size();
  auto ts = TestStats{name};
  auto expSum = 0, ixSum = 0;
  for (auto j=0;j<nRuns;j++) for (auto ix : indexes) expSum += ix;
  for (int runIx = 0; runIx < nRuns; runIx++) {
    ll st = __rdtsc();

    for (int i = 0; i < N_NUMS; i++) {
      ERR("%d,%d,", i,indexes[i]);
      int ix = f(vals[i], s);
      ixSum += ix;
      assert(ix == indexes[i]);
    }

    ll et = __rdtsc();
    ts.datum(st, et);
  }
  if (ixSum != expSum) // check even when profiling
    fprintf(stderr, "mess up\n");
  return ts;
}


struct BinStruct {
  BinStruct(const std::vector<ll>& _a) : v(_a.size() + 32,0) {
    std::copy(_a.begin(), _a.end(), v.begin() + 16);
    for (int i = 0; i < 16; i++) v[i] = std::numeric_limits<ll>::min();
    for (int i = _a.size() + 16; i < v.size(); i++) v[i] = std::numeric_limits<int64_t>::max();
  }
  auto szA() { return v.size() - 32; };
  const ll* a() { return &v[16]; };
  private:
  std::vector<ll> v;
};

struct IntStruct {
  IntStruct(const std::vector<ll>& _a) : v(_a.size() + 64,0) {
    std::copy(_a.begin(), _a.end(), v.begin() + 32);
    for (int i = 0; i < 32; i++) v[i] = std::numeric_limits<ll>::min();
    for (int i = _a.size() + 32; i < v.size(); i++) v[i] = std::numeric_limits<int64_t>::max();
    lgScale = lg(_a.size() - 1);
    ll tR = _a.size() - 1;
    ll rSc = __builtin_clz(r);
    r2 = tR << rSc;
    lg_d2 = lgl(_a.back() - _a.front()) + rSc - 64;
    lg_d = lgl(_a.back() - _a.front()) - lgScale;
    yLS = _a.front() >> lgScale;
    ll d = (_a.back() - _a.front()) >> lgScale;
    dL.one_d(d, p, lg_q);
  }

  unsigned lgScale;
  unsigned lg_d;
  unsigned lg_d2;
  ll r;
  ll r2;
  ll yLS;
  ll p; int lg_q;
  std::vector<ll> v;
  auto szA() { return v.size() - 64; };
  const ll* a() { return &v[32]; };
};

struct OracleStruct {
  // note that to ensure fast code, offset is always to the left, which means that offset searches wont have the full offset
  OracleStruct(const std::vector<ll>& _a, const std::vector<int>& _i, const int offset = 0, const bool rnd = false) : a(_a), j(0) {
    // easier to make a single long pipeline instead of a branch
    // There's blocking on the output to be produced, instead of being able to go through the steady state throughput
    // You may end up paying the latency cost multiple times instead of just once with left-deep
    // single relation plan because only one table in FROM clause
    for (int ti : _i) {
      int tti = ti;
      if (rnd && ((rand() % 2) == 0)) {
        tti += offset;
      } else {
        tti -= offset;
      }
      unsigned r = a.size() - 1;
      tti = tti < 0 ? 0 : tti;
      tti = tti > r ? r : tti;
      i.push_back(tti);
    }
  }

  const std::vector<ll>& a;
  std::vector<int> i;
  int j;
};

int bsPVKEq2(const ll x, BinStruct& s) {
  const ll* array = s.a();
  const int MIN_EQ_SZ = 2;
  long leftIndex = 0;                                                               
  int n = s.szA();
  int half;
  while ((half = n) > MIN_EQ_SZ) {
    half /= 2;
    n -= half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  while ((half = n) > 1) {
    half /= 2;
    n = array[leftIndex + half] == x ? 0 : n - half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }                                                                                
  assert(array[leftIndex] == x);  
  return leftIndex;
}

using SearchFn = int64_t(const ll*, int64_t, ll);
template < bool reverse=false, int roll=3>
int64_t linSIMD(const ll* arr, const int64_t guessIx, const ll x) {
  auto vecXX = reverse? _mm256_set1_epi64x(x): _mm256_set1_epi64x(x-1);
  auto ptr = arr;
  auto i = guessIx;
	auto misalignment = ((uintptr_t)(ptr+i) & 31)/sizeof(ll);
	for (int j = 0; j < 4*roll; j++)
    if (reverse? (arr[i-j] <= x) : arr[i+j] >= x) return reverse? i-j : i+j;
  i = reverse? (i-4*(roll-1) - misalignment) : i + 4*roll - misalignment;
  // 32-aligned main loop                                                          
  for (;;i = reverse?(i-16) : i+16) {
    assert(i<1000);
    assert(i>-32);
    auto sign = reverse?-1:1;
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i + sign*0));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + sign*4));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + sign*8));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + sign*12));
    auto cmp3 = reverse? _mm256_cmpgt_epi64(vecXX, av3) :  _mm256_cmpgt_epi64(av3, vecXX);
    auto msk3 = _mm256_movemask_epi8(cmp3);
    if (!msk3) continue;
    auto cmp0 = reverse? _mm256_cmpgt_epi64(vecXX, av0) :   _mm256_cmpgt_epi64(av0,vecXX );
    auto cmp1 = reverse? _mm256_cmpgt_epi64(vecXX, av1) :   _mm256_cmpgt_epi64(av1,vecXX );
    auto cmp2 = reverse? _mm256_cmpgt_epi64(vecXX, av2) :   _mm256_cmpgt_epi64(av2,vecXX );
    auto msk0 = _mm256_movemask_epi8(cmp0);
    auto msk1 = _mm256_movemask_epi8(cmp1);
    auto msk2 = _mm256_movemask_epi8(cmp2);
    if (msk0) return reverse? (i + 4 - _lzcnt_u32(msk0) / 8 - 0 * 4) : i + _tzcnt_u32(msk0) / 8 + 0 * 4;
    if (msk1) return reverse? (i + 4 - _lzcnt_u32(msk1) / 8 - 1 * 4) : i + _tzcnt_u32(msk1) / 8 + 1 * 4;
    if (msk2) return reverse? (i + 4 - _lzcnt_u32(msk2) / 8 - 2 * 4) : i + _tzcnt_u32(msk2) / 8 + 2 * 4;
    if (msk3) return reverse? (i + 4 - _lzcnt_u32(msk3) / 8 - 3 * 4) : i + _tzcnt_u32(msk3) / 8 + 3 * 4;
  }
}

template <bool reverse=false,int n=12>
int64_t linUnroll(const ll* a, int64_t m, ll k) {
  for (;;m = (reverse?m-n:m+n)) {
    for (int i = 0; i < n; i++) {
      assert(m+i < 1032); assert((m-i) > -32);
      if (reverse?(a[m-i]<=k):(a[m+i]>=k)) return reverse?(m-i):(m+i);
    }
  }
}

int bsLin(const ll x, BinStruct& s) {
  const ll* array = s.a();
  const int MIN_EQ_SZ = 32;
  const int roll = 12;
  auto leftIndex = 0;                                                               
  auto n = s.szA();
  while (n > MIN_EQ_SZ) {
    auto half = n / 2;
    n -= half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  auto guess = leftIndex + n/2;
  if (array[guess] < x) {
    guess++;
    for (;;guess+=roll) for(int i=0;i<roll;i++)  if (array[guess +i] >= x) return guess+i;
  } else for (;;guess-=roll) for(int i=0;i<roll;i++) if (array[guess-i] <= x) return guess-i;
}

template <int MIN_EQ_SZ, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
int bsLinT(const ll x, BinStruct& s) {
  const ll* array = s.a();
  auto leftIndex = 0;                                                               
  auto n = s.szA();
  while (n > MIN_EQ_SZ) {
    auto half = n / 2;
    n -= half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  auto guess = leftIndex + n/2;
  if (array[guess] < x) return baseForwardSearch(array,guess+1,x);
  else return baseBackwardSearch(array,guess,x);
}

int oracle(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  return s.a[i] == y ? i : -1;
}


int is(const ll y, IntStruct& s) {
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    m--;
    for (;; m -= 8)
      for (ll i = 0; i < 8; i++)
        if (a[m-i] <= y) return (int)(m - i);
    return (int)m;
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  for (;; m += 8)
    for (ll i = 0; i < 8; i++)
      if (a[m+i] >= y) return (int)(m+i);
  return (int)m;
}


template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
int is2(const ll y, IntStruct& s) {
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r); assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) return baseBackwardSearch(a,m-1,y);
  return baseForwardSearch(a,m,y);
}

template <int roll>
int is3(const ll y, IntStruct& s) {
#define N (r-l)*((y-a[l]) >> lgScale)
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  const unsigned lgScale = s.lgScale;
  ll m = l + dL.divFit(N,s.p, s.lg_q);
  for (;;) {
    assert(m <= r); assert(m >= l); assert(a[m] >= a[l]); // we know this because n would've been less than d
    if (a[m] > y) {
      for (ll i = 0; i < roll; i++)
        if (a[m-i] <= y) return (int)(m - i);
      r = m;
    } else {
      for (ll i = 0; i < roll; i++)
        if (a[m+i] >= y) return (int)(m+i);
      l = m;
    }
    ll d = (a[r] - a[l]) >> lgScale;
    m = l + dL.div(N,d);
  }
#undef N
}

int main(int argc, char *argv[]) {
  using std::cin; using std::istream_iterator;
  constexpr int seed = 42;
  
  int nSamples = 10;
  int nGets = -1;
  if (1 < argc) std::stringstream(argv[1]) >> nSamples;
  nSamples = nSamples < 0 ? std::numeric_limits<int>::max() : nSamples;
  if (2 < argc) std::stringstream(argv[2]) >> nGets;

  //do { IntStruct isS(input); \
  //tests.emplace_back( \
  //    benchmark<N_RUNS, IntStruct, is2>( \
  //      "bsLin", testKeys, testIndexes, isS)); } while (0);

#define RUN_IS \
  do { IntStruct isS(input); \
  tests.emplace_back( \
      benchmark<IntStruct, is2<linSIMD, linSIMD<true>>>( \
        "isSIMD", testKeys, testIndexes, isS)); } while (0);
#define RUN_BS \
  do { BinStruct bsS(input); \
  tests.emplace_back( \
      benchmark<BinStruct, bsLinT<32, linUnroll, linUnroll<true>>>( \
        "bsSIMD", testKeys, testIndexes, bsS)); } while (0);
#define RUN_OS \
  do { OracleStruct osS(input, testIndexes); \
    tests.emplace_back( \
      benchmark<OracleStruct, oracle>( \
        "os", testKeys, testIndexes, osS)); } while (0);
#define RUN_IS2 \
  do { IntStruct isS(input); \
  tests.emplace_back( \
      benchmark<IntStruct, is2<linUnroll, linUnroll<true>>>( \
        "is", testKeys, testIndexes, isS)); } while (0);

  int nNums;
  cin >> nNums;
  assert(nNums != 0);
  auto input = std::vector<ll>(istream_iterator<ll>(cin), istream_iterator<ll>());
  if (nGets < 1) nGets = input.size();

  // permute the items
  std::vector<int> testIndexes(nGets);
  std::iota(testIndexes.begin(), testIndexes.end(), 0);
  std::shuffle(testIndexes.begin(), testIndexes.end(), std::mt19937{seed});

  std::vector<ll> testKeys(testIndexes.size());
  for (int i = 0; i < testKeys.size(); i++)
    testKeys[i] = input[testIndexes[i]];

  std::vector<TestStats> tests;
  OracleStruct osS(input, testIndexes);
  T1(RUN_IS, RUN_BS, RUN_OS,RUN_IS2)
  T2(RUN_IS, RUN_BS, RUN_OS,RUN_IS2)
  T3(RUN_IS, RUN_BS, RUN_OS,RUN_IS2)

  // Set-up the headers and organize data
  bool first = true;
  for (auto& t : tests) {
    printf("%s%s", true != first ? "," : "", t.name.c_str());
#ifndef NSORT
    std::sort(t.cyclesByIx.begin(), t.cyclesByIx.end());
#endif
    first = false;
    assert(t.cyclesByIx.size() == N_RUNS);
  }
  for (int i = 0; i < N_RUNS && i < nSamples; i++) {
    first = true;
    printf("\n");
    for (auto& ts : tests) {
      printf("%s%.3f", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]) / (double)testKeys.size());
//      printf("%s%ld", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]) / (double)input.size());
      first = false;
    }
  }
  printf("\n");
}
