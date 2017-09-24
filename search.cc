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
DivLut4 dL4;

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

template <typename S, ll (f) (const ll, S& s)>
TestStats benchmark(
    std::string&& name,
    const std::vector<ll>& vals,
    const std::vector<int>& indexes,
    S& s) {
  constexpr int nRuns = N_RUNS;
  const auto N_NUMS = vals.size();
  auto ts = TestStats{name};
  auto expSum = 0UL, valSum = 0UL;
  for (auto j=0;j<nRuns;j++) for (auto i : indexes) expSum += vals[i];
  for (int runIx = 0; runIx < nRuns; runIx++) {
    ll st = __rdtsc();

    for (int i = 0; i < N_NUMS; i++) {
      ERR("%d,%d,", i,indexes[i]);
      auto val = f(vals[i], s);
      valSum += val;
      assert(val == vals[i]);
    }

    ll et = __rdtsc();
    ts.datum(st, et);
  }
  if (valSum != expSum) // check even when profiling
    fprintf(stderr, "mess up %lu\n", vals.back());
  return ts;
}


struct BinStruct {
  BinStruct(const std::vector<ll>& _a) : v(_a.size() + 32,0) {
    std::copy(_a.begin(), _a.end(), v.begin() + 16);
    for (int i = 0; i < 16; i++) v[i] = std::numeric_limits<ll>::min();
    for (int i = _a.size() + 16; i < v.size(); i++) v[i] = std::numeric_limits<int64_t>::max();
    n_bs = lgl(_a.size()) - 5;
    assert(n_bs > 0);
  }
  long n_bs; 
  auto szA() { return v.size() - 32; };
  const ll* a() { return &v[16]; };
  private:
  std::vector<ll> v;
};

struct IntStruct {
  auto szA() { return v.size() - 64; };
  const ll* a() { return &v[32]; };

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
    ll d = (_a.back() - _a.front()) >> lgScale; // consider getting rid of this
    dL.one_d(d, p, lg_q);
    p2 = (p >> lgScale) * (szA() - 1);
    p3 = p2 >> lg_q;
    d4 = DivLut4::Divisor(d);
    d4.p = (d4.p >> lgScale) * (szA() - 1);
  }

  unsigned lgScale;
  unsigned lg_d;
  unsigned lg_d2;
  ll r;
  ll r2;
  ll yLS;
  ll p; int lg_q;
  ll p2; int lg_q2;
  ll p3;
  DivLut4::Divisor d4;
  std::vector<ll> v;
};

// TODO programatically find the size of the pruned array
struct PruneStruct {
  auto a() { return &va[32]; };
  auto szA() { return va.size() - 64; };
  auto aBack() { return a()[szA()-1]; };
  auto b() { return &vb[32]; };
  auto szB() { return vb.size() - 64; };

  PruneStruct(const std::vector<ll>& _a) {
    // a is guaranteed to be sorted, so we can generate each element myopically
    auto vaSz = (int)(_a.size() * 0.801);
    auto slope = (double)(vaSz-1) / (_a.back() - _a.front());
    va.reserve(vaSz + 64);
    vb.reserve(_a.size());
    for (auto i=0;i<32;i++) va.push_back(std::numeric_limits<ll>::min());
    for (auto i=0;i<32;i++) vb.push_back(std::numeric_limits<ll>::min());
    va.push_back(_a[0]);
    for (auto i=1,j=1;i<vaSz;i++) {
      while (j<_a.size() && (ll)(slope * (_a[j] - _a.front())) <= i) vb.push_back(_a[j++]);
      //while (j<_a.size() && (vb.size() <= 32 || (ll)(slope * (_a[j] - _a.front())) <= i)) vb.push_back(_a[j++]);
      assert(vb.size() > 32);
      va.push_back(vb.back());
      vb.pop_back();
    }
    for (auto i=0;i<32;i++) va.push_back(std::numeric_limits<int64_t>::max());
    for (auto i=0;i<32;i++) vb.push_back(std::numeric_limits<int64_t>::max());
    std::sort(va.begin(), va.end());
    assert(std::is_sorted(vb.begin(), vb.end()));
    lgScale = lg(vaSz - 1);
    ll d = (aBack() - a()[0]) >> lgScale;
    dL.one_d(d, p, lg_q);
    
    // pre-compute max dist
    maxDistF = 0;
    maxDistB = 0;
    for (auto i=0L;i<szA();i++) {
      auto x = a()[i];
      auto dist = i-(long)(slope * (x - _a.front()));
      if (dist < maxDistB) maxDistB = dist;
      if (dist > maxDistF) maxDistF = dist;
    }
    maxDistB = -maxDistB;
    assert(maxDistB == 0);
  }

  std::vector<ll> va;
  std::vector<ll> vb;
  unsigned lgScale;
  ll p;
  int lg_q;
  int maxDistF;
  int maxDistB;
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

using SearchFn = int64_t(const ll*, int64_t, ll);
template < bool reverse=false, int roll=2>
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

template <bool reverse=false,int n=8>
int64_t linUnroll(const ll* a, int64_t m, ll k) {
  for (;;m = (reverse?m-n:m+n)) {
    for (int i = 0; i < n; i++) {
      assert(m+i < 1032); assert((m-i) > -32);
      if (reverse?(a[m-i]<=k):(a[m+i]>=k)) return reverse?(m-i):(m+i);
    }
  }
}

template <int MIN_EQ_SZ, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto binLin(const ll* a,  int leftIndex, int n, const ll x) {
  while (n > MIN_EQ_SZ) {
    auto half = n / 2;
    n -= half;
    leftIndex = a[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  auto guess = leftIndex + n/2;
  if (a[guess] < x) return baseForwardSearch(a,guess+1,x);
  else return baseBackwardSearch(a,guess,x);
}

auto bsPVKEq2(const ll x, BinStruct& s) {
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
  return array[leftIndex];
}

template <int MIN_EQ_SZ, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto bsLinT2(const ll x, BinStruct& s) {
  auto n = s.szA();
  auto a = s.a();
  auto l = 0L;
  for (auto i = 0L; i < s.n_bs; i++) {
    auto half = n/2;
    n -= half;
    l = a[l+half] <= x ? l+half : l;
  };
  auto guess = l + n/2;
  if (a[guess] < x) return a[baseForwardSearch(a,guess+1,x)];
  return a[baseBackwardSearch(a,guess,x)];
}

ll bs(const ll k, BinStruct& s) {
  auto a = s.a();
  unsigned l = 0;
  unsigned r = s.szA();

  while (r - l > 1) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = l + (r-l) / 2;
    if (a[m] < k) {
      l = m + 1;
    } else if (a[m] > k) {
      r = m;
    } else {
      l = r = m;
    }
  }
  assert(a[l] == k);
  return a[l];
}

auto oracle(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  return s.a[i] == y ? s.a[i] : -1;
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto is2(const ll y, IntStruct& s) {
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r); assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) return a[baseBackwardSearch(a,m-1,y)];
  return a[baseForwardSearch(a,m,y)];
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto is3(const ll y, IntStruct& s) {
  const ll* a = s.a();
  ll l = 0;
  assert(s.szA() - l >= 0); // assume non-empty vector
  ll yL = a[l];
  //ll n = (r-l)*((y-yL) >> lgScale);
  //ll m = l + dL.divFit(n,s.p, s.lg_q);
  ll m = l + dL.divFit(y-yL, s.p2, s.lg_q);

  assert(m <= s.szA()); assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) return a[baseBackwardSearch(a,m-1,y)];
  return a[baseForwardSearch(a,m,y)];
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto is4(const ll y, IntStruct& s) {
  constexpr ll l = 0;

  const ll* a = s.a();
  assert(s.szA() - l >= 0); // assume non-empty vector
  ll m = l + (y-a[l]) / s.d4;

  assert(m <= s.szA()); assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) return a[baseBackwardSearch(a,m-1,y)];
  return a[baseForwardSearch(a,m,y)];
}

template <SearchFn* sF, SearchFn* sB, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps5(const ll y, PruneStruct& s) {
  auto a = s.a();
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  auto m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  assert(a[m] <= y);
  m = sF(a,m,y);
  if (a[m] == y) return a[m];

  auto b = s.b();
  assert(s.szA() / s.szB() == 4);
  auto m2 = m >> 2;
  if (b[m2] > y) return b[baseBackwardSearch(b,m2-1,y)]; 
  return b[baseForwardSearch(b, m2, y)];
}

// ./x [nRuns] [subsetSize]
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
      benchmark<IntStruct, is3<linSIMD<false>, linSIMD<true>>>( \
        "is3", testKeys, testIndexes, isS)); } while (0);
#define RUN_IS2 \
  do { IntStruct isS(input); \
  tests.emplace_back( \
      benchmark<IntStruct, is4<linSIMD, linSIMD<true>>>( \
        "is4", testKeys, testIndexes, isS)); } while (0);
#if BS_T == 1
#define RUN_BS \
  do { BinStruct bsS(input); \
  tests.emplace_back( \
      benchmark<BinStruct, bsPVKEq2>( \
        "bs", testKeys, testIndexes, bsS)); } while (0);
#elif BS_T == 2
#define RUN_BS \
  do { BinStruct bsS(input); \
  tests.emplace_back( \
      benchmark<BinStruct, bsLinT2<32, linUnroll, linUnroll<true>>>( \
        "bsLin", testKeys, testIndexes, bsS)); } while (0);
#elif BS_T == 3
#define RUN_BS \
  do { BinStruct bsS(input); \
  tests.emplace_back( \
      benchmark<BinStruct, bs>( \
        "bs", testKeys, testIndexes, bsS)); } while (0);
#endif
#define RUN_OS \
  do { OracleStruct osS(input, testIndexes); \
    tests.emplace_back( \
      benchmark<OracleStruct, oracle>( \
        "os", testKeys, testIndexes, osS)); } while (0);

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
