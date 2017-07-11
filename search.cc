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

typedef uint64_t ll;

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
      cyclesByIx.push_back(std::make_tuple(dt, cyclesByIx.size()));
  }
  void datum(ll st, ll et) { datum(et-st); }
  void data(TestStats& ts) {
      for (auto t_i : ts.cyclesByIx)
          datum(std::get<0>(t_i));
  }
};

template <int nRuns, typename S, int (f) (const ll, S& s)>
TestStats benchmark(
    std::string&& name,
    const std::vector<ll>& vals,
    const std::vector<int>& indexes,
    S& s) {
  const int N_NUMS = vals.size();
  auto ts = TestStats{name};
  ll sumIx = 0;
  ll sum0_n = 0;
  for (int runIx = 0; runIx < nRuns; runIx++) {
    ll st = __rdtsc();

    for (int i = 0; i < N_NUMS; i++) {
#ifndef NDEBUG
      fprintf(stderr,"%d,%d,", i,indexes[i]);
#endif
      int ix = f(vals[i], s);
      sumIx += ix;
      assert(ix == indexes[i]);
    }

    ll et = __rdtsc();
    ts.datum(st, et);
    sum0_n += N_NUMS * (N_NUMS - 1) / 2;
  }
  if (sumIx != sum0_n) // check even when profiling
    fprintf(stderr, "mess up\n");
  return ts;
}

template <int LG_N>
struct DivLut {
  // Note: All indices are off by one to save a subtraction during division.
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut() {
    memset(pT, 0, sizeof(pT[0]) * N);
    memset(lg_qT, 0, sizeof(lg_qT[0]) * N);

    pT[1-1] = ~(T)0; // x * 1 ~= x * 2^64 - 1 / 2^64
    lg_qT[1-1] = 0;

    // do powers of 2
    // start at zero since we multiply by 2^64 / 2
    for (int i = 2, j = 0; i <= N; i*=2, j++) {
      pT[-1+i] = (T)1 << 63;
      lg_qT[-1+i] = j;
    }

    // generate fast multiply results
    for (int i = 1; i <= N; i++) {
      int d = i;
      int lg_q = 0;
      // reduce to an odd number
      while (pT[-1+d] == 0 && d % 2 == 0) {
        d /= 2;
        lg_q++;
      }
      if (pT[-1+d] != 0) {
        pT[-1+i] = pT[-1+d];
        lg_qT[-1+i] = lg_qT[-1+d] + lg_q;
        continue;
      }

      // find a 2^n + 1 multiple
      int j = 64;
      for (; j < 64+lgl(d); j++) {
        DT x = ((DT)1 << j) + 1;
        int r = x % d;
        if (r == 0) {
          pT[-1+i] = (((DT)1 << j) + d) / d;
          lg_qT[-1+i] = lg_q + (j - 64);
          break;
        }
      }
      if (pT[-1+i] == 0) {
        pT[-1+i] = (((T)1 << 63) / d) << 1;
        // guaranteed i > 1, so we will always have room to shift left
        // we always divide by 2^64, so we have to make sure that's happening
        lg_qT[-1+i] = lg_q;
      }
    }
  }

  T div(T n, T d) {
    d--; // start the d-1 early
    if (d > N-1) {
      // note that this becomes ceiling log
      const int k = lgl_flr(d) - LG_N;
      const T pIx = d >> k; // 1 + (d-1) >> k, but we started d, and -1 for index
      return divFit(n, pT[pIx], k + lg_qT[pIx]);
    } else {
      return divFit(n, pT[d], lg_qT[d]);
    }
  }

  void one_d(const T d, T& p, int& lg_q) const {
    if (d > N) {
      const int k = lgl_flr(d) - LG_N;
      ll pIx = (d-1) >> k;
      p = pT[pIx];
      lg_q = k + lg_qT[pIx];
    } else {
      p = pT[d];
      lg_q = lg_qT[d];
    }
  }

  T divFit(T n, T p, int lg_q) {
    // assuming that we're in the range
    T hi = ((DT)n * p) >> 64;
    return hi >> lg_q;
  }

  T pT[N];
  int lg_qT[N];
};

template <int LG_N>
struct DivLut3 {
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut3() {
    memset(pT, 0, sizeof(pT[0]) * N);

    // counting on lg_q being small enough that we can shift right

    // I've moved everything back by one to accomodate common case

    // multiplying by one is approximated by (2^64 - 1) / 2^64
    pT[1-1] = ~(T)0;

    // do powers of 2
    // start at zero since we multiply by 2^64 / 2
    for (int i = 2, lg_q = 0; i <= N; i*=2, lg_q++) {
      T p = (T)1 << 63;
      p = p >> lg_q;
      pT[-1+i] = p;
    }

    // generate fast multiply results
    for (int i = 1; i <= N; i++) {
      int d = i;
      int lg_q = 0;
      // reduce to an odd number
      while (pT[-1+d] == 0 && d % 2 == 0) {
        d /= 2;
        lg_q++;
      }
      if (pT[-1+d] != 0) {
        if (i == d) continue;
        pT[-1+i] = pT[-1+d] >> lg_q;
        continue;
      }

      // find a 2^n + 1 multiple
      int j = 64;
      for (; j < 64+lgl(d); j++) {
        DT x = ((DT)1 << j) + 1;
        int r = x % d;
        if (r == 0) {
          lg_q = lg_q + (j-64);
          pT[-1+i] = ((((DT)1 << j) + d) / d) >> lg_q;
          break;
        }
      }
      if (pT[-1+i] == 0) {
        pT[-1+i] = ((((T)1 << 63) / d) << 1) >> lg_q;
      }
    }
  }

  T div(T n, T d) {
    d--; // start the d-1 early
    if (d > N-1) {
      // note that this becomes ceiling log
      const int k = lgl_flr(d) - LG_N;
      const T p = d >> k; // +1 for math, -1 for offset
      return divFit(n, p) >> k;
    } else {
      // -1 for index adjustment + 1 for early adjustment
      return divFit(n, d); 
    }
  }

  T divFit(T n, T p) {
    assert(p <= N);
    // assuming that we're in the range
    return ((DT)n * pT[p]) >> 64;
  }

  T pT[N];
};

DivLut<8> dL;
DivLut3<8> dL2;

struct BinStruct {
  BinStruct(const std::vector<ll>& _a) : a(_a) {  }

  const std::vector<ll>& a;
};

struct IntStruct {
  IntStruct(const std::vector<ll>& _a) : v(_a.size() + 16,0) {
    std::copy(_a.begin(), _a.end(), v.begin() + 8);
    for (int i = 0; i < 8; i++) v[i] = std::numeric_limits<ll>::min();
    for (int i = _a.size() + 8; i < v.size(); i++) v[i] = std::numeric_limits<ll>::max();
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
  size_t szA() { return v.size() - 16; };
  const ll* a() { return &v[8]; };
};

struct IntStructT {
  IntStructT(const std::vector<ll>& _a) : array(_a) {
    lgScale = lg(a.size() - 1);
    a.push_back(_a[0] >> lgScale);
    for (int i = 1; i < _a.size(); i++) {
      a.push_back(_a[i] >> lgScale);
      assert(_a[i] == _a[i-1] || (_a[i] >> lgScale) > (_a[i-1] >> lgScale));
    }
  }

  unsigned lgScale;
  std::vector<ll> a;
  const std::vector<ll> array;
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
  const std::vector<ll>& array = s.a;
  const int MIN_EQ_SZ = 2;
  long leftIndex = 0;                                                               
  int n = array.size();                                                            
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

int bsPVKEq3(const ll x, BinStruct& s) {
  const std::vector<ll>& array = s.a;
  const int MIN_EQ_SZ = 2;
  long leftIndex = 0;                                                               
  int n = array.size();                                                            
  while (n > MIN_EQ_SZ) {
    int half = n / 2;
    n -= half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  while (n > 1) {
    int half = n / 2;
    n = array[leftIndex + half] == x ? 0 : n - half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }                                                                                
  assert(array[leftIndex] == x);  
  return leftIndex;
}

int isSIMD(const ll y, IntStruct& s) {
  typedef uint64_t v4u __attribute__ ((vector_size (32)));
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll m = l + dL.divFit(n,s.p, s.lg_q);
#ifndef NDEBUG
  fprintf(stderr,"%ld\n", m);
#endif
  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    while (m >= 4) {
      m -= 4;
      const ll* v = a + m;
      int off = 0;
      const v4u ymm1 = _mm256_cmpgt_epi64((v4u){v[off++], v[off++], v[off++], v[off++]}, (v4u){y,y,y,y});
      const int t1 = _mm256_movemask_epi8(ymm1);
      const int t2 = ~t1;
      if (t2) {
        const int t3 = __builtin_clz(t2);
        const int t4 = t3 / 8;
        const int t5 = 3 + m - t4;
        return t5;
      }
    }
    do { m--; } while (a[m] > y);
    return (int)m;
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  //while (m + 3 <= r) {
  //  const ll* v = a.data() + m;
  //  // Compare a vector of the data with a vector of the searched for element
  //  // Gather together the results of each comparison into one byte per element
  //  // is mask ever not true?
  //  int off = 0;
  //  const v4u ymm1 =
  //    _mm256_cmpgt_epi64((v4u){y,y,y,y}, (v4u){v[off++], v[off++], v[off++], v[off++]});
  //  const int t1 = _mm256_movemask_epi8(ymm1);
  //  const int t2 = ~t1;
  //  //const v4u ymm2 =
  //  //  _mm256_cmpgt_epi64((v4u){y,y,y,y}, (v4u){v[off++], v[off++], v[off++], v[off++]});
  //  //const int ea1 =
  //  //  _mm256_testc_si256(
  //  //      ymm1,
  //  //      (v4u){(ll)-1, (ll)-1, (ll)-1, (ll)-1});
  //  //const int ea2 =
  //  //  _mm256_testc_si256(
  //  //      ymm2,
  //  //      (v4u){(ll)-1, (ll)-1, (ll)-1, (ll)-1});
  //  //const int ea3 = ea1 & ea2; // look for any zero-bits

  //  if (t2) {
  //    //ll t1 = (uint32_t)_mm256_movemask_epi8(ymm1);
  //    //ll t2 = (uint32_t)_mm256_movemask_epi8(ymm2);
  //    //ll t3 = t2 << 32;
  //    //ll t4 = t3 | t1;
  //    //int t5 = __builtin_clzl(t4);
  //    //int t6 = t5 / 8;
  //    //int t7 = 8 + (int)m - t6;
  //    int t3 = __builtin_ctz(t2);
  //    int t4 = t3 / 8;
  //    int t5 = m + t4;
  //    return t5;
  //  }
  //  m += 4;
  //}

  //if (y >= yR) return (int)r;
  while (a[m] < y) m++;
  return (int)m;
}

int is2(const ll y, IntStruct& s) {
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

int oracle(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  return s.a[i] == y ? i : -1;
}

int main(int argc, char *argv[]) {
  using std::cin; using std::istream_iterator;
  constexpr int seed = 42;
  
  int nSamples = 10;
  if (argc >= 2)
    std::stringstream(argv[1]) >> nSamples;
  nSamples = nSamples < 0 ? std::numeric_limits<int>::max() : nSamples;

#define RUN_IS \
  do { IntStruct isS(input); \
  tests.emplace_back( \
      benchmark<N_RUNS, IntStruct, is2>( \
        "is", testKeys, testIndexes, isS)); } while (0);
#define RUN_BS \
  do { BinStruct bsS(input); \
  tests.emplace_back( \
      benchmark<N_RUNS, BinStruct, bsPVKEq3>( \
        "bs", testKeys, testIndexes, bsS)); } while (0);
#define RUN_OS \
  do { OracleStruct osS(input, testIndexes); \
    tests.emplace_back( \
      benchmark<N_RUNS, OracleStruct, oracle>( \
        "os", testKeys, testIndexes, osS)); } while (0);
#define RUN_IS2 \
  do { IntStruct isS(input); \
  tests.emplace_back( \
      benchmark<N_RUNS, IntStruct, isSIMD>( \
        "isSIMD", testKeys, testIndexes, isS)); } while (0);

  int nNums;
  cin >> nNums;
  assert(nNums != 0);
  auto input = std::vector<ll>(istream_iterator<ll>(cin), istream_iterator<ll>());

  // permute the items
  std::vector<ll> testKeys;
  std::vector<int> testIndexes(input.size());
  testKeys.reserve(input.size());
  testIndexes.reserve(input.size());
  std::iota(testIndexes.begin(), testIndexes.end(), 0);
  std::shuffle(testIndexes.begin(), testIndexes.end(), std::mt19937{seed});
  for (int ix : testIndexes)
    testKeys.push_back(input[ix]);

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
      printf("%s%ld", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]));
      first = false;
    }
  }
  printf("\n");
}
