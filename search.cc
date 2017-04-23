#include <cmath>
#include <tuple>
#include <numeric>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <assert.h>
#include <cstdio>
#include <vector>
#include <x86intrin.h>
#include <immintrin.h>
#include <algorithm>
#include <time.h>

// could we improve binary search by dividing the tree into fixed size chunks
// and looking at the ends of the chunks? This would add additional comparisons
// but might allow us to take advantage of cache lines and pre-fetching

// TODO remember that integer type matters. Parameterize?
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
};

struct Search {
//  int ix;
  int ix, steps1, steps2;
};

template <typename StructT, typename F, F f>
void RunBenchmark( const std::vector<ll>& vals, const std::vector<int>& indexes, StructT& s, TestStats& ts) {
//void RunBenchmark( const std::vector<std::tuple<ll, int> >& order, const StructT& s, TestStats& ts) {
  unsigned A;

  const ll nNums = vals.size();

//int nEq = 0;
  int sumSteps = 0;
  ll nIx = 0;
  ll st = __rdtsc();
  for (int i = 0; i < nNums; i++) {
#ifndef NDEBUG
    printf("\n%d", indexes[i]);
#endif
    Search r = f(vals[i], s);
    nIx += r.ix;
    //bool eq = r.ix == indexes[i];
//    sumSteps += r.steps1;
    //nEq += eq;
#ifndef NDEBUG
    if (r.ix != indexes[i]) {
      printf("a[%d]=%ld <> %d\n", indexes[i], vals[i], r.ix);
      assert(r.ix == indexes[i]);
    }
#endif
  }
  ll dt = __rdtscp(&A) - st;
  ts.runStats.push_back(RunStats{sumSteps});
  ts.cyclesByIx.push_back(std::make_tuple(dt, ts.cyclesByIx.size()));
  if (nIx != nNums * (nNums - 1) / 2)
      printf("mess up\n");
  //printf("%s %llu steps, %d %llu \n", ts.name.c_str(), sumSteps, maxSteps, maxStepV);
//  if (nEq != (int)nNums)
//    printf("%d < %llu matched\n", nEq, nNums);
}


unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
}

inline unsigned lgl(uint64_t x) {
  assert(x >= 2); // subtracting and clz < 1 undefined
  return 64 - __builtin_clzll(x-1);
}

inline int lgl_flr(uint64_t x) {
  assert(x >= 1); // clz < 1 undefined
  return 64 - __builtin_clzll(x);
}

// https://locklessinc.com/articles/sat_arithmetic/
inline ll sub_sat_u64(ll x, ll y) {
  ll res = x-y;
  res &= -(res <= x);
  return res;
}

template <int LG_N>
struct DivLut {
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut() {
    memset(pT, 0, sizeof(pT[0]) * N);
    memset(lg_qT, 0, sizeof(lg_qT[0]) * N);

    // I've moved everything back by one to accomodate common case

    // multiplying by one is approximated by (2^64 - 1) / 2^64
    pT[1-1] = ~(T)0;
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

  const std::vector<ll> a;
};

struct IntStruct {
  IntStruct(const std::vector<ll>& _a) : a(_a) {
    lgScale = lg(a.size() - 1);
    ll tR = a.size() - 1;
    ll rSc = __builtin_clz(r);
    r2 = tR << rSc;
    lg_d2 = lgl(a.back() - a.front()) + rSc - 64;
    lg_d = lgl(a.back() - a.front()) - lgScale;
    yLS = a.front() >> lgScale;
    ll d = (a.back() - a.front()) >> lgScale;
    dL.one_d(d, p, lg_q);
  }

  unsigned lgScale;
  unsigned lg_d;
  unsigned lg_d2;
  ll r;
  ll r2;
  ll yLS;
  ll p; int lg_q;
  const std::vector<ll> a;
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

  const std::vector<ll> a;
  std::vector<int> i;
  int j;
};

Search bsPVKEq2(const ll x, const BinStruct& s) {
  const std::vector<ll>& array = s.a;
  const int MIN_EQ_SZ = 2;
  int leftIndex = 0;                                                               
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
  return Search{leftIndex};
}

Search is(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll d = ((yR - yL) >> lgScale);
  ll m = l + dL.div(n,d);
#ifndef NDEBUG
  printf(" %ld", m);
#endif
  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    do { m--; } while (a[m] > y);
    return Search{(int)m};
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  if (y >= yR) return Search{(int)r};
  while (a[m] < y) m++;
  return Search{(int)m};
}

Search is2(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll m = l + dL.divFit(n,s.p, s.lg_q);

#ifndef NDEBUG
  printf(" %ld", m);
#endif
  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    do { m--; } while (a[m] > y);
    return Search{(int)m};
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  if (y >= yR) return Search{(int)r};
  while (a[m] < y) m++;
  return Search{(int)m};
}

Search is3(const ll y, const IntStruct& s) {
  typedef __uint128_t u128;

  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  ll n = ((u128)s.r2 * (y-yL)) >> 64;
  ll m = l + (n >> s.lg_d2);
#ifndef NDEBUG
  printf(" %ld", m);
#endif
  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    do { m--; } while (a[m] > y);
    return Search{(int)m};
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  if (y >= yR) return Search{(int)r};
  while (a[m] < y) m++;
  return Search{(int)m};
}

Search oracle(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  return Search{i, 0}; //, 1
}

int main() {
  // lg only supports [2, 2^32-1]
  assert(lg(2) == 1);
  assert(lg(3) == 2);
  assert(lg(4) == 2);
  assert(lg(17) == 5);
  assert(lg(31) == 5);
  assert(lg(32) == 5);
#ifndef NDEBUG
  assert(0 < printf("log tests pass\n"));
#endif

  // TESTS Div
#ifndef NDEBUG
  if (dL.div(0xFFFULL, 2) != 0x7FFULL) {
    printf("0xFFF / 2 != %lx\n", dL.div(0xFFFULL, 2));
    assert(dL.div(0xFFFULL, 2) == 0x7FFULL);
  }
#endif

  int nNums;
  std::vector<ll> input, searchVal;
  std::vector<int> searchIndex;

  std::vector<ll> n, p;

  {
    std::vector<std::tuple<ll, int> > search;

    bool loaded = 1 == scanf("%d", &nNums); (void)loaded; // silence not used
    assert(loaded);
    assert(nNums > 0);
    for (int i = 0; i < nNums; i++) {
      ll x;
      loaded = 1 == scanf("%ld", &x);
      assert(loaded);
      input.push_back(x == 0 ? 1 : x); // temporary measure to keep from div by 0
      search.emplace_back(x, i);
    }
    std::srand(10);
    std::random_shuffle(search.begin(), search.end());
    for (int i = 0; i < nNums; i++) {
      ll x; int ix;
      std::tie(x, ix) = search[i];
      searchVal.push_back(x);
      searchIndex.push_back(ix);
    }
  }

  IntStruct isS(input);
  BinStruct bsS(input);
  OracleStruct oS(input, searchIndex, 0);
  // must seed first
  OracleStruct oS5(input, searchIndex, 5, true);


  typedef TestStats TS;
  //using BsFn = Search (*)(const ll, const BinStruct&);
  using IsFn = Search (*)(const ll, const IntStruct&);
  //using OsFn = Search (*)(const ll, OracleStruct&);


  //std::vector<TS> tests = { TS{"bsPVKEq2"}, TS{"is2"}
//  std::vector<TS> tests = { TS{"is"}
  std::vector<TS> tests = {
#ifndef REVERSE
     TS{"is   "}
    ,TS{"is2  "}
#else
     TS{"is2  "}
    ,TS{"is   "}
#endif
 //      TS{"bs   "}, TS{"is   "}
 //     ,TS{"bs   "}, TS{"is2  "}
//     ,TS{"bs   "}, TS{"is2  "}
//
//     ,TS{"bs   "}, TS{"is2  "}
//     ,TS{"bs   "}, TS{"is   "}
  };
  //const int N_RUNS = 1 << 5;
#ifndef NDEBUG
  const int N_RUNS = 1 << 1;
#else
  const int N_RUNS = 1000 * (1 << 16) / nNums;
#endif
  time_t lastTime = time(NULL);

  for (int i = 0; i < N_RUNS; i++) {
    time_t nowTime = time(NULL);
    if (nowTime != lastTime) {
      fprintf(stderr, "%d ", i);
      lastTime = nowTime;
    }
    int testIx = 0;
#ifndef REVERSE
    RunBenchmark<IntStruct, IsFn, is >(searchVal, searchIndex, isS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is2 >(searchVal, searchIndex, isS, tests[testIx++]);
#else
    RunBenchmark<IntStruct, IsFn, is2 >(searchVal, searchIndex, isS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is >(searchVal, searchIndex, isS, tests[testIx++]);
#endif

    assert((size_t)testIx == tests.size());
  }
  fprintf(stderr, "\n");
  bool first = true;
  for (auto& t : tests) {
    printf("%s%s", true != first ? "," : "", t.name.c_str());
    std::sort(t.cyclesByIx.begin(), t.cyclesByIx.end());
    first = false;
    assert(t.cyclesByIx.size() == N_RUNS);
    //printf("*** %s,%d ****\n", t.name.c_str(), t.runStats[0].sumSteps1);
  }
  for (int i = 0; i < N_RUNS && i < 10; i++) {
    first = true;
    printf("\n");
    for (auto& ts : tests) {
      printf("%s%ld", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]));
      first = false;
    }
  }
  printf("\n");
}
