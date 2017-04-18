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
#include <algorithm>
#include <time.h>

// could we improve binary search by dividing the tree into fixed size chunks
// and looking at the ends of the chunks? This would add additional comparisons
// but might allow us to take advantage of cache lines and pre-fetching

// TODO remember that integer type matters. Parameterize?
typedef unsigned long long ll;

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
  int ix, steps1, steps2;
};

template <typename StructT, typename F, F f>
void RunBenchmark(const std::vector<ll>& input, const std::vector<std::tuple<ll, int> >& order, ll nNums, TestStats& ts) {
  unsigned A;

  std::vector<ll> array(input);
  StructT s(array);
  int nEq = 0;
  int sumSteps = 0;
  //ll nIx = 0;
  ll st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
#ifndef NDEBUG
    printf("\n%d", std::get<1>(order[i]));
#endif
    Search r = f(std::get<0>(order[i]), array, s);
    //nIx += r.ix;
    bool eq = r.ix == std::get<1>(order[i]);
//    sumSteps += r.steps1;
    nEq += eq;
    assert(eq);
  }
  ll dt = __rdtscp(&A) - st;
  ts.runStats.push_back(RunStats{sumSteps});
  ts.cyclesByIx.push_back(std::make_tuple(dt, ts.cyclesByIx.size()));
//  if (nIx != nNums * (nNums - 1) / 2)
//      printf("mess up\n");
  //printf("%s %llu steps, %d %llu \n", ts.name.c_str(), sumSteps, maxSteps, maxStepV);
  if (nEq != (int)nNums)
    printf("%d < %llu matched\n", nEq, nNums);
}


unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
}

inline unsigned lgl(ll x) {
  assert(x >= 2); // subtracting and clz < 1 undefined
  return 64 - __builtin_clzll(x-1);
}

inline int lgl_flr(ll x) {
  assert(x >= 1); // clz < 1 undefined
  return 64 - __builtin_clzll(x);
}

// https://locklessinc.com/articles/sat_arithmetic/
inline ll sub_sat_u64(ll x, ll y) {
  ll res = x-y;
  res &= -(res <= x);
  return res;
}

// http://stackoverflow.com/questions/19016099/lookup-table-with-constexpr
template <int LG_W, int LG_L>
struct DivLut {
  const static ll W = 1ULL << LG_W;
  const static ll L = 1 + (1ULL << LG_L);
  const static ll OVERFLOW_CHK = 1ULL << (64 - LG_W);

  DivLut() : rNT(), lg_rDT() {
    /*
     * 1/d ~= 1 / r / 2^k ~= rN / rD / 2^k
     */
    // skip zero since divide by zero makes no sense
    for (ll r = 1; r < L; r++) {
      ll rN = (1ULL << 63) / r; // 1/r ~= rN / 2^63
      assert(rN > 0); // multiply by zero doesn't make sense.

      int lz = __builtin_clzll(rN);
      lg_rDT[r] = lz + LG_W-1;
      rNT[r] = rN >> (63 - lg_rDT[r]);
      assert(rNT[r] < W);
#ifndef NDEBUG
      printf("1/%llu ~ %llx / 2^%d\n", r, rNT[r], lg_rDT[r]);
      ll rQ = (0xFFFFFFFF * rNT[r]) >> lg_rDT[r];
      ll q = 0xFFFFFFFF / r;
      assert(rQ <= q);
#endif
    }
  }
  ll rNT[L];
  int lg_rDT[L]; // 2^8 * lg(rD)
  
  ll div(const ll n, const ll d) const {
    assert(n >= d);
//    if (n < d) return 0; // can't handle shifts > 63
    // approximate divisor
    ll rN;
    int totalShift;
    one_d(d, rN, totalShift);

    // calculate quotient
    ll q = quot(n, rN, totalShift);
    return q;
  }

  void one_d(const ll d, ll &p, int &lg_q) const {
    if (d >= L) {
      // is scaling needed
      const int LG_D = lgl_flr(d);
      const int k = LG_D - LG_L;
      const ll r = 1 + ((d-1) >> k);
      p = rNT[r];
      lg_q = lg_rDT[r] + k;
    } else {
      // do lookup directly
      p = rNT[d];
      lg_q = lg_rDT[d];
    }
  }

  static inline ll quot(const ll n, const ll rN, const int lg_d) {
    if (n <= OVERFLOW_CHK) {
      ll q = n;
      q *= rN;
      q >>= lg_d;
      return q;
    } else {
      ll q = n;
      q >>= LG_W;
      q *= rN;
      q >>= (lg_d - LG_W);
      return q;
    }
  }
};

template <int LG_W, int LG_L>
struct DivLut2 {
  const static ll W = 1ULL << LG_W;
  const static ll L = 1 + (1ULL << LG_L);
  const static ll OVERFLOW_CHK = 1ULL << (64 - LG_W);

  DivLut2() {
    /*
     * 1/d ~= 1 / r / 2^k ~= rN / rD / 2^k
     */
    // skip zero since divide by zero makes no sense
    for (ll r = 1; r < L; r++) {
      ll rN = (1ULL << 63) / r; // 1/r ~= rN / 2^63
      assert(rN > 0); // multiply by zero doesn't make sense.

      int lz = __builtin_clzll(rN);
      int lg_q = lz + LG_W-1;
      ll p = rN >> (63 - lg_q);
      assert(p < W);

      t[r] = { p, lg_q };
    }
  }
  
  ll div(const ll n, const ll d) const {
    assert(n >= d);
//    if (n < d) return 0; // can't handle shifts > 63
    // approximate divisor
    ll p;
    int lg_q;
    one_d(d, p, lg_q);

    // calculate quotient
    ll r = quot(n, p, lg_q);
    return r;
  }

  void one_d(const ll d, ll &p, int &lg_q) const {
    // TODO simplify into two separate functions
    if (d >= L) {
      // is scaling needed
      const int lg_d = lgl_flr(d);
      const int k = lg_d - LG_L;
      const ll r = 1 + ((d-1) >> k);
      p = t[r].p;
      lg_q = t[r].lg_q + k;
    } else {
      // do lookup directly
      p = t[d].p;
      lg_q = t[d].lg_q;
    }
  }

  static inline ll quot(const ll n, const ll rN, const int lg_d) {
    // TODO simplify into two functions
    if (n <= OVERFLOW_CHK) {
      ll q = n;
      q *= rN;
      q >>= lg_d;
      return q;
    } else {
      ll q = n;
      q >>= LG_W;
      q *= rN;
      q >>= (lg_d - LG_W);
      return q;
    }
  }

  struct {
      ll p;
      int lg_q;
  } t[L];
};


DivLut<22, 8> dL;
DivLut2<22, 8> dL2;

struct BinStruct {
  BinStruct(const std::vector<ll>& a) { (void)a; }
};

struct IntStruct {
  unsigned lgScale;

  IntStruct(const std::vector<ll>& a) {
    lgScale = lg(a.size() - 1);
  }
};

Search bsPVKEq2(const ll x, const std::vector<ll>& array, BinStruct& s) {
  (void)s;
  const int MIN_EQ_SZ = 2;
  int leftIndex = 0;                                                               
  int n = array.size();                                                            
  int half;
  if ((half = n) > MIN_EQ_SZ) {
    do {
        half /= 2;
        n -= half;
        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    } while ((half = n) > MIN_EQ_SZ);
  }                                                                                
  if ((half = n) > 1) {
    do {
        half /= 2;
        n = array[leftIndex + half] == x ? 0 : n - half;
        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    } while ((half = n) > 1);
  }                                                                                
  assert(array[leftIndex] == x);  
  return Search{leftIndex};
}

Search is2(const ll y, const std::vector<ll>& a, IntStruct& s) {
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    if (n < d) break; // if n < d, n / d == 0, so no progress made

    //steps++;
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %llu", m);
#endif
    assert(m <= r);
    assert(m >= l);
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      return Search{(int)m, steps};
    }
  }

  while (a[l] < y && l < r) l++;
  return Search{(int)l, steps};
}

Search is3(const ll y, const std::vector<ll>& a, IntStruct& s) {
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    if (n < d) break; // if n < d, n / d == 0, so no progress made

    //steps++;
    ll scOff = dL2.div(n,d);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %llu", m);
#endif
    assert(m <= r);
    assert(m >= l);
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      return Search{(int)m, steps};
    }
  }

  while (a[l] < y && l < r) l++;
  return Search{(int)l, steps};
}

int main() {
  using BsFn = Search (*)(const ll, const std::vector<ll>&, BinStruct&);
  using IsFn = Search (*)(const ll, const std::vector<ll>&, IntStruct&);

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
    printf("0xFFF / 2 != %llx\n", dL.div(0xFFFULL, 2));
    assert(dL.div(0xFFFULL, 2) == 0x7FFULL);
  }
#endif

  int nNums;
  std::vector<ll> input;
  std::vector<std::tuple<ll, int> > search;

  bool loaded = 1 == scanf("%d", &nNums); (void)loaded; // silence not used
  assert(loaded);
  assert(nNums > 0);
  for (int i = 0; i < nNums; i++) {
    ll x;
    loaded = 1 == scanf("%llu", &x);
    assert(loaded);
    input.push_back(x == 0 ? 1 : x); // temporary measure to keep from div by 0
    search.push_back(std::tuple<ll, int>(x, i));
  }
  std::vector<ll> t_lb(nNums), t_bs(nNums);

  std::srand(10);
  std::random_shuffle(search.begin(), search.end());

  typedef TestStats TS;

  //std::vector<TS> tests = { TS{"bsPVKEq2"}, TS{"is2"}
  std::vector<TS> tests = {
    TS{"bsPVKEq2"}, TS{"is2"}
    ,TS{"bsPVKEq2"}, TS{"isDL"}
    ,TS{"bsPVKEq2"}, TS{"is2"}
    ,TS{"bsPVKEq2"}, TS{"isDL"}
  };
  //const int N_RUNS = 1 << 5;
#ifndef NDEBUG
  const int N_RUNS = 1 << 1;
#else
  const int N_RUNS = 1 << 13;
#endif
  time_t lastTime = time(NULL);

  for (int i = 0; i < N_RUNS; i++) {
    time_t nowTime = time(NULL);
    if (nowTime != lastTime) {
      fprintf(stderr, "%d ", i);
      lastTime = nowTime;
    }
    int testIx = 0;
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is3>(input, search, nNums, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is2 >(input, search, nNums, tests[testIx++]);

    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is2 >(input, search, nNums, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is3>(input, search, nNums, tests[testIx++]);
    assert((size_t)testIx == tests.size());
  }
  fprintf(stderr, "\n");
  bool first = true;
  for (auto& t : tests) {
    printf("%s%s", true != first ? "," : "", t.name.c_str());
    std::sort(t.cyclesByIx.begin(), t.cyclesByIx.end());
    first = false;
    assert(t.cyclesByIx.size() == N_RUNS);
//    printf("%s,%d\n", t.name.c_str(), t.runStats[0].sumSteps1);
  }
  for (int i = 0; i < N_RUNS && i < 10; i++) {
    first = true;
    printf("\n");
    for (auto& ts : tests) {
      printf("%s%llu", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]));
      first = false;
    }
  }
  printf("\n");
}
