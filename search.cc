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
typedef int64_t ll;

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
  for (ll i = 0; i < nNums; i++) {
#ifndef NDEBUG
    printf("\n%d", indexes[i]);
#endif
    Search r = f(vals[i], s);
    nIx += r.ix;
    //bool eq = r.ix == indexes[i];
//    sumSteps += r.steps1;
    //nEq += eq;
    assert(r.ix == indexes[i]);
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
      uint64_t rN = (1ULL << 63) / r; // 1/r ~= rN / 2^63
      assert(rN > 0); // multiply by zero doesn't make sense.

      int lz = __builtin_clzll(rN);
      lg_rDT[r] = lz + LG_W-1;
      rNT[r] = rN >> (63 - lg_rDT[r]);
      assert(rNT[r] < W);
#ifndef NDEBUG
      printf("1/%ld ~ %lx / 2^%d\n", r, rNT[r], lg_rDT[r]);
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

DivLut<22, 8> dL;

struct BinStruct {
  BinStruct(const std::vector<ll>& _a) : a(_a) {  }

  const std::vector<ll> a;
};

struct IntStruct {
  IntStruct(const std::vector<ll>& _a) : a(_a) {
    lgScale = lg(a.size() - 1);
  }

  unsigned lgScale;
  const std::vector<ll> a;
};

struct OracleStruct {
  OracleStruct(const std::vector<ll>& _a, const std::vector<int>& _i) : a(_a), i(_i), j(0) { }

  const std::vector<ll> a;
  const std::vector<int> i;
  int j;
};

Search bsPVKEq2(const ll x, const BinStruct& s) {
  (void)s;
  const std::vector<ll>& array = s.a;
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

Search is(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll d = ((yR - yL) >> lgScale);
  if (n >= d) {
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %ld", m);
#endif
    assert(m <= r);
    assert(m >= l);
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      while (a[r] > y && l < r) r--;
      return Search{(int)r};
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
    } else {
      return Search{(int)m};
    }
    // TODO see if we can reduce this to an if statement and a conditional
  }

  // n < d implies that we should start from the left
  while (a[l] < y && l < r) l++;
  return Search{(int)l};
}

Search is3(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll d = ((yR - yL) >> lgScale);
  ll m = l;
  if (n >= d) {
    ll scOff = dL.div(n,d);
    m = m + scOff;
#ifndef NDEBUG
    printf(" %ld", m);
#endif
    assert(m <= r);
    assert(m >= l);

    if (y < a[m]) {
      // over estimate
      r = m - 1;
      while (a[r] > y && l < r) r--;
      return Search{(int)r};
    }
    // l = m
//    else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//    } else {
//      return Search{(int)m};
//    }
    // TODO see if we can reduce this to an if statement and a conditional
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  while (a[m] < y && m < r) m++;
  return Search{(int)m};
}

Search is2(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  //int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    if (n < d) break; // if n < d, n / d == 0, so no progress made

    //steps++;
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %ld", m);
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
      return Search{(int)m};
    }
  }

  while (a[l] < y && l < r) l++;
  return Search{(int)l};
}

Search oracle(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  return Search{i, 1};
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
      search.emplace_back(std::make_tuple(x, i));
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
  OracleStruct oS(input, searchIndex);


  typedef TestStats TS;
  using BsFn = Search (*)(const ll, const BinStruct&);
  using IsFn = Search (*)(const ll, const IntStruct&);
  using OsFn = Search (*)(const ll, OracleStruct&);


  //std::vector<TS> tests = { TS{"bsPVKEq2"}, TS{"is2"}
  std::vector<TS> tests = {
    TS{"bsPVKEq2"}, TS{"is3"}
    ,TS{"bsPVKEq2"}, TS{"is"}
    ,TS{"bsPVKEq2"}, TS{"oracle"}

    ,TS{"bsPVKEq2"}, TS{"oracle"}
    ,TS{"bsPVKEq2"}, TS{"is"}
    ,TS{"bsPVKEq2"}, TS{"is3"}
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
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is3 >(searchVal, searchIndex, isS, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is >(searchVal, searchIndex, isS, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<OracleStruct, OsFn, oracle>(searchVal, searchIndex, oS, tests[testIx++]);

    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<OracleStruct, OsFn, oracle>(searchVal, searchIndex, oS, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is >(searchVal, searchIndex, isS, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(searchVal, searchIndex, bsS, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is3 >(searchVal, searchIndex, isS, tests[testIx++]);
    /*
     *
bsPVKEq2,is2
64472,48707
64486,48730
64489,48745
64489,48754
64492,48771
64497,48776
64500,48779
64501,48783
64504,48783
64506,48785


bsPVKEq2,is2,bsT
63764,61076,50571
63805,61155,50586
63811,61164,50589
63820,61190,50601
63823,61193,50601
63825,61198,50607
63826,61201,50607
63828,61234,50609
63828,61236,50610
63828,61245,50610

4426 
bsPVKEq2,is2,bsPVKEq2,bsT,bsPVKEq2,bsT,bsPVKEq2,is2
63732,59718,63724,64291,63765,64215,63706,60212
63744,59837,63729,64294,63770,64227,63724,60224
63750,59945,63756,64317,63773,64230,63729,60238
63773,59962,63767,64326,63776,64239,63733,60305
63788,59988,63776,64331,63784,64241,63741,60320
63793,59997,63782,64332,63784,64242,63744,60369
63794,60009,63794,64332,63785,64242,63764,60378
63799,60012,63794,64335,63790,64244,63765,60384
63802,60032,63802,64343,63793,64244,63767,60387
63805,60035,63805,64346,63793,64245,63770,60390

not using tuple
5211 
bsPVKEq2,is2,bsPVKEq2,bsT,bsPVKEq2,bsT,bsPVKEq2,is2
64483,56512,64520,64989,64523,64931,64486,56768
64527,56515,64541,65038,64547,64975,64547,56835
64530,56561,64544,65044,64559,64977,64550,56846
64532,56576,64544,65062,64573,64980,64558,56847
64538,56643,64556,65065,64576,64980,64562,56896
64544,56687,64556,65065,64579,64992,64567,56905
64544,56701,64559,65076,64582,64992,64576,56934
64544,56733,64564,65082,64582,64992,64576,56934
64547,56739,64564,65088,64582,64992,64576,56943
64550,56759,64567,65088,64584,64992,64582,56969

292 5142 
bsPVKEq2,is2,bsPVKEq2,bsT,bsPVKEq2,oracle,bsPVKEq2,oracle,bsPVKEq2,bsT,bsPVKEq2,is2
58563,61347,58569,58752,58540,30630,58630,30595,58601,58760,58531,61905
58581,61347,58572,58752,58542,30645,58633,30604,58607,58772,58537,61973
58584,61359,58575,58766,58546,30653,58633,30606,58607,58778,58537,61993
58586,61367,58578,58767,58549,30653,58635,30607,58610,58778,58537,62033
58592,61370,58581,58772,58557,30653,58636,30607,58610,58778,58539,62056
58592,61388,58583,58773,58557,30656,58636,30607,58613,58781,58539,62069
58592,61391,58584,58781,58557,30656,58636,30610,58613,58781,58540,62080
58592,61408,58586,58781,58560,30656,58639,30612,58615,58781,58540,62080
58592,61422,58589,58781,58560,30656,58639,30613,58615,58781,58540,62086
58595,61425,58589,58784,58560,30659,58642,30613,58618,58784,58540,62092

1219 6072 
bsPVKEq2,is2,bsPVKEq2,bsT,bsPVKEq2,oracle,bsPVKEq2,oracle,bsPVKEq2,bsT,bsPVKEq2,is2
58575,61996,58589,58731,58531,30342,58534,30336,58525,58645,58534,62388
58581,61998,58592,58732,58540,30345,58546,30342,58534,58665,58534,62397
58583,62004,58592,58735,58540,30351,58549,30345,58534,58668,58537,62414
58586,62005,58592,58740,58540,30354,58551,30345,58537,58670,58539,62417
58592,62039,58592,58740,58540,30354,58551,30347,58540,58679,58540,62420
58592,62045,58595,58741,58540,30354,58554,30348,58540,58685,58540,62421
58595,62051,58595,58749,58540,30356,58557,30348,58546,58685,58542,62426
58595,62051,58595,58749,58543,30357,58560,30348,58548,58685,58543,62426
58595,62065,58595,58749,58543,30357,58560,30350,58552,58688,58543,62444
58598,62071,58595,58749,58543,30357,58560,30351,58554,58688,58545,62461


89 5075 
bsPVKEq2,is2,bsPVKEq2,bsT,bsPVKEq2,oracle,bsPVKEq2,oracle,bsPVKEq2,bsT,bsPVKEq2,is2
58574,52873,58574,58723,58539,30357,58531,30350,58531,58659,58519,53309
58575,52884,58580,58732,58545,30360,58531,30351,58531,58662,58531,53315
58577,52916,58583,58737,58545,30362,58534,30351,58534,58668,58531,53318
58577,52922,58586,58741,58546,30362,58546,30353,58534,58676,58543,53324
58578,52922,58589,58743,58549,30362,58548,30353,58537,58679,58543,53327
58583,52925,58589,58743,58551,30362,58549,30359,58540,58679,58546,53332
58586,52948,58592,58746,58552,30365,58551,30359,58545,58679,58546,53344
58586,52957,58592,58746,58554,30365,58551,30360,58546,58682,58546,53344
58586,53000,58592,58749,58554,30365,58551,30360,58546,58685,58546,53350
58586,53006,58592,58749,58554,30365,58551,30362,58546,58685,58546,53353
*/
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
      printf("%s%ld", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]));
      first = false;
    }
  }
  printf("\n");
}
