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

// could we improve binary search by dividing the tree into fixed size chunks
// and looking at the ends of the chunks? This would add additional comparisons
// but might allow us to take advantage of cache lines and pre-fetching

// TODO remember that integer type matters. Parameterize?
typedef unsigned long long ll;

unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

inline unsigned lgl(ll x) {
  assert(x >= 2); // subtracting and clz < 1 undefined
  return 64 - __builtin_clzll(x-1);
}

inline unsigned lgl_flr(ll x) {
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
//  const static ll LG_W = 16;
  const static ll W = 1ULL << LG_W;
//  const static ll LG_L = 8;
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
      printf("1/%llu ~ %llu / 2^%llu\n", r, rNT[r], lg_rDT[r]);
      ll rQ = (0xFFFFFFFF * rNT[r]) >> lg_rDT[r];
      ll q = 0xFFFFFFFF / r;
      assert(rQ <= q);
#endif
    }
  }
  ll rNT[L];
  ll lg_rDT[L]; // 2^8 * lg(rD)

  ll div(const ll n, const ll d) const {
    if (n < d) return 0; // can't handle shifts > 63
    // approximate divisor
    ll rN, totalShift, r;
    if (d >= L) {
      const ll LG_D = lgl_flr(d);
      const ll k = LG_D - LG_L;
      r = 1 + ((d-1) >> k);
#ifndef NDEBUG
      if (r >= L) {
        printf("d %llu k %llu r %llu \n", d, k, r);
        assert(r < L);
      }
#endif
      rN = rNT[r];
      const ll lg_rD = lg_rDT[r];
      totalShift = k + lg_rD;
    } else {
      r = d;
      rN = rNT[r];
      totalShift = lg_rDT[r];
    }

    // calculate quotient
    ll q = quot(n, rN, totalShift);
#ifndef NDEBUG
    ll diff = (n/d) - q;
    if (diff > 10) {
      printf("%llu / %llu ~= %llu * %llu / 2^%llu (2^%llu)\n", n, d, n, rN, totalShift, lg_rDT[r]);
      printf("%llu = %llu - %llu r %llu\n", diff, n/d, q, r);
      assert(q * rN >= q);
      //assert(diff < 100);
    }
#endif
    return q;
  }

  static inline ll quot(const ll n, const ll rN, const ll lg_d) {
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

struct Search {
  int ix, steps;
};

struct BinStruct {
  BinStruct(const std::vector<ll>& a) { (void)a; }
};

struct IntStruct {
  unsigned lgScale;

  IntStruct(const std::vector<ll>& a) {
    lgScale = lg(a.size() - 1);
  }
};

Search bs(const ll k, const std::vector<ll>& a, BinStruct& s) {
  (void)s;
  unsigned l = 0;
  unsigned r = a.size();
  int steps = 0;

  while (r - l > 1) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    steps++;
    if (a[m] < k) {
      l = m + 1;
    } else if (a[m] > k) {
      r = m;
    } else {
      l = r = m;
      //ret.ix = m;
      //break;
    }
  }
  assert(a[l] == k);
  return Search{(int)l, steps};
}

Search bsNoEq(const ll k, const std::vector<ll>& a, BinStruct& s) {
  (void)s;
  unsigned l = 0;
  unsigned r = a.size();
  int steps = 0;

  while (r - l > 1) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    steps++;
    if (a[m] <= k) {
      l = m;
    } else if (a[m] > k) {
      r = m;
    }
  }
  assert(a[l] == k);
  return Search{(int)l, steps};
}

// PVK : https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
Search bsPVK(const ll x, const std::vector<ll>& array, BinStruct& s) {
  (void)s;
  int leftIndex = 0;                                                               
  int n = array.size();                                                            
  int half;
  int steps = 0;
  if ((half = n) > 1) {
    do {
        half /= 2;
        steps++;
        n = array[leftIndex + half] == x ? 0 : n - half;
        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    } while ((half = n) > 1);
  }                                                                                
  assert(array[leftIndex] == x);  
  return Search{leftIndex, steps};
}

Search bsPVKEq2(const ll x, const std::vector<ll>& array, BinStruct& s) {
  (void)s;
  const int MIN_EQ_SZ = 2;
  int leftIndex = 0;                                                               
  int n = array.size();                                                            
  int half;
  int steps = 0;
  if ((half = n) > MIN_EQ_SZ) {
    do {
        half /= 2;
        steps++;
        n -= half;
        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    } while ((half = n) > MIN_EQ_SZ);
  }                                                                                
  if ((half = n) > 1) {
    do {
        half /= 2;
        steps++;
        n = array[leftIndex + half] == x ? 0 : n - half;
        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    } while ((half = n) > 1);
  }                                                                                
  assert(array[leftIndex] == x);  
  return Search{leftIndex, steps};
}

template<int MIN_RANK_SZ>
Search bsRank(const ll y, const std::vector<ll>& a, BinStruct& s) {
  (void)s;
  int l = 0;
  int n = a.size();
  int steps = 0;
  // try also doing it with half and quarter
  int m;
  while ((m = n / 4) >= MIN_RANK_SZ) {
    int rank = (a[l+m] <= y) + (a[l+2*m] <= y) + (a[l+3*m] <= y);
    steps += 3;
    n = a[l+3*m] <= y ? n - 3 * m : m;
    l += m * rank;
  }
  while ((m = n / 2) > 0) {
    steps++;
    n -= a[l+m] == y ? n : m;
    l += a[l+m] <= y ? m : 0;
  }
  assert(a[l] == y);
  return Search{l, steps};
}

template<int MIN_RANK_SZ, int K>
Search bsRank2(const ll y, const std::vector<ll>& a, BinStruct& s) {
  (void)s;
  const int MIN_M = MIN_RANK_SZ / K;
  int l = 0;
  int n = a.size();
  int steps = 0;
  // try also doing it with half and quarter
  int m;
  while ((m = n / K) >= MIN_M) {
    int rank = a[l+m] <= y;
    for (int i = 2; i < K; i++) {
        rank += a[l+i*m] <= y;
    }
    steps += (K - 1);
    n = a[l + (K - 1) * m] <= y ? n - (K-1) * m : m;
    l += m * rank;
  }
  while ((m = n / 2) > 0) {
    steps++;
    n -= a[l+m] == y ? n : m;
    l += a[l+m] <= y ? m : 0;
  }
  assert(a[l] == y);
  return Search{l, steps};
}

Search bsArrBranchless(const ll y, const std::vector<ll>& a, BinStruct& s) {
  // can we generate a branchless binary search?
  // maybe use an integer to represent branching state. 
  // could even put speculation in if there's enough depth
  int steps = 0;
  ll i[] = { 0, 0, a.size()};
  while (i[2]-i[0] > 0) {
    i[1] = (i[0] + i[2]) / 2;
    ll x = a[i[1]] > i[1];
    assert(x == 0 || x == 1);
    i[0] = i[0+x];
    i[2] = i[1+x];
  }
  if (a[i[0]] == y) {
    return Search{(int)i[0], steps};
  } else {
    return Search{-1, steps};
  }
}

Search leapSearch(const ll k, const std::vector<ll>& a, IntStruct& s) {
  // 1. Interpolate
  // 2. If the remaining bound is small, and not much progess is made, break;
  // 3. Exponentially search for other bound.
  // 4. Repeat
  // 5. The bound is small and interpolation doesn't help much. Binary search.
  //    No reason to go to interpolation search because bound is small.
  // 
  // can we change this to fix pt arithm?
  // first priority is reducing number of lookups
  // two ideas about systematically underestimating and overestimating
  // 1. each repeated time you under-estimate, over-estimate by more the next time.
  //    if you go on the other end, then reset and repeat.
  // 2. Take the error by which you were off, and multiply your next estimate by that.
  // 3. Look at euclid's alg and stuff like for for using error info to guide search
  //
  // have a notion of 'making progress'. If we interpolated and
  // under-estimated, what is the distance from L to C? If we interpolated and
  // over-estimated, what's the distance from R to C? If that distance is
  // small, interpolation won't make good progress.
  // 
  // This can be extended to any cut. Interpolation is just one way of making
  // a cut. We need to weigh the precision of attaining the other bound against
  // the number of evaluations.
  //
  // If leap search hits half of the array, then we know that interpolation
  // was extremely inaccurate. We should not repeat multiple bounded leap
  // searches. 
  //
  // leaping expemplifies trust in the interpolation, so we can't measure
  // progress.
  //
  // If we have to do log leaps, we shouldn't trust the interpolation, so do
  // binary search.
  int l = 0, r = a.size()-1;
  Search ret = Search{-1, 0};
  assert(r - l >= 0); // assume non-empty vector
  const int MIN_SZ = 16;

  while (r - l > 0) {
    double p = 1.0 / (a[r] - a[l]) * (k - a[l]);
    int m = l + (int)(p * (r-l));
    ret.steps++;
    // 4381798903807338309, 5921090195583789457, 6937039822968985763
#ifndef NDEBUG
    if (k == 0) {
      printf("%d < %d < %d | %d i\n", l, m, r, r-l);
    }
#endif
    if (k < a[m]) {
      // over estimate
      r = m - 1;
      double progress = 1.0 - p; 
      if (progress < 0.5) {
        if ((r-l) <= MIN_SZ) break;
        // potentially do x2
        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)((1.0 - (1.0 / (a[r] - a[l]) * (k - a[l]))) * 2 * (r-l))));
        int newR = r;
        while (leap < (r-l)/2) {
          m = r- leap;
          ret.steps++;
#ifndef NDEBUG
          if (k == 0) printf("%d < %d < %d | %d li \n", l, m, newR, newR-l);
#endif
          if (k >= a[m]) {
            break;
          }
          newR = m;
          leap *=2;
        }
        if (leap >= (r-l)/2) {
          r = newR;
          break;
        } else {
          r = newR;
          assert(m > l);
          assert(m < r);
          l = m;
        }
      }
    } else if (k > a[m]) {
      // under estimate
      l = m + 1;
      double progress = p;
      if (progress < 0.5) {
        if ((r-l) <= MIN_SZ) break;
        //int leap = MIN_SZ;
        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)(1.0 / (a[r] - a[l]) * (k - a[l]) * 2 * (r-l))));
        int newL = l;
        while (leap < (r-l)/2) {
          m = l + leap;
          ret.steps++;
#ifndef NDEBUG
          if (k == 0) printf("%d < %d < %d | %d ri \n", newL, m, r, r-newL);
#endif
          if (k <= a[m]) {
            break;
          }
          newL = m;
          leap *= 2;
        }
        if (leap >= (r-l)/2) {
          l = newL;
          break;
        } else {
          l = newL;
          assert(m > l);
          assert(m < r);
          r = m;
        }
      }
    } else {
      ret.ix = m;
      return ret;
    }
  }
  while (r - l > 0) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    ret.steps++;
#ifndef NDEBUG
    if (k == 0) printf("%d < %d < %d | %d b \n", l, m, r, r-l);
#endif
    if (a[m] < k) {
      l = m + 1;
    } else if (a[m] > k) {
      r = m;
    } else {
      ret.ix = m;
      break;
    }
  }
  if (k == a[l]) {
    ret.ix = l;
  }
  return ret;
}

// L + (R - L)(y - yL) / (yR - yL)
Search isIntDiv(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  while (r - l > 0) {
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    ll d = ((yR - yL) >> s.lgScale);
    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
    ll m = l + scOff;
    assert(m <= r);
    assert(m >= l);
    ret.steps++;
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      ret.ix = m;
      return ret;
    }
  }
  if (y == a[l]) {
    ret.ix = l;
  }

  return ret;
}

// L + (R - L)(y - yL) / (yR - yL)
Search is(const ll y, const std::vector<ll>& a, IntStruct& s) {
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
    assert(m <= r);
    assert(m >= l);
    steps++;
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      // TODO try removing eq check
      return Search{(int)m, steps};
    }
  }

  return Search{(int)l, steps};
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

    ll scOff = dL.div(n,d);
    ll m = l + scOff;
    assert(m <= r);
    assert(m >= l);
    steps++;
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

  // this too slow if not in there
  while (l < r) {
    steps++;
    r = y == a[l] ? l : r;
    l++;
  }

  return Search{(int)r, steps};
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

    ll scOff = dL.div(n,d);
    if (scOff == 0) break; // if n < d, n / d == 0, so no progress made
    ll m = l + scOff;
    assert(m <= r);
    assert(m >= l);
    steps++;
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m;
//      yL = a[l];
    } else {
      l = m;
      yL = a[l];
    }
  }

  // this too slow if not in there
  while (l < r) {
    //steps++;
    r = y == a[l] ? l : r;
    l++;
  }

  return Search{(int)r, steps};
}

Search intSearch(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0;
  int r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
#ifndef NDEBUG
  int bad = 0;
#endif
  
  while (r - l > 0) {
    ll yR = a[r], yL = a[l];
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    unsigned k = lgl(yR - yL) - s.lgScale;
    ll scOff = ((r-l) * ((y - yL) >> s.lgScale)) >> k;
    assert(scOff < 1ULL << 32);
    int m = l + scOff;
#ifndef NDEBUG
    ll off = (y - yL) / ((yR - yL) / (r-l));
    int m2 = l + off;
    if ((m - m2) * (m-m2) > 1) {
      printf("[k l m r] --> [%llu %d %d %d]\n", y, l, m-m2, r);
      assert(++bad < 100);
    }
#endif
    assert(m <= r);
    assert(m >= l);
    ret.steps++;
    if (y < a[m]) {
      // over estimate
      r = m - 1;
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
    } else {
      ret.ix = m;
      return ret;
    }
  }
  if (y == a[l]) {
    ret.ix = l;
  }

  return ret;
}

struct PerfStats {
  std::string name;
  std::vector<double> t;
  std::vector<ll> v;
};

PerfStats perfStats(std::vector<ll>& time, const std::vector<ll>& array, const std::vector<std::tuple<ll, ll> >& order, const std::string& name) {
  ll nNums = array.size();
  std::vector<std::tuple<ll, int> > time_ix;
  ll sum = 0;
  for (ll i = 0; i < nNums; i++) {
    time_ix.push_back(std::tuple<ll, ll>(time.back(), std::get<1>(order[i])));
    //time_ix.push_back(std::make_tuple(time.back(), std::get<1>(order[i])));
    sum += time.back();
    time.pop_back();
  }

  std::sort(time_ix.begin(), time_ix.end());
  int ix90 = 0.9 * nNums, ix99 = 0.99 * nNums, ix999 = 0.999 * nNums;

  PerfStats ret;
  ret.name = name;
  ret.t.push_back((double)sum / nNums);
  ret.t.push_back(std::get<0>(time_ix[ix90]));
  ret.t.push_back(std::get<0>(time_ix[ix99]));
  ret.t.push_back(std::get<0>(time_ix[ix999]));
  ret.t.push_back(std::get<0>(time_ix.back()));
  ret.v.push_back(std::get<1>(time_ix[ix90]));
  ret.v.push_back(std::get<1>(time_ix[ix99]));
  ret.v.push_back(std::get<1>(time_ix[ix999]));
  ret.v.push_back(std::get<1>(time_ix.back()));
  return ret;
}

using BenchFn = Search (*)(const ll, const std::vector<ll>&, IntStruct&);

struct RunStats {
};

struct TestStats {
  std::string name;
  std::vector<RunStats> runStats;
  std::vector<std::tuple<ll, int> > cyclesByIx;
};

template <typename StructT, typename F>
void RunBenchmark(const std::vector<ll>& input, const std::vector<std::tuple<ll, int> >& order, ll nNums, F f, TestStats& ts) {
  unsigned A;

  std::vector<ll> array(input);
  StructT s(array);
//  int maxSteps = -1;
//  ll maxStepV = -1;
//  ll sumSteps = 0;
  int nEq = 0;
  //ll nIx = 0;
  ll st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    Search r = f(std::get<0>(order[i]), array, s);
    //nIx += r.ix;
    bool eq = r.ix == std::get<1>(order[i]);
    nEq += eq;
    assert(eq);
//    sumSteps += r.steps;
//    if (r.steps > maxSteps) {
//      maxSteps = r.steps;
//      maxStepV = std::get<0>(order[i]);
//    }
  }
  ll dt = __rdtscp(&A) - st;
  ts.runStats.push_back(RunStats{});
  ts.cyclesByIx.push_back(std::make_tuple(dt, ts.cyclesByIx.size()));
//  if (nIx != nNums * (nNums - 1) / 2)
//      printf("mess up\n");
  //printf("%s %llu steps, %d %llu \n", ts.name.c_str(), sumSteps, maxSteps, maxStepV);
  if (nEq != (int)nNums)
    printf("%d < %llu matched\n", nEq, nNums);
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

//  ll const * const inPtr = input.data();
//  unsigned A;
//  ll sum = 0;
//  printf("name,cycles\n");
//  ll st;
//  for (int i = 0; i < 10; i++) {
//    st = __rdtsc();
//    for (int i = 0; i < nNums; i++) {
//      for (int j = 0; j < nNums; j++) {
//        sum += input[i] / input[j];
//      }
//    }
//    ll it = __rdtscp(&A) - st;
//    printf("intDiv,%llu\n", it);
//    st = __rdtsc();
//    for (int i = 0; i < nNums; i++) {
//      for (int j = 0; j < nNums; j++) {
//        sum += dL.div(inPtr[i], inPtr[j]);
//      }
//    }
//    ll lt = __rdtscp(&A);
//    printf("LUT2048,%llu\n", lt-st);
//    st = __rdtsc();
//    for (int i = 0; i < nNums; i++) {
//      for (int j = 0; j < nNums; j++) {
//        sum += dL.div(inPtr[i], inPtr[j]);
//      }
//    }
//    ll lt2 = __rdtscp(&A);
//    printf("LUT2048,%llu\n", lt2-st);
//    st = __rdtsc();
//    for (int i = 0; i < nNums; i++) {
//      for (int j = 0; j < nNums; j++) {
//        sum += dL2.div(inPtr[i], inPtr[j]);
//      }
//    }
//    ll lt3 = __rdtscp(&A);
//    printf("LUT256,%llu\n", lt3-st);
//  }
//  printf("sum,%llu\n", sum);


  std::vector<TestStats> tests = { TestStats{"bs"}, TestStats{"bsNoEq"}, TestStats{"bsPVK"}, TestStats{"bsPVKEq2"}, TestStats{"is"}, TestStats{"is2"}, TestStats{"is3"}
    //TestStats{"isLUTDiv", (void*)isLUTDiv}
  };
  const int N_RUNS = 1 << 15;
  for (int i = 0; i < N_RUNS; i++) {
    RunBenchmark<BinStruct>(input, search, nNums, bs, tests[0]);
    RunBenchmark<BinStruct>(input, search, nNums, bsNoEq, tests[1]);
    RunBenchmark<BinStruct>(input, search, nNums, bsPVK, tests[2]);
    RunBenchmark<BinStruct>(input, search, nNums, bsPVKEq2, tests[3]);
    RunBenchmark<IntStruct>(input, search, nNums, is, tests[4]);
    RunBenchmark<IntStruct>(input, search, nNums, is2, tests[5]);
    RunBenchmark<IntStruct>(input, search, nNums, is3, tests[6]);
//    RunBenchmark<BinStruct>(input, search, nNums, "bsRank16_4", bsRank2<16, 4>);
//    RunBenchmark<BinStruct>(input, search, nNums, "bsRank16_6", bsRank2<16, 4>);
    //RunBenchmark<IntStruct>(input, search, nNums, "isIntDiv", isIntDiv);
    //RunBenchmark<IntStruct>(input, search, nNums, "isLUTDiv", isLUTDiv);
    //RunBenchmark<IntStruct>(input, search, nNums, "isLUTDiv2", isLUTDiv2);
  }
  bool first = true;
  for (auto& t : tests) {
    printf("%s%s", true != first ? "," : "", t.name.c_str());
    std::sort(t.cyclesByIx.begin(), t.cyclesByIx.end());
    first = false;
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
    //RunBenchmark<BinStruct>(input, search, nNums, "bs", binSearch);
    //RunBenchmark<BinStruct>(input, search, nNums, "bs4", bs4);
    //RunBenchmark<BinStruct>(input, search, nNums, "bs3", bs3);
    //RunBenchmark<IntStruct>(input, search, nNums, "isFp", isFp);
    //RunBenchmark<IntStruct>(input, search, nNums, "isDiv", isDiv);
    //RunBenchmark<IntStruct>(input, search, nNums, "isDiv2", isDiv2);
    //RunBenchmark<IntStruct>(input, search, nNums, "isScLn", isScLn);
    //RunBenchmark<IntStruct>(input, search, nNums, "isExc", isExc);
    //RunBenchmark<IntStruct>(input, search, nNums, "is", intSearch);
    //RunBenchmark<IntStruct>(input, search, nNums, "ls", leapSearch);
}
