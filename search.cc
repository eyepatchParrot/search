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
    if (n < d) return 0; // can't handle shifts > 63
    // approximate divisor
    ll rN, r;
    int totalShift;
    if (d >= L) {
      const int LG_D = lgl_flr(d);
      const int k = LG_D - LG_L;
      r = 1 + ((d-1) >> k);
#ifndef NDEBUG
      if (r >= L) {
        printf("d %llu k %d r %llu \n", d, k, r);
        assert(r < L);
      }
#endif
      rN = rNT[r];
      const int lg_rD = lg_rDT[r];
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
      printf("%llu / %llu ~= %llu * %llu / 2^%d (2^%d)\n", n, d, n, rN, totalShift, lg_rDT[r]);
      printf("%llu = %llu - %llu r %llu\n", diff, n/d, q, r);
      assert(q * rN >= q);
      //assert(diff < 100);
    }
#endif
    return q;
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
//  const static ll LG_W = 16;
  const static ll W = 1ULL << LG_W;
//  const static ll LG_L = 8;
  const static ll L = 1 + (1ULL << LG_L);
  //const static ll OVERFLOW_CHK = 1ULL << (64 - LG_W);

  DivLut2() : rNT(), lg_rDT() {
    /*
     * 1/d ~= 1 / r / 2^k ~= rN / rD / 2^k
     */
    // skip zero since divide by zero makes no sense
    for (ll r = 1; r < L; r++) {
      ll rN = (1ULL << 63) / r; // 1/r ~= rN / 2^63
      assert(rN > 0); // multiply by zero doesn't make sense.

      // start counting LG_W bits from the first non-zero bit
      int lz = __builtin_clzll(rN);
      int lg_q = lz + LG_W-1;
      rN = rN >> (63 - lg_q);
      
      // we've extracted LG_W bits
      assert(rN != 0);
      int tz = __builtin_ctzll(rN);
      lg_q -= tz;
      rN = rN >> tz;

      assert(rN < W);
      assert(rN != 0);

      rNT[r] = rN;
      lg_rDT[r] = lg_q;

      if (rN == 1) {
        lg_pT[r] = 1;
      } else {
        lg_pT[r] = lgl(rN);
      }

#ifndef NDEBUG
      printf("1/%llu ~ %llu / 2^%d\n", r, rNT[r], lg_rDT[r]);
      ll rQ = (0xFFFFFFFF * rNT[r]) >> lg_rDT[r];
      ll q = 0xFFFFFFFF / r;
      assert(rQ <= q);
#endif
    }
  }
  ll rNT[L];
  int lg_rDT[L]; // 2^8 * lg(rD)
  int lg_pT[L];   // lg rNT[i]

  ll div(const ll n, const ll d) const {
    ll p;
    int lg_p, lg_q;
    one_d(d, p, lg_p, lg_q);
    return timesFrac(n, p, lg_p, lg_q);
  }

  void one_d(const ll d, ll &p, int &lg_p, int &lg_q) const {
    if (d >= L) {
      // is scaling needed
      const int LG_D = lgl_flr(d);
      const int k = LG_D - LG_L;
      const ll r = 1 + ((d-1) >> k);
      p = rNT[r];
      lg_p = lg_pT[r];
      lg_q = lg_rDT[r] + k;
    } else {
      // do lookup directly
      p = rNT[d];
      lg_p = lg_pT[d];
      lg_q = lg_rDT[d];
    }
  }

  static ll abc_d(const ll a, const int lg_a, const ll b, const int lg_b, const ll c, const int lg_c, const int lg_d) {
    int k = lg_a + lg_b + lg_c - 64;
    if (k > 0) {
      // shift right until multiplied together they sum to less than 64
      a >>= k;
      lg_d -= k;
      a = a * b * c;
      return (a * b * c) >> lg_d;
    } else {
      return (a * b * c) >> lg_d;
    }
  }


  //static ll timesShift(const ll n, const int lg_n

  static ll timesFrac(const ll n, const ll p, const int lg_p, const int lg_q) {
    // I got stuck on the overflow math before.
    // a * b overflows if lg a >= 64 - lg b
    // we might be able to approximate this by taking ~0 and shifting right
    const ll OVERFLOW_CHK = 1ULL << (64 - lg_p);
    //printf("%d \n", lg_p);
    if (n <= OVERFLOW_CHK) {
      ll q = n;
      q *= p;
      q >>= lg_q;
      return q;
    } else {
      ll q = n;
      q >>= lg_p;
      q *= p;
      q >>= (lg_q - lg_p);
      return q;
    }
  }
};

DivLut<22, 8> dL;
DivLut2<22,8> dL2;

struct Search {
  int ix, steps1, steps2;
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

//Search bs(const ll k, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  unsigned l = 0;
//  unsigned r = a.size();
//  int steps = 0;
//
//  while (r - l > 1) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    steps++;
//    if (a[m] < k) {
//      l = m + 1;
//    } else if (a[m] > k) {
//      r = m;
//    } else {
//      l = r = m;
//      //ret.ix = m;
//      //break;
//    }
//  }
//  assert(a[l] == k);
//  return Search{(int)l, steps};
//}

//Search bsNoEq(const ll k, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  unsigned l = 0;
//  unsigned r = a.size();
//  int steps = 0;
//
//  while (r - l > 1) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    steps++;
//    if (a[m] <= k) {
//      l = m;
//    } else if (a[m] > k) {
//      r = m;
//    }
//  }
//  assert(a[l] == k);
//  return Search{(int)l, steps};
//}

// PVK : https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
//Search bsPVK(const ll x, const std::vector<ll>& array, BinStruct& s) {
//  (void)s;
//  int leftIndex = 0;                                                               
//  int n = array.size();                                                            
//  int half;
//  int steps = 0;
//  if ((half = n) > 1) {
//    do {
//        half /= 2;
//        steps++;
//        n = array[leftIndex + half] == x ? 0 : n - half;
//        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
//    } while ((half = n) > 1);
//  }                                                                                
//  assert(array[leftIndex] == x);  
//  return Search{leftIndex, steps};
//}

//template<int MIN_RANK_SZ>
//Search bsRank(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  int l = 0;
//  int n = a.size();
//  int steps = 0;
//  // try also doing it with half and quarter
//  int m;
//  while ((m = n / 4) >= MIN_RANK_SZ) {
//    int rank = (a[l+m] <= y) + (a[l+2*m] <= y) + (a[l+3*m] <= y);
//    steps += 3;
//    n = a[l+3*m] <= y ? n - 3 * m : m;
//    l += m * rank;
//  }
//  while ((m = n / 2) > 0) {
//    steps++;
//    n -= a[l+m] == y ? n : m;
//    l += a[l+m] <= y ? m : 0;
//  }
//  assert(a[l] == y);
//  return Search{l, steps};
//}

//template<int MIN_RANK_SZ, int K>
//Search bsRank2(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  const int MIN_M = MIN_RANK_SZ / K;
//  int l = 0;
//  int n = a.size();
//  int steps = 0;
//  // try also doing it with half and quarter
//  int m;
//  while ((m = n / K) >= MIN_M) {
//    int rank = a[l+m] <= y;
//    for (int i = 2; i < K; i++) {
//        rank += a[l+i*m] <= y;
//    }
//    steps += (K - 1);
//    n = a[l + (K - 1) * m] <= y ? n - (K-1) * m : m;
//    l += m * rank;
//  }
//  while ((m = n / 2) > 0) {
//    steps++;
//    n -= a[l+m] == y ? n : m;
//    l += a[l+m] <= y ? m : 0;
//  }
//  assert(a[l] == y);
//  return Search{l, steps};
//}

//Search bsArrBranchless(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  // can we generate a branchless binary search?
//  // maybe use an integer to represent branching state. 
//  // could even put speculation in if there's enough depth
//  int steps = 0;
//  ll i[] = { 0, 0, a.size()};
//  while (i[2]-i[0] > 0) {
//    i[1] = (i[0] + i[2]) / 2;
//    ll x = a[i[1]] > i[1];
//    assert(x == 0 || x == 1);
//    i[0] = i[0+x];
//    i[2] = i[1+x];
//  }
//  if (a[i[0]] == y) {
//    return Search{(int)i[0], steps};
//  } else {
//    return Search{-1, steps};
//  }
//}

//Search leapSearch(const ll k, const std::vector<ll>& a, IntStruct& s) {
//  // 1. Interpolate
//  // 2. If the remaining bound is small, and not much progess is made, break;
//  // 3. Exponentially search for other bound.
//  // 4. Repeat
//  // 5. The bound is small and interpolation doesn't help much. Binary search.
//  //    No reason to go to interpolation search because bound is small.
//  // 
//  // can we change this to fix pt arithm?
//  // first priority is reducing number of lookups
//  // two ideas about systematically underestimating and overestimating
//  // 1. each repeated time you under-estimate, over-estimate by more the next time.
//  //    if you go on the other end, then reset and repeat.
//  // 2. Take the error by which you were off, and multiply your next estimate by that.
//  // 3. Look at euclid's alg and stuff like for for using error info to guide search
//  //
//  // have a notion of 'making progress'. If we interpolated and
//  // under-estimated, what is the distance from L to C? If we interpolated and
//  // over-estimated, what's the distance from R to C? If that distance is
//  // small, interpolation won't make good progress.
//  // 
//  // This can be extended to any cut. Interpolation is just one way of making
//  // a cut. We need to weigh the precision of attaining the other bound against
//  // the number of evaluations.
//  //
//  // If leap search hits half of the array, then we know that interpolation
//  // was extremely inaccurate. We should not repeat multiple bounded leap
//  // searches. 
//  //
//  // leaping expemplifies trust in the interpolation, so we can't measure
//  // progress.
//  //
//  // If we have to do log leaps, we shouldn't trust the interpolation, so do
//  // binary search.
//  int l = 0, r = a.size()-1;
//  Search ret = Search{-1, 0};
//  assert(r - l >= 0); // assume non-empty vector
//  const int MIN_SZ = 16;
//
//  while (r - l > 0) {
//    double p = 1.0 / (a[r] - a[l]) * (k - a[l]);
//    int m = l + (int)(p * (r-l));
//    ret.steps++;
//    // 4381798903807338309, 5921090195583789457, 6937039822968985763
//#ifndef NDEBUG
//    if (k == 0) {
//      printf("%d < %d < %d | %d i\n", l, m, r, r-l);
//    }
//#endif
//    if (k < a[m]) {
//      // over estimate
//      r = m - 1;
//      double progress = 1.0 - p; 
//      if (progress < 0.5) {
//        if ((r-l) <= MIN_SZ) break;
//        // potentially do x2
//        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)((1.0 - (1.0 / (a[r] - a[l]) * (k - a[l]))) * 2 * (r-l))));
//        int newR = r;
//        while (leap < (r-l)/2) {
//          m = r- leap;
//          ret.steps++;
//#ifndef NDEBUG
//          if (k == 0) printf("%d < %d < %d | %d li \n", l, m, newR, newR-l);
//#endif
//          if (k >= a[m]) {
//            break;
//          }
//          newR = m;
//          leap *=2;
//        }
//        if (leap >= (r-l)/2) {
//          r = newR;
//          break;
//        } else {
//          r = newR;
//          assert(m > l);
//          assert(m < r);
//          l = m;
//        }
//      }
//    } else if (k > a[m]) {
//      // under estimate
//      l = m + 1;
//      double progress = p;
//      if (progress < 0.5) {
//        if ((r-l) <= MIN_SZ) break;
//        //int leap = MIN_SZ;
//        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)(1.0 / (a[r] - a[l]) * (k - a[l]) * 2 * (r-l))));
//        int newL = l;
//        while (leap < (r-l)/2) {
//          m = l + leap;
//          ret.steps++;
//#ifndef NDEBUG
//          if (k == 0) printf("%d < %d < %d | %d ri \n", newL, m, r, r-newL);
//#endif
//          if (k <= a[m]) {
//            break;
//          }
//          newL = m;
//          leap *= 2;
//        }
//        if (leap >= (r-l)/2) {
//          l = newL;
//          break;
//        } else {
//          l = newL;
//          assert(m > l);
//          assert(m < r);
//          r = m;
//        }
//      }
//    } else {
//      ret.ix = m;
//      return ret;
//    }
//  }
//  while (r - l > 0) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    ret.steps++;
//#ifndef NDEBUG
//    if (k == 0) printf("%d < %d < %d | %d b \n", l, m, r, r-l);
//#endif
//    if (a[m] < k) {
//      l = m + 1;
//    } else if (a[m] > k) {
//      r = m;
//    } else {
//      ret.ix = m;
//      break;
//    }
//  }
//  if (k == a[l]) {
//    ret.ix = l;
//  }
//  return ret;
//}

// L + (R - L)(y - yL) / (yR - yL)
//Search isIntDiv(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  Search ret = Search{-1, 0};
//  ll l = 0, r = a.size() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yR = a[r], yL = a[l];
//  while (r - l > 0) {
//    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
//    ll d = ((yR - yL) >> s.lgScale);
//    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
//    ll m = l + scOff;
//    assert(m <= r);
//    assert(m >= l);
//    ret.steps++;
//    if (y < a[m]) {
//      // over estimate
//      r = m - 1;
//      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//      yL = a[l];
//    } else {
//      ret.ix = m;
//      return ret;
//    }
//  }
//  if (y == a[l]) {
//    ret.ix = l;
//  }
//
//  return ret;
//}

// L + (R - L)(y - yL) / (yR - yL)
//Search is(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  ll l = 0, r = a.size() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yR = a[r], yL = a[l];
//  const unsigned lgScale = s.lgScale;
//  int steps = 0;
//  while (r - l > 0) {
//    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
//    ll n = (r-l)*((y-yL) >> lgScale);
//    ll d = ((yR - yL) >> lgScale);
//    ll scOff = dL.div(n,d);
//    ll m = l + scOff;
//    assert(m <= r);
//    assert(m >= l);
//    steps++;
//    if (y < a[m]) {
//      // over estimate
//      r = m - 1;
//      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//      yL = a[l];
//    } else {
//      return Search{(int)m, steps};
//    }
//  }
//
//  return Search{(int)l, steps};
//}

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

Search is4(const ll y, const std::vector<ll>& a, IntStruct& s) {
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    if (n < d) break;

    //steps++;
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
//    if (m == r) {
//      while (a[r] > y && r > l) r--;
//      return Search{(int)r, steps};
//    }
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
  //ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  //  printf("y %llu\n", y);
  ll d = (a[r] - a[l]) / (r-l);
  ll m = l + (y - a[l]) / d;
  while (r - l > 0) {
#ifndef NDEBUG
    printf(" %llu", m);
#endif
    if (y < a[m]) {
      r = m-1;
      ll n = a[r] - y;
      if (n < d) {
        int oldR = r;
        while (a[r] > y) r--;
        steps += oldR - r;
        return Search{(int)r, steps};
      }
      m = r - dL.div(n, d);
      steps++;
      if (m <= l) {
        int oldL = l;
        while (a[l] < y) l++;
        steps += l - oldL;
        return Search{(int)l, steps};
      }
    } else if (y > a[m]) {
      l = m+1;
      ll n = y - a[l];
      if (n < d) {
        int oldL = l;
        while (a[l] < y) l++;
        steps += l - oldL;
        return Search{(int)l, steps};
      }
      m = l + dL.div(n,  d);
      steps++;
      if (m >= r) {
        int oldR = r;
        while (a[r] > y) r--;
        steps += oldR - r;
        return Search{(int)r, steps};
      }
    } else {
      return Search{(int)m, steps};
    }
  }

  return Search{(int)l, steps};
}


//Search is3(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  int l = 0, n = a.size()-1;
//  const unsigned lgScale = s.lgScale;
//  int steps = 0;
//  while (n > 1) {
//    const int lg_n = lg_flr(n) - 1;
//    const int R = l + (1 << lg_n);
//    ll p = y-a[l];
//    ll q = (a[R] - a[l]) >> lg_n;
//    if (p < q) break;
//
//    steps++;
//    int scOff = dL.div(p,q);
//    int m = l + scOff;
//    if (m > R) {
//      l = n - R;
//      continue;
//    }
//    assert(m <= l + n);
//    assert(m >= l);
//    if (y < a[m]) {
//      // over estimate
//      n = m - l - 1;
//      assert(l + n < m);
//    } else if (y > a[m]) {
//      // under estimate
//      const int tR = l + n;
//      n -= scOff + 1;
//      l = m+1;
//      assert(l + n == tR);
//    } else {
//      return Search{(int)m, steps};
//    }
//  }
//
//  // this too slow if not in there
//  int r = l + n;
//  while (l < r) {
//    r = y == a[l] ? l : r;
//    l++;
//  }
//
//  return Search{(int)r, steps};
//}

using BsFn = Search (*)(const ll, const std::vector<ll>&, BinStruct&);
using IsFn = Search (*)(const ll, const std::vector<ll>&, IntStruct&);

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
  if (dL2.div(0xFFFULL, 2) != 0x7FFULL) {
    printf("0xFFF / 2 != %llx\n", dL2.div(0xFFFULL, 2));
    assert(dL2.div(0xFFFULL, 2) == 0x7FFULL);
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
  std::vector<TS> tests = { TS{"bsPVKEq2"}, TS{"is2"}, TS{"bsPVKEq2"}, TS{"is3"}, TS{"bsPVKEq2"}, TS{"is4"}
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
    RunBenchmark<IntStruct, IsFn, is2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is3>(input, search, nNums, tests[testIx++]);
    RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    RunBenchmark<IntStruct, IsFn, is4>(input, search, nNums, tests[testIx++]);
    //RunBenchmark<BinStruct, BsFn, bsPVKEq2>(input, search, nNums, tests[testIx++]);
    //RunBenchmark<IntStruct, IsFn, is2>(input, search, nNums, tests[testIx++]);
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
