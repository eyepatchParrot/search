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

unsigned lgl(ll x) {
  assert(x >= 2); // subtracting and clz < 1 undefined
  return 64 - __builtin_clzll(x-1);
}

struct Search {
  int ix, steps;
};

struct BinStruct {
  BinStruct(const std::vector<ll>& a) { }
};

struct IntStruct {
  unsigned lgScale;

  IntStruct(const std::vector<ll>& a) {
    lgScale = lg(a.size() - 1);
  }
};

Search binSearch(const ll k, const std::vector<ll>& a, BinStruct& s) {
  unsigned l = 0;
  unsigned r = a.size();
  Search ret = Search{-1, 0};

  while (r - l > 0) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    ret.steps++;
    if (a[m] < k) {
      l = m + 1;
    } else if (a[m] > k) {
      r = m;
    } else {
      ret.ix = m;
      break;
    }
  }

  return ret;
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


Search isFp(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0;
  int r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  double yR = a[r], yL = a[l], yY = y;
  
  while (r - l > 0) {
    double p = 1.0 / (yR - yL) * (yY - yL);
    int m = l + (int)(p * (r-l));
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


Search isDiv(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0;
  int r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  while (r - l > 0) {
    ll yR = a[r], yL = a[l];
    ll off = (y - yL) / ((yR - yL) / (r-l));
    int m = l + off;
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

// L + (R - L)(y - yL) / (yR - yL)
Search isSc1(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  while (r - l > 0) {
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    ll d = ((yR - yL) >> s.lgScale);
    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
    int m = l + scOff;
    assert(m <= r);
    assert(m >= l);
//#ifndef NDEBUG
//    ll yS = y >> s.lgScale, yLS = yL >> s.lgScale, yRS = yR >> s.lgScale;
//    ll scOff2 = r * ((y - yL) >> s.lgScale) + l * ((yR - y) >> s.lgScale);
//    int m2 = scOff2 / d;
//    ll scOff3 = r * yS - r * yLS + l * yRS - l * yS;
//    int m3 = scOff3 / (yRS - yLS);
//    int m4 = r * ((y - yL) >> s.lgScale) / d + l * ((yR - y) >> s.lgScale) / d;
//    int m5 = (r * yS / d) - (r * yLS / d) + (l * yRS / d) - (l * yS / d);
//    if ((m - m2) * (m-m2) > 1) printf("%d m2\n", m - m2);
//    if ((m - m3) * (m-m3) > 1) printf("%d m3\n", m - m3);
//    if ((m - m4) * (m-m4) > 1) printf("%d m4\n", m - m4);
//    if ((m - m5) * (m-m5) > 1) printf("%d m5\n", m - m5);
//#endif
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

Search isSc2(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  while (r - l > 0) {
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    ll d = ((yR - yL) >> s.lgScale);
    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
    int m = l + scOff;
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
Search isScLn(const ll y, const std::vector<ll>& a, IntStruct& s) {
  Search ret = Search{-1, 0};
  int l = 0;
  int r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  while (r - l > 2) {
    ll yR = a[r], yL = a[l];
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    ll d = ((yR - yL) >> s.lgScale);
    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
    int m = l + scOff;
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
  while (r-l >= 0) {
    ret.steps++;
    if (y > a[l]) {
      l++;
    } else {
      ret.ix = l;
      return ret;
    }
  }
  return ret;
}

// C = (R - L - 2)(y - yL)/(yR - yL)
Search isExc(const ll y, const std::vector<ll>& a, IntStruct& s) {
  // TODO try using r, l, m as long longs
  int l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yL = a[l], yR = a[r];
  if (yL == y) {
    return Search{l, 2};
  } else if (yR == y) {
    return Search{r, 2};
  }
  // TODO try removing ret
  Search ret = Search{-1, 2};
  while (r - l > 2) {
    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
    ll n = ((y - yL) >> s.lgScale), d = ((yR - yL) >> s.lgScale), w = (r-l-2);
    int m = l + 1 + w * n / d;
    assert(m < r);
    assert(m > l);
    ret.steps++;
    if (y < a[m]) {
      // over estimate
      r = m;
      yR = a[m];
    } else if (y > a[m]) {
      // under estimate
      l = m;
      yL = a[m];
    } else {
      ret.ix = m;
      return ret;
    }
  }
  return ret;
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

template <typename StructT, typename F>
void RunBenchmark(std::vector<ll>& input, std::vector<std::tuple<ll, ll> >& order, ll nNums, const char* name, F f) {
  unsigned A;

  std::vector<ll> array(input);
  StructT s(array);
  int maxSteps = -1;
  ll maxStepV = -1;
  ll sumSteps = 0;
  ll st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    int steps = f(std::get<0>(order[i]), array, s).steps;
    sumSteps += steps;
    if (steps > maxSteps) {
      maxSteps = steps;
      maxStepV = std::get<0>(order[i]);
    }
  }
  printf("%llu %llu %d %llu %s\n", __rdtscp(&A) - st, sumSteps, maxSteps, maxStepV, name);
}

int main() {
  // lg only supports [2, 2^32-1]
  assert(lg(2) == 1);
  assert(lg(3) == 2);
  assert(lg(4) == 2);
  assert(lg(17) == 5);
  assert(lg(31) == 5);
  assert(lg(32) == 5);
  assert(0 < printf("log tests pass\n"));
  ll nNums;
  std::vector<ll> input;
  std::vector<std::tuple<ll, ll> > search;

  bool loaded = 1 == scanf("%llu", &nNums); (void)loaded; // silence not used
  assert(loaded);
  assert(nNums > 0);
  for (ll i = 0; i < nNums; i++) {
    ll x;
    loaded = 1 == scanf("%llu", &x);
    assert(loaded);
    input.push_back(x);
    search.push_back(std::tuple<ll, ll>(x, i));
  }
  std::vector<ll> t_lb(nNums), t_bs(nNums);

  std::srand(10);
  std::random_shuffle(search.begin(), search.end());

  /*
//  {
//    std::vector<ll> array(input);
//    ll st = __rdtsc();
//    for (ll i = 0; i < nNums; i++) {
//      std::lower_bound(array.begin(), array.end(), std::get<0>(search[i]));
//    }
//    printf("%llu\n", __rdtscp(&A) - st);
//  }
  */
  RunBenchmark<BinStruct>(input, search, nNums, "bs", binSearch);
  RunBenchmark<IntStruct>(input, search, nNums, "isFp", isFp);
  RunBenchmark<IntStruct>(input, search, nNums, "isDiv", isDiv);
  //RunBenchmark<IntStruct>(input, search, nNums, "isDiv2", isDiv2);
  RunBenchmark<IntStruct>(input, search, nNums, "isSc1", isSc1);
  RunBenchmark<IntStruct>(input, search, nNums, "isSc2", isSc2);
  RunBenchmark<IntStruct>(input, search, nNums, "isScLn", isScLn);
  RunBenchmark<IntStruct>(input, search, nNums, "isExc", isExc);
  //RunBenchmark<IntStruct>(input, search, nNums, "is", intSearch);
  //RunBenchmark<IntStruct>(input, search, nNums, "ls", leapSearch);
}
