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

struct Search {
  int ix, steps;
};

Search binSearch(const ll k, const std::vector<ll>& a) {
  unsigned l = 0;
  unsigned r = a.size();
  Search s = Search{-1, 0};

  while (r - l > 0) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    s.steps++;
    if (a[m] < k) {
      l = m + 1;
    } else if (a[m] > k) {
      r = m;
    } else {
      s.ix = m;
      break;
    }
  }

  return s;
}

struct IntStruct {
  unsigned l, r;
  double vL, vR, slope;

  IntStruct(const std::vector<ll>& a) {
    l = 0;
    r = a.size() - 1;
    vL = a[l];
    vR = a[r];
    slope = 1.0 / (vR - vL);
  }
};

Search intSearch(const ll k, const std::vector<ll>& a, IntStruct& s) {
  // can we change this to fix pt arithm?
  // first priority is reducing number of lookups
  // two ideas about systematically underestimating and overestimating
  // 1. each repeated time you under-estimate, over-estimate by more the next time.
  //    if you go on the other end, then reset and repeat.
  // 2. Take the error by which you were off, and multiply your next estimate by that.
  // 3. Look at euclid's alg and stuff like for for using error info to guide search
  int l = s.l, r = s.r;
  double vL = s.vL, vR = s.vR, vK = k;
  Search ret = Search{-1, 0};
  assert(r - l >= 0); // assume non-empty vector
  
  while (r - l > 0) {
    double p = 1.0 / (vR - vL) * (vK - vL);
    int m = l + (int)(p * (r-l));
    ret.steps++;
    if (k < a[m]) {
      // over estimate
      r = m - 1;
      vR = a[r];
    } else if (k > a[m]) {
      // under estimate
      l = m + 1;
      vL = a[l];
    } else {
      ret.ix = m;
      return ret;
    }
  }
  if (k == a[l]) {
    ret.ix = l;
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
    if (k == 6937039822968985763) {
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
          if (k == 6937039822968985763) printf("%d < %d < %d | %d li \n", l, m, newR, newR-l);
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
          if (k == 6937039822968985763) printf("%d < %d < %d | %d ri \n", newL, m, r, r-newL);
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
    if (k == 6937039822968985763) printf("%d < %d < %d | %d b \n", l, m, r, r-l);
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

//      int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)((1.0 - (1.0 / (a[r] - a[l]) * (k - a[l]))) * 2 * (r-l))));
//      m = r - leap;
//      ret.steps++;
//      if (k == 8496300069059941) printf("%d < %d < %d | %d l1 \n", l, m, r, r-l);
//      int newR = r;
//      // exponential search only works to halfway. We'd have to have been very
//      // inaccurate to get there.
//      while (leap < (r-l)/2 && k < a[m]) {
//        newR = m;
//        leap = std::min(leap * 2, (r-l)/2);
//        m = r - leap;
//        // 8496300069059941
//        if (k == 8496300069059941) printf("%d < %d < %d | %d li \n", l, m, newR, newR-l);
//        ret.steps++;
//      }
//      r = newR;
//      if (k >= a[m] && m > l) {
//        l = m;
//      }
//      assert(k <= a[r] && k >= a[l]);
//    } else if (k > a[m]) {
//      // under estimate
//      l = m + 1;
//
//      int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)(1.0 / (a[r] - a[l]) * (k - a[l]) * 2 *(r-l))));
//      m = l + leap;
//      ret.steps++;
//      if (k == 8496300069059941) printf("%d < %d < %d | %d r1 \n", l, m, r, r-l);
//
//      int newL = l;
//      // see earlier note about exp search
//      while (leap < (r-l)/2 && k > a[m]) {
//        newL = l;
//        leap = std::min(leap * 2, (r-l)/2);
//        m = l + leap;
//        ret.steps++;
//      }
//      l = newL;
//    
//      // TODO record steps?
//      if (k <= a[m] && m < r) {
//        r = m;
//      }
//      assert(k <= a[r] && k >= a[l]);
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

int main() {
  ll nNums;
  unsigned A;
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

//  {
//    std::vector<ll> array(input);
//    ll st = __rdtsc();
//    for (ll i = 0; i < nNums; i++) {
//      std::lower_bound(array.begin(), array.end(), std::get<0>(search[i]));
////      t_lb.push_back(__rdtscp(&A) - st);
//    }
//    printf("%llu\n", __rdtscp(&A) - st);
//  }
//
  {
    std::vector<ll> array(input);
    int steps = -1;
    ll st = __rdtsc();
    for (ll i = 0; i < nNums; i++) {
      steps = std::max(steps, binSearch(std::get<0>(search[i]), array).steps);
    }
    printf("%llu %d bs\n", __rdtscp(&A) - st, steps);
  }

  //{
  //  std::vector<ll> array(input);
  //  IntStruct s(array);
  //  int maxSteps = -1;
  //  ll maxStepV = -1;
  //  ll st = __rdtsc();
  //  for (ll i = 0; i < nNums; i++) {
  //    int steps = intSearch(std::get<0>(search[i]), array, s).steps;
  //    if (steps > maxSteps) {
  //      maxSteps = steps;
  //      maxStepV = std::get<0>(search[i]);
  //    }
  //  }
  //  printf("%llu %d %llu is\n", __rdtscp(&A) - st, maxSteps, maxStepV);
  //}
  {
    std::vector<ll> array(input);
    IntStruct s(array);
    int maxSteps = -1;
    ll maxStepV = -1;
    ll st = __rdtsc();
    for (ll i = 0; i < nNums; i++) {
      int steps = leapSearch(std::get<0>(search[i]), array, s).steps;
      if (steps > maxSteps) {
        maxSteps = steps;
        maxStepV = std::get<0>(search[i]);
      }
//      t_it.push_back();
    }
    printf("%llu %d %llu ls\n", __rdtscp(&A) - st, maxSteps, maxStepV);
  }

//  {
//    std::vector<ll> array(input);
//    IntStruct s(array);
//    int maxSteps = -1;
//    int maxStepIx = -1;
//    ll st = __rdtsc();
//    for (ll i = 0; i < nNums; i++) {
//      int steps = intSearch(std::get<0>(search[i]), array, s).steps;
//      if (steps > maxSteps) {
//        maxSteps = steps;
//        maxStepIx = std::get<1>(search[i]);
//      }
////      t_it.push_back();
//    }
//    printf("%llu %d %d is\n", __rdtscp(&A) - st, maxSteps, maxStepIx);
//  }

//  std::vector<PerfStats> ps;
//  ps.push_back(perfStats(t_lb, input, search, "lb"));
//  ps.push_back(perfStats(t_bs, input, search, "bs"));
//  ps.push_back(perfStats(t_it, input, search, "it"));
//
//  std::vector<std::string> fields {
//    "name", "avg", "t90", "t99", "t999", "tMax", "v90", "v99", "v999", "vMax" };
//  for (unsigned i = 0; i < fields.size(); i++) {
//    printf("%s%s", i == 0 ? "" : ",", fields[i].c_str());
//  }
//  printf("\n");
//  for (unsigned i = 0; i < ps.size(); i++) {
//    printf("%s", ps[i].name.c_str());
//    for (double t : ps[i].t) {
//      printf(",%e", t);
//    }
//    for (ll v : ps[i].v) {
//      printf(",%llu", v);
//    }
//    printf("\n");
//  }
}
