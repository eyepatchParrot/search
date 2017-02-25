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
  ll l = s.l, r = s.r;
  double vL = s.vL, vR = s.vR, vK = k;
  Search ret = Search{-1, 0};
  assert(r - l >= 0); // assume non-empty vector
  while (r - l > 0) {
    double p =  1.0 / (vR - vL) * (vK - vL);
    unsigned m = std::min(r, l + (unsigned)(p * (r-l)));
//    if (k == a[12185329]) {
//      printf("%u\n", m);
//    }
    ret.steps++;
    if (k < a[m]) {
      r = m - 1;
      vR = a[r];
    } else if (k > a[m]) {
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
//  {
//    std::vector<ll> array(input);
//    int steps = -1;
//    ll st = __rdtsc();
//    for (ll i = 0; i < nNums; i++) {
//      steps = std::max(steps, binSearch(std::get<0>(search[i]), array).steps);
////      t_bs.push_back(__rdtscp(&A) - st);
//    }
//    printf("%llu %d bs\n", __rdtscp(&A) - st, steps);
//  }

  {
    std::vector<ll> array(input);
    IntStruct s(array);
    int maxSteps = -1;
    int maxStepIx = -1;
    ll st = __rdtsc();
    for (ll i = 0; i < nNums; i++) {
      int steps = intSearch(std::get<0>(search[i]), array, s).steps;
      if (steps > maxSteps) {
        maxSteps = steps;
        maxStepIx = std::get<1>(search[i]);
      }
//      t_it.push_back();
    }
    printf("%llu %d %d is\n", __rdtscp(&A) - st, maxSteps, maxStepIx);
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
