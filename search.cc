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

int binSearch(const ll k, const std::vector<ll>& a) {
  unsigned l = 0;
  unsigned r = a.size();
  int idx = -1;

  while (r - l > 0) {
    assert(l < r);    // ordering check
    assert(l+r >= r); // overflow check
    unsigned m = (l+r) / 2;
    if (a[m] == k) {
      idx = m;
      break;
    } else if (a[m] < k) {
      l = m + 1;
    } else {
      r = m;
    }
  }

  return idx;
}

int intSearch(const ll k, const std::vector<ll>& a) {
  // r is inclusive
  unsigned l = 0, r = a.size() - 1;
  ll vL = a[l], vR = a[r];
  int ix = -1;
  while (r - l >= 0) {
    assert(l <= r);    // ordering check
    double p = (double)(k-vL) / (double)(vR - vL);
    unsigned m = l + p * (r-l);
    if (k < a[m]) {
      r = m - 1;
      vR = a[r];
    } else if (k > a[m]) {
      l = m + 1;
      vL = a[l];
    } else {
      ix = m;
      break;
    }
  }
  return ix;
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

  {
    std::vector<ll> array(input);
    for (ll i = 0; i < nNums; i++) {
      ll st = __rdtsc();
      std::lower_bound(array.begin(), array.end(), std::get<0>(search[i]));
      t_lb.push_back(__rdtscp(&A) - st);
    }
  }

  {
    std::vector<ll> array(input);
    for (ll i = 0; i < nNums; i++) {
      ll st = __rdtsc();
      binSearch(std::get<0>(search[i]), array);
      t_bs.push_back(__rdtscp(&A) - st);
    }
  }

  std::vector<ll> t_it;
  {
    std::vector<ll> array(input);
    for (ll i = 0; i < nNums; i++) {
      ll st = __rdtsc();
      intSearch(std::get<0>(search[i]), array);
      t_it.push_back(__rdtscp(&A) - st);
    }
  }

  std::vector<PerfStats> ps;
  ps.push_back(perfStats(t_lb, input, search, "lb"));
  ps.push_back(perfStats(t_bs, input, search, "bs"));
  ps.push_back(perfStats(t_it, input, search, "it"));

  std::vector<std::string> fields {
    "name", "avg", "t90", "t99", "t999", "tMax", "v90", "v99", "v999", "vMax" };
  for (unsigned i = 0; i < fields.size(); i++) {
    printf("%s%s", i == 0 ? "" : ",", fields[i].c_str());
  }
  printf("\n");
  for (unsigned i = 0; i < ps.size(); i++) {
    printf("%s", ps[i].name.c_str());
    for (double t : ps[i].t) {
      printf(",%e", t);
    }
    for (ll v : ps[i].v) {
      printf(",%llu", v);
    }
    printf("\n");
  }
}
