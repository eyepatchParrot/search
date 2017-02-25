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

int binSearchBlk(const ll k, const std::vector<ll>& a) {
  const int BLK = 8;
  assert(BLK % 2 == 0);

  unsigned l = 0, r = (BLK - 1 + a.size()) / BLK;
  int ix = -1;

  while (r - l > 1) {
    assert(l < r);
    assert(l+r >= r); // overflow check
    unsigned mBlk = (l+r) / 2;
    unsigned m = mBlk * BLK;
    unsigned m2 = std::min(m + BLK - 1, (unsigned)a.size());
    if (k < a[m]) {
      r = mBlk;
    } else if (k > a[m2]) {
      l = mBlk+1;
    } else {
      l = r = mBlk;
    }
  }

  l *= BLK;
  r = std::min((unsigned)a.size(), l + BLK);
  while (r - l > 1) {
    unsigned m = (l + r) / 2;
    if (k < a[m]) {
      r = m;
    } else if (k > a[m]) {
      l = m+1;
    } else {
      ix = m;
      break;
    }
  }
  return ix;
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
  const int FLUSH_SZ = 1 << 28;
  ll nNums;
  unsigned A;
  std::vector<ll> array;
  std::vector<std::tuple<ll, ll> > search;
  char* flush1 = new char[FLUSH_SZ];
  char* flush2 = new char[FLUSH_SZ];

  bool loaded = 1 == scanf("%llu", &nNums); (void)loaded; // silence not used
  assert(loaded);
  assert(nNums > 0);
  for (ll i = 0; i < nNums; i++) {
    ll x;
    loaded = 1 == scanf("%llu", &x);
    assert(loaded);
    array.push_back(x);
    search.push_back(std::tuple<ll, ll>(x, i));
    //search.push_back(std::make_tuple(x, i));
  }
  std::vector<ll> t_lb(nNums), t_bs(nNums);

  std::srand(10);
  std::random_shuffle(search.begin(), search.end());

  // search in w.c. n, so search for everything should be better than n^2, so doable for 1e5 elements
  // don't forget to try cycles in cycles
  memcpy(flush1, flush2, FLUSH_SZ);
  for (ll i = 0; i < nNums; i++) {
    ll st = __rdtsc();
    // can I beat STL bin search?
    std::lower_bound(array.begin(), array.end(), std::get<0>(search[i]));
    t_lb.push_back(__rdtscp(&A) - st);
  }

  memcpy(flush1, flush2, FLUSH_SZ);
  for (ll i = 0; i < nNums; i++) {
    ll st = __rdtsc();
    binSearch(std::get<0>(search[i]), array);
    t_bs.push_back(__rdtscp(&A) - st);
  }

  std::vector<ll> t_bk;
  memcpy(flush1, flush2, FLUSH_SZ);
  for (ll i = 0; i < nNums; i++) {
    ll st = __rdtsc();
    binSearchBlk(std::get<0>(search[i]), array);
    t_bk.push_back(__rdtscp(&A) - st);
  }

  std::vector<ll> t_it;
  memcpy(flush1, flush2, FLUSH_SZ);
  for (ll i = 0; i < nNums; i++) {
    ll st = __rdtsc();
    intSearch(std::get<0>(search[i]), array);
    t_it.push_back(__rdtscp(&A) - st);
  }

  std::vector<PerfStats> ps;
  ps.push_back(perfStats(t_lb, array, search, "lb"));
  ps.push_back(perfStats(t_bs, array, search, "bs"));
  ps.push_back(perfStats(t_bk, array, search, "bk"));
  ps.push_back(perfStats(t_it, array, search, "it"));

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

  delete[] flush1;
  delete[] flush2;
}
