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

int main() {
  const int FLUSH_SZ = 1 << 28;
  ll nNums;
  unsigned A;
  std::vector<ll> a1, a2, search;
  char* flush1 = new char[FLUSH_SZ];
  char* flush2 = new char[FLUSH_SZ];

  bool loaded = 1 == scanf("%llu", &nNums); (void)loaded; // silence not used
  assert(loaded);
  assert(nNums > 0);
  for (ll i = 0; i < nNums; i++) {
    ll x;
    loaded = 1 == scanf("%llu", &x);
    assert(loaded);
    a1.push_back(x);
    search.push_back(x);
  }

  std::srand(10);
  std::random_shuffle(search.begin(), search.end());

  printf("lower_bound,binSearch\n");

  // search in w.c. n, so search for everything should be better than n^2, so doable for 1e5 elements
  // don't forget to try cycles in cycles
  memcpy(flush1, flush2, FLUSH_SZ);
  ll st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    // can I beat STL bin search?
    std::lower_bound(a1.begin(), a1.end(), search[i]);
    __rdtscp(&A);
  }
  double t1 = (double)(__rdtscp(&A) - st);
  printf("%e", t1);

  memcpy(flush1, flush2, FLUSH_SZ);
  st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    binSearch(search[i], a1);
    __rdtscp(&A);
  }
  double t2 = (double)(__rdtscp(&A) - st);
  printf(",%e", t2);

  memcpy(flush1, flush2, FLUSH_SZ);
  ll mit = 1000000, mat = 0;
  // I have the feeling that the serialization is the issue here.
  st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    ll st2 = __rdtsc();
    binSearch(search[i], a1);
    ll et2 = __rdtscp(&A) - st2;
    mit = std::min(mit, et2);
    mat = std::max(mat, et2);
  }
  double t3 = (double)(__rdtscp(&A) - st);
  printf("\nindividual %e %e %e", t3, (double)mit, (double)mat);

  printf("\n");
  delete[] flush1;
  delete[] flush2;
}
