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

typedef unsigned long long ll;

int main() {
  const int FLUSH_SZ = 1 << 28;
  ll nNums;
  std::vector<ll> a1, a2;
  std::vector<ll> search;
  char* flush1 = new char[FLUSH_SZ];
  char* flush2 = new char[FLUSH_SZ];

  assert(1 == scanf("%llu", &nNums));
  for (ll i = 0; i < nNums; i++) {
    ll x;
    assert(1 == scanf("%llu", &x));
    a1.push_back(x);
    search.push_back(x);
  }

  std::srand(10);
  std::random_shuffle(search.begin(), search.end());

  // search in w.c. n, so search for everything should be better than n^2, so doable for 1e5 elements
  printf("Start Binary Search\n");
  // don't forget to try cycles in cycles
  unsigned A;
  memcpy(flush1, flush2, FLUSH_SZ);
  ll st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    // can I beat STL bin search?
    std::lower_bound(a1.begin(), a1.end(), search[i]);
  }
  ll t1 = __rdtscp(&A) - st;

  memcpy(flush1, flush2, FLUSH_SZ);
  st = __rdtsc();
  for (ll i = 0; i < nNums; i++) {
    std::lower_bound(a1.begin(), a1.end(), search[i]);
  }
  ll t2 = __rdtscp(&A) - st;
  printf("End Binary Search\n");
  printf("%e %e cycles to search every entry in array\n", (double)t1, (double)t2);
//  double binSearchCycles = (double)(et - st) / nNums;
  delete[] flush1;
  delete[] flush2;
}
