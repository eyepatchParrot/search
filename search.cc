#include <time.h>
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
  const int REPS = 1;
  ll nNums;
  std::vector<ll> array;
  std::vector<ll> search;

  assert(1 == scanf("%llu", &nNums));
  for (ll i = 0; i < nNums; i++) {
    ll x;
    assert(1 == scanf("%llu", &x));
    array.push_back(x);
    search.push_back(x);
  }

  struct timespec ts;

  // search in w.c. n, so search for everything should be better than n^2, so doable for 1e5 elements
  printf("Start Binary Search\n");
  // don't forget to try cycles in cycles
  std::lower_bound(array.begin(), array.end(), search[30]);
  ll mit = 1e18;
  ll mins = 1e18;
  ll mat = 0;
  ll mans = 0;
  for (ll i = 0; i < nNums; i++) {
    // can I beat STL bin search?
    unsigned A;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    ll sns = (ll)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    ll st = __rdtsc();
    for (int j = 0; j < REPS; j++) {
      std::lower_bound(array.begin(), array.end(), search[30]);
    }
    ll et = __rdtscp(&A);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    ll ens = (ll)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    mit = std::min(mit, et-st);
    mins = std::min(mins, ens-sns);
    mat = std::max(mat, et-st);
    mans = std::max(mans, ens-sns);
  }
  printf("End Binary Search\n");
  printf("min %llu max %llu diff %f\n", mit, mat, (mat - mit) / (double)mat);
  printf("min %llu max %llu diff %f\n", mins, mans, (mans - mins) / (double)mans);
//  double binSearchCycles = (double)(et - st) / nNums;
}
