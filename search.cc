#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <immintrin.h>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <time.h>
#include <tuple>
#include <vector>
#include <x86intrin.h>
#include <limits>

#include "util.h"
#include "div.h"

typedef uint64_t ll;

struct TestStats {
  std::string name;
  std::vector<std::tuple<ll, int> > cyclesByIx;

  void datum(ll dt) { cyclesByIx.emplace_back(dt, cyclesByIx.size()); }
  void datum(ll st, ll et) { datum(et-st); }
  void data(TestStats& ts) {
      for (auto t_i : ts.cyclesByIx) datum(std::get<0>(t_i));
  }
};

template <class S>
TestStats benchmark(
    const std::vector<ll>& input,
    const std::vector<int>& indexes,
    S& search) {
  constexpr int nRuns = N_RUNS;
  auto ts = TestStats{search.name()};

  // get verification info
  auto expSum = 0UL, valSum = 0UL;
  for (auto j=0;j<nRuns;j++) for (auto i : indexes) expSum += input[i];

  // vals is shuffled, so can't use it. Maybe shuffle indices and use that
  // next time
  std::vector<ll> vals(indexes.size());
  for (int i = 0; i < vals.size(); i++)
    vals[i] = input[indexes[i]];


  for (int runIx = 0; runIx < nRuns; runIx++) {
    ll st = __rdtsc();

    for (int i = 0; i < vals.size(); i++) {
      ERR("%d,%d,", i,indexes[i]);
      auto val = search(vals[i]);
      valSum += val;
      assert(val == vals[i]);
    }

    ll et = __rdtsc();
    ts.datum(st, et);
  }
  if (valSum != expSum) // check even when profiling
    fprintf(stderr, "mess up %lu\n", vals.back());
  return ts;
}

template <class S>
TestStats benchmark(
    const std::vector<ll>& input,
    const std::vector<int>& indexes) {
  S s(input);
  return benchmark(input, indexes, s);
}

using SearchFn = int64_t(const ll*, int64_t, ll);
template < bool reverse=false, int roll=2>
int64_t linSIMD(const ll* arr, const int64_t guessIx, const ll x) {
  auto vecXX = reverse? _mm256_set1_epi64x(x): _mm256_set1_epi64x(x-1);
  auto ptr = arr;
  auto i = guessIx;
	auto misalignment = ((uintptr_t)(ptr+i) & 31)/sizeof(ll);
	for (int j = 0; j < 4*roll; j++)
    if (reverse? (arr[i-j] <= x) : arr[i+j] >= x) return reverse? i-j : i+j;
  i = reverse? (i-4*(roll-1) - misalignment) : i + 4*roll - misalignment;
  // 32-aligned main loop                                                          
  for (;;i = reverse?(i-16) : i+16) {
    assert(i<1000);
    assert(i>-32);
    auto sign = reverse?-1:1;
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i + sign*0));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + sign*4));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + sign*8));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + sign*12));
    auto cmp3 = reverse? _mm256_cmpgt_epi64(vecXX, av3) :  _mm256_cmpgt_epi64(av3, vecXX);
    auto msk3 = _mm256_movemask_epi8(cmp3);
    if (!msk3) continue;
    auto cmp0 = reverse? _mm256_cmpgt_epi64(vecXX, av0) :   _mm256_cmpgt_epi64(av0,vecXX );
    auto cmp1 = reverse? _mm256_cmpgt_epi64(vecXX, av1) :   _mm256_cmpgt_epi64(av1,vecXX );
    auto cmp2 = reverse? _mm256_cmpgt_epi64(vecXX, av2) :   _mm256_cmpgt_epi64(av2,vecXX );
    auto msk0 = _mm256_movemask_epi8(cmp0);
    auto msk1 = _mm256_movemask_epi8(cmp1);
    auto msk2 = _mm256_movemask_epi8(cmp2);
    if (msk0) return reverse? (i + 4 - _lzcnt_u32(msk0) / 8 - 0 * 4) : i + _tzcnt_u32(msk0) / 8 + 0 * 4;
    if (msk1) return reverse? (i + 4 - _lzcnt_u32(msk1) / 8 - 1 * 4) : i + _tzcnt_u32(msk1) / 8 + 1 * 4;
    if (msk2) return reverse? (i + 4 - _lzcnt_u32(msk2) / 8 - 2 * 4) : i + _tzcnt_u32(msk2) / 8 + 2 * 4;
    if (msk3) return reverse? (i + 4 - _lzcnt_u32(msk3) / 8 - 3 * 4) : i + _tzcnt_u32(msk3) / 8 + 3 * 4;
  }
}

template <bool reverse=false,int n=8>
int64_t linUnroll(const ll* a, int64_t m, ll k) {
  for (;;m = (reverse?m-n:m+n)) {
    for (int i = 0; i < n; i++) {
      assert(m+i < 1032); assert((m-i) > -32);
      if (reverse?(a[m-i]<=k):(a[m+i]>=k)) return reverse?(m-i):(m+i);
    }
  }
}

// use an enum to get search type
enum BsFn {
  BS_EQ,
  BS_LIN
};
template <BsFn f
         ,int MIN_EQ_SZ = 32
         ,SearchFn* baseForwardSearch = linUnroll
         ,SearchFn* baseBackwardSearch = linUnroll<true>
         >
class Binary {
  std::vector<ll> v;

  auto szA() { return v.size() - 32; };
  const ll* a() { return &v[16]; };

  ll bsEq(const ll x) {
    auto a = this->a();
    unsigned l = 0;
    unsigned r = szA();

    while (r - l > 1) {
      assert(l < r);    // ordering check
      assert(l+r >= r); // overflow check
      unsigned m = l + (r-l) / 2;
      if (a[m] < x) {
        l = m + 1;
      } else if (a[m] > x) {
        r = m;
      } else {
        l = r = m;
      }
    }
    assert(a[l] == x);
    return a[l];
  }

  // https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
  auto bsLin(const ll x) {
    auto n = szA();
    auto a = this->a();
    auto leftIndex = 0L;
    while (n > MIN_EQ_SZ) {
      auto half = n / 2;
      n -= half;
      leftIndex = a[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    }
    if constexpr(MIN_EQ_SZ == 1) return a[leftIndex];

    auto guess = leftIndex + n/2;
    if (a[guess] < x) return a[baseForwardSearch(a,guess+1,x)];
    else return a[baseBackwardSearch(a,guess,x)];
  }

  public:
  Binary(const std::vector<ll>& _a) : v(_a.size() + 32,0) {
    // put barriers to allow fast linear search
    std::copy(_a.begin(), _a.end(), v.begin() + 16);
    std::fill(v.begin(), v.begin() + 16, std::numeric_limits<ll>::min());
    std::fill(v.end()-16, v.end(), std::numeric_limits<ll>::max());
  }
  ll operator()(const ll x) {
    switch (f) {
      case BsFn::BS_EQ:
        return bsEq(x);
        break;
      case BsFn::BS_LIN:
        return bsLin(x);
        break;
    };
    return 0;
  }
  std::string name() {
    std::stringstream ss;
    switch (f) {
      case BsFn::BS_EQ:
        ss << "bsEq"; break;
      case BsFn::BS_LIN:
        if (MIN_EQ_SZ == 1) {
          ss << "bs";
        } else {
          ss << "bsLin_" << MIN_EQ_SZ;
        }
        break;
      default:
        ss << "???";
    }
    return ss.str();
  }
};

enum IsFn {
  IS_LIN
 ,IS_RECURSE
};
template <IsFn f = IS_LIN
         ,int nIter = 1
         ,SearchFn* baseForwardSearch = linSIMD
         ,SearchFn* baseBackwardSearch = linSIMD<true>
         >
struct Interpolation {
  static constexpr int pad = 32;
  int lgScale;
  DivLut::Divisor d;
  // describe all of the operations included
  DivLut tableDL;
  std::vector<ll> v;

  auto szA() { return v.size() - 2*pad; };
  const ll* a() { return &v[pad]; };

  auto isRecurse(const ll x) {
    ll l = 0;
    ll r = szA() - 1;
    assert(szA() - l >= 0); // assume non-empty vector
    auto a = this->a();

    auto m = l + (x-a[l]) / d;
    assert(m < szA()); assert(m >= l);
    assert(a[m] >= a[l]); // we know this because n would've been less than d

    for (;;) {
      if (a[m] < x) l = m + 1;
      else if (a[m] > x) r = m - 1;
      else return a[m];

      auto n = ((x-a[l]) >> lgScale) * (r-l);
      auto d = tableDL[(a[r] - a[l]) >> lgScale];
      m = l + n/d;
      assert(m <= szA()); assert(m >= l);
      if (m >= r) return a[baseBackwardSearch(a,r,x)];
      else if (m <= l) return a[baseForwardSearch(a,l,x)];
    }
  }

  auto isLin(const ll x) {
    ll l = 0;
    ll r = szA() - 1;

    const auto a = this->a();
    assert(r - l >= 0); // assume non-empty vector
    auto m = l + (x-a[l]) / d;

    for (int i = 1; i < nIter; i++) {
      assert(m <= r); assert(m >= l);
      if (a[m] <= x) l = m;
      if (a[m] >= x) r = m;
      // maybe test for equality
      auto n = ((x-a[l]) >> lgScale) * (r-l);
      auto d = tableDL[(a[r] - a[l]) >> lgScale];
      m = l + n/d;
      assert(m <= szA()); assert(m >= l);
    }

    if (a[m] > x) return a[baseBackwardSearch(a,m-1,x)];
    return a[baseForwardSearch(a,m,x)];
  }

  public:
  // can put lgScale internally if change indexing logic, doesn't help because
  // you'd need a subraction instead.
  // can avoid it with hoisted because don't need to do lookup
  // lookup pays for additional shift
  // I have to use my buffer space since my division is off by one
  Interpolation(const std::vector<ll>& _a) : v(_a.size() + 2*pad) {
    // put barriers to allow fast linear search
    std::copy(_a.begin(), _a.end(), v.begin() + pad);
    std::fill(v.begin(), v.begin() + pad, std::numeric_limits<ll>::min());
    std::fill(v.end()-pad, v.end(), std::numeric_limits<ll>::max());

    // maybe want a.size() since we truncate
    lgScale = lg(_a.size() -1);
    //(tableDL <<= lgScale) *= (szA() - 1); 
    //tableDL <<= lgScale; 
    //auto d2 = tableDL[(_a.back() - _a.front()) >> lgScale];
    d = (DivLut::Divisor((_a.back() - _a.front()) >>  lgScale) << lgScale) * (szA() - 1);
    //auto d3 = d2;
  }
  ll operator()(const ll x) {
    switch (f) {
      case IsFn::IS_LIN:
        return isLin(x);
        break;
      case IsFn::IS_RECURSE:
        return isRecurse(x);
        break;
    };
    return 0;
  }
  std::string name() {
    std::stringstream ss;
    switch (f) {
      case IsFn::IS_LIN:
        ss << "isLin_" << nIter; break;
      case IsFn::IS_RECURSE:
        ss << "isRecurse"; break;
      default:
        ss << "is???"; break;
    }
    return ss.str();
  }
};

class Oracle {
  const std::vector<ll>& a;
  std::vector<int> i;
  int j;

  public:
  // note that to ensure fast code, offset is always to the left, which means that offset searches wont have the full offset
  Oracle(const std::vector<ll>& _a, const std::vector<int>& _i, const int offset = 0, const bool rnd = false) : a(_a), j(0) {
    // easier to make a single long pipeline instead of a branch
    // There's blocking on the output to be produced, instead of being able to go through the steady state throughput
    // You may end up paying the latency cost multiple times instead of just once with left-deep
    // single relation plan because only one table in FROM clause
    for (int ti : _i) {
      int tti = ti;
      if (rnd && ((rand() % 2) == 0)) {
        tti += offset;
      } else {
        tti -= offset;
      }
      unsigned r = a.size() - 1;
      tti = tti < 0 ? 0 : tti;
      tti = tti > r ? r : tti;
      i.push_back(tti);
    }
  }

  ll operator()(const ll x) {
    int k = i[j++];
    j = j >= i.size() ? 0 : j;
    return a[k] == x ? a[k] : -1;
  }
  std::string name() { return "oracle"; }
};

// ./x [outputLength] [subsetSize]
int main(int argc, char *argv[]) {
  using std::cin; using std::istream_iterator;
  constexpr int seed = 42;
  
  int nNums;
  cin >> nNums;
  assert(nNums != 0);
  auto input = std::vector<ll>(istream_iterator<ll>(cin), istream_iterator<ll>());
  int nGets = SUBSET_SIZE < 1 ? input.size() : SUBSET_SIZE;

  // permute the items
  std::vector<int> testIndexes(nGets);
  {
    std::vector<int> allIndexes(input.size());
    std::iota(allIndexes.begin(), allIndexes.end(), 0);
    std::shuffle(allIndexes.begin(), allIndexes.end(), std::mt19937{seed});
    std::copy_n(allIndexes.begin(), nGets, testIndexes.begin());
  }
  Oracle o(input, testIndexes);

  std::vector<TestStats> tests;
  for (int i = 1; i < argc; i++) {
    TestStats ts; 
    std::string s = argv[i];
    if (s == "bsEq") ts = benchmark<Binary<BS_EQ>>(input,testIndexes);
    else if (s == "bs") ts = benchmark<Binary<BS_LIN,1>>(input,testIndexes);
    else if (s == "bsLin_32") ts = benchmark<Binary<BS_LIN,32>>(input,testIndexes);
    else if (s == "isRecurse") ts = benchmark<Interpolation<IS_RECURSE,-1,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "isLin_1") ts = benchmark<Interpolation<>>(input, testIndexes);
    else if (s == "isLin_2") ts = benchmark<Interpolation<IS_LIN,2>>(input, testIndexes);
    else if (s == "oracle") ts = benchmark(input, testIndexes,o);
    tests.push_back(ts);
  }

  // Set-up the headers and organize data
  bool first = true;
  for (auto& t : tests) {
    printf("%s%s", true != first ? "," : "", t.name.c_str());
#ifndef NSORT
    std::sort(t.cyclesByIx.begin(), t.cyclesByIx.end());
#endif
    first = false;
    assert(t.cyclesByIx.size() == N_RUNS);
  }
  for (int i = 0; i < N_RUNS && (N_SAMPLES < 0 ? true : i < N_SAMPLES); i++) {
    first = true;
    printf("\n");
    for (auto& ts : tests) {
      printf("%s%.3f", true != first ? "," : "", std::get<0>(ts.cyclesByIx[i]) / (double)testIndexes.size());
      first = false;
    }
  }
  printf("\n");
}
