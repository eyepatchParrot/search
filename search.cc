#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <time.h>
#include <tuple>
#include <vector>
#include <x86intrin.h>
#include <limits>
#include "omp.h"

#include "util.h"
#include "div.h"

typedef int64_t ll;

template <int pad=32>
class PaddedVector {
  std::vector<ll> v;

  public:
  PaddedVector(const std::vector<ll>& v) : v(v.size() + 2*pad) {
    std::copy(v.begin(), v.end(), this->v.begin() + pad);
    std::fill(this->v.begin(), this->v.begin() + pad,
        std::numeric_limits<ll>::min());
    std::fill(this->v.end() - pad, this->v.end(),
        std::numeric_limits<ll>::max());
  }
  ll& operator[](long ix) {
    // allow some inaccuracy to reduce needed precision
    assert(ix >= -pad); assert(ix <= size() + pad);
    return v[ix+pad]; 
  }
  size_t size() { return v.size() - 2*pad; }
  ll back() { return (*this)[size()-1]; }
};


struct TestStats {
  std::string name;
  std::vector<double> ns;
  bool ok;
};

template <class S>
TestStats benchmark(
    const std::vector<ll>& input,
    const std::vector<int>& indexes,
    S& search) {
  constexpr int nRuns = N_RUNS;
  auto ts = TestStats{search.name()};


  // vals is shuffled, so can't use it. Maybe shuffle indices and use that
  // next time
  std::vector<ll> vals(indexes.size());
  for (int i = 0; i < vals.size(); i++)
    vals[i] = input[indexes[i]];
  
  // get verification info
  auto expSum = 0UL, valSum = 0UL;
  for (auto j=0;j<nRuns;j++) for (auto i : indexes) expSum += vals[i];
  // https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector

  for (int runIx = 0; runIx < nRuns; runIx++) {
    double t0 = omp_get_wtime();

    for (int i = 0; i < vals.size(); i++) {
      //ERR("%d,%d,", i,indexes[i]);
      auto val = search(vals[i]);
      valSum += val;
      assert(val == vals[i]);
    }
    double t1 = omp_get_wtime();
    ts.ns.push_back(1e9*(t1-t0));
  }
  ts.ok = valSum == expSum;
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
    assert(i<1032);
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

enum BsFn {
  BS_EQ,
  BS_LIN
};
template <BsFn f
         ,int MIN_EQ_SZ = 16
         ,bool useFor = true
         ,SearchFn* baseForwardSearch = linSIMD
         ,SearchFn* baseBackwardSearch = linSIMD<true>
         >
class Binary {
  std::vector<ll> v;
  int lg_v;

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
//        return a[m];
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
    if (useFor) {
      for (int i = 0; i < lg_v; i++) {
        auto half = n / 2;
        n -= half;
        leftIndex = a[leftIndex + half] <= x ? leftIndex + half : leftIndex;
      }
    } else {
      while (n > MIN_EQ_SZ) {
        auto half = n / 2;
        n -= half;
        leftIndex = a[leftIndex + half] <= x ? leftIndex + half : leftIndex;
      }
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
    //auto n = szA();
    lg_v=0;
    for (auto n = szA();n > MIN_EQ_SZ; n -= (n/2)) lg_v++;

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
        ss << "binary-naive"; break;
      case BsFn::BS_LIN:
        if (MIN_EQ_SZ == 1) {
          ss << "binary-size";
        } else {
          ss << "binary-linear";
          //ss << "bsLin_" << MIN_EQ_SZ;
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
 ,IS_IDIV
 ,IS_FP
};
template <IsFn f = IS_LIN
         ,int nIter = 1
         ,bool fastFirst = true
         ,SearchFn* baseForwardSearch = linSIMD
         ,SearchFn* baseBackwardSearch = linSIMD<true>
         >
class Interpolation {
  static constexpr int pad = 32;
  int lgScale;
  DivLut::Divisor d_range_width;
  // describe all of the operations included
  DivLut tableDL;
  PaddedVector<> A;
  ll i_range_width;
  double f_aL, f_width_range;

  ll getMid(const ll x, const ll left, const ll right) {
    switch (f) {
      case IS_FP:
        return left + ((double)x - (double)(A[left])) * (double)(right-left) /
          (double)(A[right]-A[left]);
        break;
      case IS_IDIV:
        return left + (x-A[left]) / ((A[right]-A[left]) / (right-left));
        break;
      case IS_LIN:
        return left + (ll)(((x - A[left]) >> lgScale) * (right-left)) /
          tableDL[(A[right] - A[left]) >> lgScale];
        break;
    }
  }

  ll getMid(const ll x) {
    switch (f) {
      case IS_FP:
        return (ll)(((double)x - f_aL) * f_width_range);
        break;
      case IS_IDIV:
        return (x - A[0]) / i_range_width;
        break;
      case IS_LIN:
        return (uint64_t)(x - A[0]) / d_range_width;
        break;
    }
  }

  auto is(const ll x) {
    ll left = 0;
    ll right = A.size() - 1;
    assert(A.size() >= 1);

    ll mid = fastFirst? getMid(x) : getMid(x, left, right);
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      if (A[mid] < x) left = mid+1;
      else if (A[mid] > x) right = mid-1;
      else return A[mid];
      if (left==right) return A[left];
      
      assert(left<right);
      mid = getMid(x, left, right);
      if (nIter < 0) {
        if (mid >= right) return A[baseBackwardSearch(&A[0], right, x)];
        else if (mid <= left) return A[baseForwardSearch(&A[0], left, x)];
      }
      assert(mid >= left); assert(mid <= right);
    }

    if (A[mid] > x) return A[baseBackwardSearch(&A[0], mid - 1, x)];
    return A[baseForwardSearch(&A[0], mid, x)];
  }

  public:
  // can put lgScale internally if change indexing logic, doesn't help because
  // you'd need a subraction instead.
  // can avoid it with hoisted because don't need to do lookup
  // lookup pays for additional shift
  // I have to use my buffer space since my division is off by one
  Interpolation(const std::vector<ll>& v) : A(v) {
    // maybe want a.size() since we truncate
    lgScale = lg(A.size() - 1);
    //(tableDL <<= lgScale) *= (szA() - 1); 
    //tableDL <<= lgScale; 
    //auto d2 = tableDL[(_a.back() - _a.front()) >> lgScale];
    d_range_width = (DivLut::Divisor((A.back() - A[0]) >>  lgScale) << lgScale)
      / (A.size() - 1);
    i_range_width = (A.back() - A[0]) / (A.size() - 1);
    f_aL = A[0];
    f_width_range =  (double)(A.size() - 1) / (double)(A.back() - A[0]);
  }
  ll operator()(const ll x) {
    switch (f) {
      case IsFn::IS_LIN:
      case IsFn::IS_FP:
      case IsFn::IS_IDIV:
        return is(x);
        break;
    };
    return 0;
  }
  std::string name() {
    std::stringstream ss;
    switch (f) {
      case IsFn::IS_LIN:
        if (nIter < 0) ss << "isRecurse";
        else ss << "interpolation-linear";
        //else ss << "isLin_" << nIter << '_' << fastFirst;
        break;
      case IsFn::IS_IDIV:
        ss << "isIDiv" << '_' << fastFirst; break;
      case IsFn::IS_FP:
        if (nIter < 0) {
          if (!fastFirst) ss << "interpolation-naive";
          else ss << "interpolation-recurse";
        } else {
          ss << "interpolation-linear-fp";
          //else ss << "isFp" << '_' << fastFirst; break;
        }
        break;
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

// ./x fileName [benchmarks ...]
int main(int argc, char *argv[]) {
  using std::istream_iterator;
  constexpr int seed = 42;

  std::cerr << argv[1] << '\n';

  std::ifstream f(argv[1]);
  
  int nNums;
  f >> nNums;
  assert(nNums != 0);
  auto input = std::vector<ll>(istream_iterator<ll>(f), istream_iterator<ll>());
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
  for (int i = 2; i < argc; i++) {
    TestStats ts; 
    std::string s = argv[i];
    if (s == "binary-naive") ts = benchmark<Binary<BS_EQ>>(input,testIndexes);
    else if (s == "binary-size") ts = benchmark<Binary<BS_LIN,1>>(input,testIndexes);
    else if (s == "binary-linear") ts = benchmark<Binary<BS_LIN>>(input,testIndexes);
    else if (s == "interpolation-naive") ts = benchmark<Interpolation<IS_FP,-1,false,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-recurse") ts = benchmark<Interpolation<IS_FP,-1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-linear-fp") ts = benchmark<Interpolation<IS_FP,1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "interpolation-linear") ts = benchmark<Interpolation<>>(input, testIndexes);
    else if (s == "isRecurse") ts = benchmark<Interpolation<IS_LIN,-1,true,linUnroll,linUnroll<true>>>(input, testIndexes);
    else if (s == "isIDiv") ts = benchmark<Interpolation<IS_IDIV>>(input, testIndexes);
    else if (s == "isLin_1_slow") ts = benchmark<Interpolation<IS_LIN,1,false>>(input, testIndexes);
    else if (s == "isLin_2") ts = benchmark<Interpolation<IS_LIN,2>>(input, testIndexes);
    else if (s == "isSub") ts = benchmark<Interpolation<IS_LIN,1,true>>(input, testIndexes);
    else if (s == "oracle") ts = benchmark(input, testIndexes,o);
    if (!ts.ok)
      std::cerr << "mess up " << argv[1] << ' ' << s << '\n';
    tests.push_back(ts);
  }


  // Set-up the headers and organize data
  std::vector<std::vector<size_t>> testsIxs(tests.size(),
      std::vector<size_t>(N_RUNS));
  for (int i = 0; i < tests.size(); i++) {
    auto& t = tests[i];
    assert(t.cyclesByIx.size() == N_RUNS);
    printf("%s%s", 0 != i ? "," : "", t.name.c_str());
    std::iota(testsIxs[i].begin(), testsIxs[i].end(), 0);
#ifndef NSORT
    std::sort(testsIxs[i].begin(), testsIxs[i].end(), [t](size_t a, size_t b) {
        return t.ns[a] < t.ns[b]; });
#endif
  }
  for (int runIx = 0; runIx < N_RUNS && (N_SAMPLES < 0 ? true : runIx < N_SAMPLES); runIx++) {
    printf("\n");
    for (int testIx = 0; testIx < tests.size(); testIx++) {
      auto& t = tests[testIx];
      auto nsIx = testsIxs[testIx][runIx];
      printf("%s%.3f", 0 != testIx ? "," : "",
          t.ns[nsIx] / (double)testIndexes.size());
    }
  }
  printf("\n");
}
