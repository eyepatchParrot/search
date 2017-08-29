#include <benchmark/benchmark.h>

#include <numeric>
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <immintrin.h>
#include <cinttypes>

#define REVERSE

typedef int64_t V;
typedef std::vector<V> vi;

using SearchFn = int(const V*, int, V);
template <SearchFn* f, bool reverse=false>
static void BM_Search(benchmark::State& state) {
  using std::cerr; using std::endl;
  const int vSz = 1000, kSz = 64, dist = state.range(0), seed = 42;
  vi v(vSz+64, std::numeric_limits<V>::max());
  auto l = v.begin()+32;
  auto p_l = &v[32];
  auto r = l + vSz;
  std::iota(l, r, 0);
  std::fill(v.begin(),l, std::numeric_limits<V>::min());
  vi v2(l + (reverse?0:dist), r - (reverse?dist:0));
  std::shuffle(v2.begin(), v2.end(), std::default_random_engine(seed));
  const vi k(v2.begin(), v2.begin() + kSz);
  for (int i = 0; i < kSz; i++) {
		if (k[i] != f(p_l,k[i] + (reverse? dist : -dist),k[i])) {
      f(p_l,k[i]+(reverse?dist:-dist),k[i]);
			exit(42);
    }
  }
	while (state.KeepRunning()) {
		for (int i = 0; i < kSz; i++)
			benchmark::DoNotOptimize(f(p_l, k[i] + (reverse?dist:-dist), k[i]));
	}
}
//

//#define R Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(7)->Arg(8)  \
//  ->Arg(9)->Arg(10)->Arg(11)->Arg(12)->Arg(13)->Arg(14)->Arg(15)->Arg(16)->Arg(17)->Arg(18)->Arg(19) \
//  ->Arg(20)->Arg(24)->Arg(28)->Arg(32)->Arg(36)->Arg(40)
//  ->Arg(44)->Arg(48)->Arg(52)->Arg(56)->Arg(60)->Arg(64)->Arg(68)->Arg(72)->Arg(76)
#define R RangeMultiplier(3)->Ranges({{2,40}})

template < bool reverse=false, int roll=1>
int linSIMD2(const V* arr, const int guessIx, const V x) {
  auto vecXX = reverse? _mm256_set1_epi64x(x): _mm256_set1_epi64x(x-1);
  auto ptr = arr;
  auto i = guessIx;
	int misalignment;
	for (int j = 0; j < 4*roll; j++) {
    if (j == 4*(roll-1)) misalignment = ((uintptr_t)(ptr+i) & 31)/sizeof(V);
    if (reverse? (arr[i-j] <= x) : arr[i+j] >= x) return reverse? i-j : i+j;
  }
  i = reverse? (i-4*(roll-1) - misalignment) : i + 4*roll - misalignment;
  // 32-aligned main loop                                                          
  for (;;i = reverse?(i-16) : i+16) {
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
BENCHMARK_TEMPLATE(BM_Search, linSIMD2)->R;
BENCHMARK_TEMPLATE(BM_Search, linSIMD2<false,3>)->R;
template < bool reverse=false, int roll=1>
int linSIMD(const V* arr, const int guessIx, const V x) {
  auto vecXX = reverse? _mm256_set1_epi64x(x): _mm256_set1_epi64x(x-1);
  auto ptr = arr;
  auto i = guessIx;
	int misalignment = ((uintptr_t)(ptr+i) & 31)/sizeof(V);
	for (int j = 0; j < 4*roll; j++)
    if (reverse? (arr[i-j] <= x) : arr[i+j] >= x) return reverse? i-j : i+j;
  i = reverse? (i-4*(roll-1) - misalignment) : i + 4*roll - misalignment;
  // 32-aligned main loop                                                          
  for (;;i = reverse?(i-16) : i+16) {
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
BENCHMARK_TEMPLATE(BM_Search, linSIMD)->R;
BENCHMARK_TEMPLATE(BM_Search, linSIMD<false,3>)->R;


template <int n,bool reverse=false>
int linUnroll(const V* a, int m, V k) {
  for (;;m += (reverse?-n:n)) {
    for (int i = 0; i < n; i++) {
      assert(m+i < 1032); assert(m-i > -32);
      if (reverse?(a[m-i]<=k):(a[m+i]>=k)) return reverse?(m-i):(m+i);
    }
  }
}

//BENCHMARK_TEMPLATE(BM_Search, linUnroll<8>)->R;
//BENCHMARK_TEMPLATE(BM_Search, linUnroll<4>)->R;
//BENCHMARK_TEMPLATE(BM_Search, linUnroll<8,true>,true)->R;


/*
int linSIMDa(const vi& arr, const int guessIx, const int x) {                      
  __m256i reorder = _mm256_set_epi32(7,3,6,2,5,1,4,0);                             
  auto vecXX = _mm256_set1_epi32(x - 1);                                           
  const int *ptr = arr.data();
	int i = guessIx;
	int misalignment = ((uintptr_t)(ptr+i) & 31)/4;
	for (int j = 0; j < 8; j++)
		if (arr[i+j] >= x) return i+j;
	i += 8 - misalignment;
  // 32-aligned main loop                                                          
  for (; ; i += 32) {
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + 8));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + 16));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + 24));
    auto cmp3 = _mm256_cmpgt_epi32(av3, vecXX);
    auto msk3 = _mm256_movemask_epi8(cmp3);
    if (!msk3) continue;
    auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    auto msk2 = _mm256_movemask_epi8(cmp2);
    if (!msk2) return __builtin_ctz(msk3) / 4 + i + 3 * 8;
    auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    auto msk1 = _mm256_movemask_epi8(cmp1);
    if (!msk1) return __builtin_ctz(msk2) / 4 + i + 2 * 8;
    auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX); // bail early?
    auto msk0 = _mm256_movemask_epi8(cmp0);
    if (!msk0) return __builtin_ctz(msk1) / 4 + i + 1 * 8;
    return __builtin_ctz(msk0) / 4 + i + 0 * 8;
    //auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX); // bail early?
    //auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    //auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    //auto cmp = (_mm256_packs_epi16(_mm256_packs_epi32(cmp0, cmp1), _mm256_packs_epi32(cmp2, cmp3)));
    //return i + __builtin_ctz(_mm256_movemask_epi8(_mm256_permutevar8x32_epi32(cmp,reorder)));
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMDa)->R;

long lin2(const vi& a, long m, V k) {
  if (a.back() < k) return a.size();
  for (;a[m]<k;m++) ;
  return m;
}

long lin3(const vi& a, long m, V k) {
  if (a.back() < k) return a.size();
  for (; m + 3 < a.size(); m += 4) {
    for (int i = 0; i < 4; i++) {
      if (a[m+i] >= k) return m + i;
    }
  }
  for (;a[m]<k;m++) ;
  return m;
}

long lin5(const vi& a, long m, const V k) {
  typedef V v4 __attribute__ ((vector_size (4*4)));
  if (a.back() < k) return a.size();
  v4 kv = {k,k,k,k};
  while (m + 3 < a.size()) {
    v4 av;
    for (long i = 0; i < 4; i++) {
      av[i] = a[m+i];
    }
    v4 bv = av < kv;
    v4 cv = {bv[2], bv[3], bv[0], bv[1]};
    bv += cv;
    for (long i = 0; i < 2; i++) {
        m -= bv[i];
    }
    if (a[m] >= k) break;
  }
  for (;a[m]<k;m++) ;
  return m;
}

long lin9(const vi& a, long m, const V k) {
  typedef V v4 __attribute__ ((vector_size (4*4)));
  v4 kv = {k,k,k,k};
  for (; m+3 < a.size(); m += 4) {
			const V* u = &a[m];
      v4 av = {u[0],u[1],u[2],u[3]};
      v4 bv = av >= kv;
      if (bv[0] | bv[1] | bv[2] | bv[3]) break;
  }
  while (a[m] < k) m++;
  return m;
}
//BENCHMARK_TEMPLATE(BM_Search, lin2)->RangeMultiplier(M)->Ranges(R);
//BENCHMARK_TEMPLATE(BM_Search, lin3)->RangeMultiplier(M)->Ranges(R);
//BENCHMARK_TEMPLATE(BM_Search, lin4)->RangeMultiplier(M)->Ranges(R);
BENCHMARK_TEMPLATE(BM_Search, lin9)->RangeMultiplier(M)->Ranges(R);

long lin2(const vi& a, long m, V k) {
  if (a.back() < k) return a.size();
  for (;a[m]<k;m++) ;
  return m;
}

long lin3(const vi& a, long m, V k) {
  if (a.back() < k) return a.size();
  for (; m + 3 < a.size(); m += 4) {
    for (int i = 0; i < 4; i++) {
      if (a[m+i] >= k) return m + i;
    }
  }
  for (;a[m]<k;m++) ;
  return m;
}

long lin4(const vi& a, long m, V k) {
  if (a.back() < k) return a.size();
  for (; m + 3 < a.size(); m += 4) {
    if (a[m+0] >= k || a[m+1] >= k || a[m+2] >= k || a[m+3] >= k)
        break;
  }
  for (;a[m]<k;m++) ;
  return m;
}

long lin5(const vi& a, long m, const V k) {
  typedef V v4 __attribute__ ((vector_size (4*4)));
  if (a.back() < k) return a.size();
  v4 kv = {k,k,k,k};
  long loops = (a.size() - m) / 4 * 4;
  for (;m + loops < a.size(); m++) if (a[m] >= k) return m;
  while (m + 3 < a.size()) {
    long j = 0;
    v4 av;
    for (long i = 0; i < 4; i++) {
      av[i] = a[m+i];
    }
    v4 bv = av < kv;
    v4 cv = {bv[2], bv[3], bv[0], bv[1]};
    bv += cv;
    for (long i = 0; i < 2; i++) {
        m -= bv[i];
    }
    if (a[m] >= k) break;
  }
  return m;
}



long lin8(const vi& a, long m, const V k) {
  for (; m+3 < a.size(); m += 4) {
      bool matches = false;
      for (int i = 0 ; i < 4; i++) {
          if (a[m+i] >= k) matches = true;
      }
      if (matches) break;
  }
  while (a[m] < k) m++;
  return m;
}
BENCHMARK_TEMPLATE(BM_Search, lin8)->RangeMultiplier(M)->Ranges(R);

long lin9(const vi& a, long m, const V k) {
  typedef V v4 __attribute__ ((vector_size (4*4)));
  v4 kv = {k,k,k,k};
  for (; m+3 < a.size(); m += 4) {
      const V* u = &a[m];
      v4 av = {u[0],u[1],u[2],u[3]};
      v4 bv = av >= kv;
      if (bv[0] | bv[1] | bv[2] | bv[3]) break;
  }
  while (a[m] < k) m++;
  return m;
}
BENCHMARK_TEMPLATE(BM_Search, lin9)->RangeMultiplier(M)->Ranges(R);

long lin6(const vi& a, long m, const V k) {
  while (a[m] < k) m++;
  return m;
}
BENCHMARK_TEMPLATE(BM_Search, lin6)->RangeMultiplier(M)->Ranges(R);

long lin7(const vi& a, long m, const V k) {
  for (; m+3 < a.size(); m += 4) {
      for (int i = 0 ; i < 4; i++) {
          if (a[m+i] >= k) return m+i;
      }
  }
  while (a[m] < k) m++;
  return m;
}
BENCHMARK_TEMPLATE(BM_Search, lin7)->RangeMultiplier(M)->Ranges(R);

template <int n>
long lina(const vi& a, long m, V k) {
  for (;; m += n) {
    for (int i = 0; i < n; i++)
        if (a[m+i] >= k) return m+i;
  }
  return m;
}

template <int n1, int n2>
long linb(const vi& a, long m, V k) {
  for (int i = 0; i < n1; i++)
    if (a[m+i] >= k) return m+i;
  for (;; m += n2)
    for (int i = 0; i < n2; i++)
        if (a[m+i] >= k) return m+i;
  return m;
}

template <int n>
long linc(const vi& a, long m, V k) {
  for (int i = 0; i < n; i++)
    if (a[m+i] >= k) return m+i;
  goto l2;
  for (;; m += 2*n) {
    for (int i = 0; i < n; i++)
        if (a[m+i] >= k) return m+i;
l2:
    for (int i = n; i < n+n; i++)
        if (a[m+i] >= k) return m+i;
  }
  return m;
}

long linSIMD(const vi& arr, const long guessIx, const V x) {
  using w8 = V __attribute__ ((vector_size (4*8)));
  using d4 = int64_t __attribute__ ((vector_size (4*8)));
  constexpr int roll = 8;
  constexpr union {
    int32_t i32[4];
    int64_t i64[2];
  } skip = {-4,-4,-4,-4};
  w8 xVec = {x,x,x,x,x,x,x,x};
  for (int i = guessIx;; i += roll) {
    w8 arrVec;
    for (long j = 0; j < roll; j++) arrVec[j] = arr[i+j];
    union {
        w8 i32;
        d4 i64;
    } cmpVec = {arrVec < xVec};
    w8 addVec;
    for (long j=0;j<roll;j++) addVec[j]=cmpVec.i32[7-j];
    cmpVec.i32 += addVec;
    for (long j=0;j<roll;j++) addVec[j]=cmpVec.i32[(j+2)%roll];
    cmpVec.i32 += addVec;
    if (cmpVec.i64[0] == skip.i64[0]) continue;
    return i - cmpVec.i32[0] - cmpVec.i32[1];
  }
}

BENCHMARK_TEMPLATE(BM_Search, linSIMD)->RangeMultiplier(M)->Ranges(R);
int linSIMD(const vi& arr, const int guessIx, const V x) {
  using v4 = V __attribute__ ((vector_size (4*4)));
  using dv2 = int64_t __attribute__ ((vector_size (4*4)));
  constexpr int roll = 4;
  constexpr union {
    int32_t i32[2];
    int64_t i64;
  } skip = {{-2,-2}};
  v4 xVec = {x,x,x,x};
  for (int i = guessIx;; i += roll) {
    v4 arrVec;
    for (int j = 0; j < 4; j++) arrVec[j] = arr[i+j];
    union {
        v4 i32;
        dv2 i64;
    } cmpVec = {arrVec < xVec};
    v4 cmpVec2 = {cmpVec.i32[3], cmpVec.i32[2], cmpVec.i32[1],cmpVec.i32[0]};
    cmpVec.i32 += cmpVec2;
    if (cmpVec.i64[0] == skip.i64) continue;
    return i - cmpVec.i32[0] - cmpVec.i32[1];
  }
}
int linSO(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x - 1);
  const int *ptr = arr.data();
  for (int i = guessIx; true; i += 4) {
    auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));
    auto cmp = _mm_cmpgt_epi32(arrVec, vecX);
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));
    if (mask)
      return i + __builtin_ctz(mask);
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSO)->RangeMultiplier(M)->Ranges(R);


int linSIMD3(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x - 1);
  const int *ptr = arr.data();
  for (int i = guessIx; true; i += 16) {
    auto av0 = _mm_loadu_si128((__m128i*)(ptr + i));
    auto av1 = _mm_loadu_si128((__m128i*)(ptr + i + 4));
    auto av2 = _mm_loadu_si128((__m128i*)(ptr + i + 8));
    auto av3 = _mm_loadu_si128((__m128i*)(ptr + i + 12));
    auto cmp0 = _mm_cmpgt_epi32(av0, vecX);
    auto cmp1 = _mm_cmpgt_epi32(av1, vecX);
    auto cmp2 = _mm_cmpgt_epi32(av2, vecX);
    auto cmp3 = _mm_cmpgt_epi32(av3, vecX);
    auto cmp = _mm_packs_epi16(_mm_packs_epi32(cmp0, cmp1), _mm_packs_epi32(cmp2, cmp3));
    int mask = _mm_movemask_epi8(cmp);
    if (mask)
      return i + __builtin_ctz(mask);
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD3)->RangeMultiplier(M)->Ranges(R);
BENCHMARK_TEMPLATE(BM_Search, linSIMD)->RangeMultiplier(M)->Ranges(R);
int linVP(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x - 1);
  constexpr int roll = 4;
  const int *ptr = arr.data();
  for (int i = guessIx;;i+=4*roll) {
    for (int j=0;j<roll;j++) {
      auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i + j*roll));
      auto cmp = _mm_cmpgt_epi32(arrVec, vecX);
      int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));
      if (mask)
        return i + j*roll + __builtin_ctz(mask);
    }
  }
}

//BENCHMARK_TEMPLATE(BM_Search, linUnroll<8>)->RangeMultiplier(M)->Ranges(R);
BENCHMARK_TEMPLATE(BM_Search, linVP)->RangeMultiplier(M)->Ranges(R);
template <int n>
int linUnroll(const vi& a, int m, V k) {
  for (;;m += n)
    for (int i = 0; i < n; i++)
        if (a[m+i] >= k) return m+i;
}
//BENCHMARK_TEMPLATE(BM_Search, linUnroll<32>)->RangeMultiplier(M)->Ranges(R);



int linSIMD5(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x);
  const int *ptr = arr.data();
  int i = guessIx;
  // unaligned start
  int misalignment = (uintptr_t)(ptr + i) & 15;
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));
  auto cmp = _mm_cmpgt_epi32(vecX, arrVec);
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp)) ^ 0xF;
  if (mask)
    return i + __builtin_ctz(mask);
  // continue with aligned part
  i += (16 - misalignment) / 4;
  for (; ; i += 16) {
    auto av0 = _mm_load_si128((__m128i*)(ptr + i));
    auto av1 = _mm_load_si128((__m128i*)(ptr + i + 4));
    auto av2 = _mm_load_si128((__m128i*)(ptr + i + 8));
    auto av3 = _mm_load_si128((__m128i*)(ptr + i + 12));
    auto cmp0 = _mm_cmpgt_epi32(vecX, av0);
    auto cmp1 = _mm_cmpgt_epi32(vecX, av1);
    auto cmp2 = _mm_cmpgt_epi32(vecX, av2);
    auto cmp3 = _mm_cmpgt_epi32(vecX, av3);
    auto cmp = _mm_packs_epi16(_mm_packs_epi32(cmp0, cmp1), _mm_packs_epi32(cmp2, cmp3));
    int mask = _mm_movemask_epi8(cmp);
    if (mask != 0xFFFF)
      return i + __builtin_ctz(~mask);
  }
}
static void BM_SIMD5(benchmark::State& state) {
  using std::cerr; using std::endl;
  const int vSz = 1000, kSz = 64, dist = state.range(0), seed = 42;
  vi v(vSz+4, std::numeric_limits<V>::max());
  std::iota(v.begin(), v.begin() + vSz, 0);
  vi v2(v.begin() + dist, v.begin() + vSz);
  std::shuffle(v2.begin(), v2.end(), std::default_random_engine(seed));
  const vi k(v2.begin(), v2.begin() + kSz);
  for (int i = 0; i < kSz; i++) { if (k[i] != linSIMD5(v,k[i]-dist,k[i])) exit(1);}
  while (state.KeepRunning()) {
    for (int i = 0; i < kSz; i++) {
      benchmark::DoNotOptimize(linSIMD5(v, k[i]-dist, k[i]));
    }
  }
}
BENCHMARK(BM_SIMD5)->RangeMultiplier(M)->Ranges(R);

int linSIMD9(const vi& arr, const int guessIx, const int x) {                      
  auto vecX = _mm_set1_epi32(x - 1);                                               
  auto vecXX = _mm256_set1_epi32(x - 1);                                           
  const int *ptr = arr.data();                                                     
  int i = guessIx;                                                                 
  // unaligned start                                                               
  int misalignment = (uintptr_t)(ptr + i) & 15;                                    
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));                              
  auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                        
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                               
  if (mask)                                                                        
    return i + __builtin_ctz(mask);                                                
  // continue with 16-aligned part                                                 
  i += (16 - misalignment) / 4;                                                    
  // if 32-unaligned, do an other 16-byte block                                    
  if ((uintptr_t)(ptr + i) & 31) {                                                 
    auto arrVec = _mm_load_si128((__m128i*)(ptr + i));                            
    auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                      
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                             
    if (mask)                                                                      
      return i + __builtin_ctz(mask);                                              
    i += 4;                                                                        
  }                                                                                
  // 32-aligned main loop                                                          
  for (; ; i += 32) {
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + 8));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + 16));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + 24));
    auto cmp3 = _mm256_cmpgt_epi32(av3, vecXX);
    auto msk3 = _mm256_movemask_epi8(cmp3);
    if (!msk3) continue;
    auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX); // bail early?
    auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    auto msk0 = _mm256_movemask_epi8(cmp0);
    auto msk1 = _mm256_movemask_epi8(cmp1);
    auto msk2 = _mm256_movemask_epi8(cmp2);
    if (msk0) return i + __builtin_ctz(msk0) / 4;
    if (msk1) return i + __builtin_ctz(msk1) / 4 + 1 * 8;
    if (msk2) return i + __builtin_ctz(msk2) / 4 + 2 * 8;
    if (msk3) return i + __builtin_ctz(msk3) / 4 + 3 * 8;
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD9)->RangeMultiplier(M)->Ranges(R);
int linSIMD8(const vi& arr, const int guessIx, const int x) {                      
  auto vecX = _mm_set1_epi32(x - 1);                                               
  auto vecXX = _mm256_set1_epi32(x - 1);                                           
  const int *ptr = arr.data();                                                     
  int i = guessIx;                                                                 
  // unaligned start                                                               
  int misalignment = (uintptr_t)(ptr + i) & 15;                                    
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));                              
  auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                        
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                               
  if (mask)                                                                        
    return i + __builtin_ctz(mask);                                                
  // continue with 16-aligned part                                                 
  i += (16 - misalignment) / 4;                                                    
  // if 32-unaligned, do an other 16-byte block                                    
  if ((uintptr_t)(ptr + i) & 31) {                                                 
    auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));                            
    auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                      
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                             
    if (mask)                                                                      
      return i + __builtin_ctz(mask);                                              
    i += 4;                                                                        
  }                                                                                
  // 32-aligned main loop                                                          
  for (; ; i += 32) {
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + 8));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + 16));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + 24));
    auto cmp3 = _mm256_cmpgt_epi32(av3, vecXX);
    auto msk3 = _mm256_movemask_epi8(cmp3);
    if (!msk3) continue;
    auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    auto msk2 = _mm256_movemask_epi8(cmp2);
    if (!msk2) return __builtin_ctz(msk3) / 4 + i + 3 * 8;
    auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    auto msk1 = _mm256_movemask_epi8(cmp1);
    if (!msk1) return __builtin_ctz(msk2) / 4 + i + 2 * 8;
    auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX); // bail early?
    auto msk0 = _mm256_movemask_epi8(cmp0);
    if (!msk0) return __builtin_ctz(msk1) / 4 + i + 1 * 8;
    return __builtin_ctz(msk0) / 4 + i + 0 * 8;
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD8)->RangeMultiplier(M)->Ranges(R);

int linSIMD7(const vi& arr, const int guessIx, const int x) {                      
  // 0,1,2,3,8,9,a,b,00,11,22,33,88,99,aa,bb,4,5,6,7,c,d,e,f,44,55,66,77,cc,dd,ee,ff
  // 0-------1-------2-----------3-----------4-------5-------6-----------7----------
  // 0,4,1,5,2,6,3,7                                                               
  __m256i reorder = _mm256_set_epi32(7,3,6,2,5,1,4,0);                             
  auto vecX = _mm_set1_epi32(x - 1);                                               
  auto vecXX = _mm256_set1_epi32(x - 1);                                           
  const int *ptr = arr.data();                                                     
  int i = guessIx;                                                                 
  // unaligned start                                                               
  int misalignment = (uintptr_t)(ptr + i) & 15;                                    
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));                              
  auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                        
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                               
  if (mask)                                                                        
    return i + __builtin_ctz(mask);                                                
  // continue with 16-aligned part                                                 
  i += (16 - misalignment) / 4;                                                    
  // if 32-unaligned, do an other 16-byte block                                    
  if ((uintptr_t)(ptr + i) & 31) {                                                 
    auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));                            
    auto cmp = _mm_cmpgt_epi32(arrVec, vecX);                                      
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));                             
    if (mask)                                                                      
      return i + __builtin_ctz(mask);                                              
    i += 4;                                                                        
  }                                                                                
  // 32-aligned main loop                                                          
  for (; ; i += 32) {
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + 8));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + 16));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + 24));
    auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX);
    auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    auto cmp3 = _mm256_cmpgt_epi32(av3, vecXX);
    auto cmp = (_mm256_packs_epi16(_mm256_packs_epi32(cmp0, cmp1), _mm256_packs_epi32(cmp2, cmp3)));
    if (_mm256_movemask_epi8(cmp))
      return i + __builtin_ctz(_mm256_movemask_epi8(_mm256_permutevar8x32_epi32(cmp,reorder)));
  }
}

BENCHMARK_TEMPLATE(BM_Search, linSIMD7)->RangeMultiplier(M)->Ranges(R);
*/
/*

int linSIMD6(const vi& arr, const int guessIx, const int x) {
  // 0,1,2,3,8,9,a,b,00,11,22,33,88,99,aa,bb,4,5,6,7,c,d,e,f,44,55,66,77,cc,dd,ee,ff
  // 0-------1-------2-----------3-----------4-------5-------6-----------7----------
	// 0,1,2,3,4,5,6,7 0,7 0,7
  //   0,1,2,3,4,5,6 1   2
  //     0,1,2,3,4,5 2   4
	//       0,1,2,3,4 3   6
  // 1,2,3,4,5,6,7   6   5
	// 2,3,4,5,6,7     5   3
  // 3,4,5,6,7       4   1
  // 0,4,1,5,2,6,3,7
  __m256i reorder = _mm256_set_epi32(7,3,6,2,5,1,4,0);
  auto vecX = _mm_set1_epi32(x - 1);
  auto vecXX = _mm256_set1_epi32(x - 1);
  const int *ptr = arr.data();
  int i = guessIx;
  // unaligned start
  int misalignment = (uintptr_t)(ptr + i) & 15;
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));
  auto cmp = _mm_cmpgt_epi32(arrVec, vecX);
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));
  if (mask)
    return i + __builtin_ctz(mask);
  // continue with 16-aligned part
  i += (16 - misalignment) / 4;
  // if 32-unaligned, do an other 16-byte block
  if ((uintptr_t)(ptr + i) & 31) {
    auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));
    auto cmp = _mm_cmpgt_epi32(arrVec, vecX);
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));
    if (mask)
      return i + __builtin_ctz(mask);
    i += 4;
  }
  // 32-aligned main loop
  for (; ; i += 32) {
    auto av0 = _mm256_load_si256((__m256i*)(ptr + i));
    auto av1 = _mm256_load_si256((__m256i*)(ptr + i + 8));
    auto av2 = _mm256_load_si256((__m256i*)(ptr + i + 16));
    auto av3 = _mm256_load_si256((__m256i*)(ptr + i + 24));
    auto cmp0 = _mm256_cmpgt_epi32(av0, vecXX);
    auto cmp1 = _mm256_cmpgt_epi32(av1, vecXX);
    auto cmp2 = _mm256_cmpgt_epi32(av2, vecXX);
    auto cmp3 = _mm256_cmpgt_epi32(av3, vecXX);
    auto cmp = _mm256_permutevar8x32_epi32(_mm256_packs_epi16(_mm256_packs_epi32(cmp0, cmp1), _mm256_packs_epi32(cmp2, cmp3)), reorder);
    int mask = _mm256_movemask_epi8(cmp);
    if (mask)
      return i + __builtin_ctz(mask);
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD6)->RangeMultiplier(M)->Ranges(R);
int linSIMD4(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x - 1);
  const int *ptr = arr.data();
  int i = guessIx;
  // unaligned start
  int misalignment = (uintptr_t)(ptr + i) & 15;
  auto arrVec = _mm_loadu_si128((__m128i*)(ptr + i));
  auto cmp = _mm_cmpgt_epi32(arrVec, vecX);
  int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));
  if (mask)
    return i + __builtin_ctz(mask);
  // continue with aligned part
  i += (16 - misalignment) / 4;
  for (; ; i += 16) {
    auto av0 = _mm_load_si128((__m128i*)(ptr + i));
    auto av1 = _mm_load_si128((__m128i*)(ptr + i + 4));
    auto av2 = _mm_load_si128((__m128i*)(ptr + i + 8));
    auto av3 = _mm_load_si128((__m128i*)(ptr + i + 12));
    auto cmp0 = _mm_cmpgt_epi32(av0, vecX);
    auto cmp1 = _mm_cmpgt_epi32(av1, vecX);
    auto cmp2 = _mm_cmpgt_epi32(av2, vecX);
    auto cmp3 = _mm_cmpgt_epi32(av3, vecX);
    auto cmp = _mm_packs_epi16(_mm_packs_epi32(cmp0, cmp1), _mm_packs_epi32(cmp2, cmp3));
    int mask = _mm_movemask_epi8(cmp);
    if (mask)
      return i + __builtin_ctz(mask);
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD4)->RangeMultiplier(M)->Ranges(R);

int linSIMD5(const vi& arr, const int guessIx, const int x) {
  auto vecX = _mm_set1_epi32(x - 1);
  const int *ptr = arr.data();
	int i = guessIx;
	int misalignment = ((uintptr_t)(ptr+i) & 15)/4;
	for (int j = 0; j < 5; j++)
		if (arr[i+j] >= x) return i+j;
	i += 4 - misalignment;
  for (; ; i += 16) {
    auto av0 = _mm_load_si128((__m128i*)(ptr + i));
    auto av1 = _mm_load_si128((__m128i*)(ptr + i + 4));
    auto av2 = _mm_load_si128((__m128i*)(ptr + i + 8));
    auto av3 = _mm_load_si128((__m128i*)(ptr + i + 12));
    auto cmp0 = _mm_cmpgt_epi32(av0, vecX);
    auto cmp1 = _mm_cmpgt_epi32(av1, vecX);
    auto cmp2 = _mm_cmpgt_epi32(av2, vecX);
    auto cmp3 = _mm_cmpgt_epi32(av3, vecX);
    auto cmp = _mm_packs_epi16(_mm_packs_epi32(cmp0, cmp1), _mm_packs_epi32(cmp2, cmp3));
    int mask = _mm_movemask_epi8(cmp);
    if (mask)
      return i + __builtin_ctz(mask);
  }
}
BENCHMARK_TEMPLATE(BM_Search, linSIMD5)->RangeMultiplier(M)->Ranges(R);
*/

BENCHMARK_MAIN();
