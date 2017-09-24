#define NDEBUG

#include <cstdio>
#include <array>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <assert.h>
#include <string.h>

typedef uint64_t l;
typedef __uint128_t ll;

inline int lgl_flr(uint64_t x) {
  return 64 - __builtin_clzll(x);
}
inline unsigned lgl(uint64_t x) {
  return 64 - __builtin_clzll(x-1);
}

template <int LG_N>
struct DivLut {
  // Note: All indices are off by one to save a subtraction during division.
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut() {
    memset(pT, 0, sizeof(pT[0]) * N);
    memset(lg_qT, 0, sizeof(lg_qT[0]) * N);

    pT[1-1] = ~(T)0; // x * 1 ~= x * 2^64 - 1 / 2^64
    lg_qT[1-1] = 0;

    // do powers of 2
    // start at zero since we multiply by 2^64 / 2
    for (int i = 2, j = 0; i <= N; i*=2, j++) {
      pT[-1+i] = (T)1 << 63;
      lg_qT[-1+i] = j;
    }

    // generate fast multiply results
    for (int i = 1; i <= N; i++) {
      int d = i;
      int lg_q = 0;
      // reduce to an odd number
      while (pT[-1+d] == 0 && d % 2 == 0) {
        d /= 2;
        lg_q++;
      }

      // don't recalculate what we've already calculated
      if (pT[-1+d] != 0) {
        pT[-1+i] = pT[-1+d];
        lg_qT[-1+i] = lg_qT[-1+d] + lg_q;
        continue;
      }

      // find a 2^n + 1 multiple
      int j = 64;
      for (; j < 64+lgl(d); j++) {
        DT x = ((DT)1 << j) + 1;
        int r = x % d;
        if (r == 0) {
          pT[-1+i] = (((DT)1 << j) + d) / d;
          lg_qT[-1+i] = lg_q + (j - 64);
          break;
        }
      }
      if (pT[-1+i] == 0) {
        pT[-1+i] = (((T)1 << 63) / d) << 1;
        // guaranteed i > 1, so we will always have room to shift left
        // we always divide by 2^64, so we have to make sure that's happening
        lg_qT[-1+i] = lg_q;
      }
    }
  }

  T div(T n, T d) {
    d--; // start the d-1 early
    if (d > N-1) {
      // note that this becomes ceiling log
      const int k = lgl_flr(d) - LG_N;
      const T pIx = d >> k; // 1 + (d-1) >> k, but we started d, and -1 for index
      return divFit(n, pT[pIx], k + lg_qT[pIx]);
    } else {
      return divFit(n, pT[d], lg_qT[d]);
    }
  }

  void one_d(const T d, T& p, int& lg_q) const {
    if (d > N) {
      const int k = lgl_flr(d) - LG_N;
      T pIx = (d-1) >> k;
      p = pT[pIx];
      lg_q = k + lg_qT[pIx];
    } else {
      p = pT[d];
      lg_q = lg_qT[d];
    }
  }

  T divFit(T n, T p, int lg_q) {
    // assuming that we're in the range
    T hi = ((DT)n * p) >> 64;
    return hi >> lg_q;
  }

  T pT[N];
  int lg_qT[N];
};

class DivLut4 {
  using Numerator = unsigned long;
  using Denominator = int;
  using uint128_t = __uint128_t;
  constexpr static int lg_TableSize = 8;
  constexpr static int TableSize = 1 << lg_TableSize;
  class Divisor {
    public:
    Numerator p;
    Denominator lg_q;
    constexpr Divisor() : p(0), lg_q(0) {}
    constexpr Divisor(Numerator d) :
      p(1 == d ? (~0) : ((1UL << 63) + d - 1) / d * 2), // ceiling 2^k/d
      lg_q(0) { } // k-64 = 0
    friend Numerator operator/(uint128_t n, const Divisor& d) { return n * d.p >> 64 >> d.lg_q; }
    auto operator>>(int k) const {
      auto d = *this;
      d.lg_q += k;
      return d;
    }
  };

  std::array<Divisor, TableSize> t;
public:
  constexpr DivLut4() { for (auto i = 1; i < TableSize; i++) t[i] = i; }
  auto operator[](Numerator d) const {
    if (d > TableSize) {
      const Denominator k = lgl_flr(d) - lg_TableSize;
      return t[d >> k] >> k;
    } else {
      return t[d]; 
    }
  }
};

using namespace std;
int main() {
  const int LG_N = 8;
  const int N = 1 << LG_N;

  DivLut<8> dL;
  //DivLut3<8> dL3;
	DivLut4 dL3;

//#ifndef NDEBUG
  // TODO need to test on large numbers
  const int N_TEST = 1 << 10;
  int max2 = 0;
  int sum2 = 0;
  int cnt2 = 0;
  int max3 = 0;
  int sum3 = 0;
  int cnt3 = 0;

  for (int d = 1; d < N; d++) {
    for (l k = 0; k < N_TEST; k++) {
      l n = k;
      //l n = rand();
      l r = n / d;
      l r2 = dL.div(n, d);
      l r3 = n / dL3[d];
      int d2 = (int)r - (int)r2;
      max2 = max(max2, d2);
      sum2 += d2;
      if ((int)r - (int)r2 > 0) {
        cnt2 ++;
      }
      int d3 = (int)r - (int)r3;
      max3 = max(max3, d3);
      sum3 += d3;
      if ((int)r - (int)r3 > 0) {
        cnt3++;
      }
      if (d3 > d2) {
        auto dd = dL3[d];
        printf("n/d = %lu/%d = %lu ~ %lu = n * %lx / 2^%d\n", n,d, r, n / dd, dd.p, dd.lg_q);
      }
    }
  }
  printf("max d2 %d sum d2 %d cnt d2 %d\n", max2, sum2, cnt2);
  printf("max d3 %d sum d3 %d cnt d3 %d\n", max3, sum3, cnt3);
//#endif


  /*
  printf("p = { 0x%lx", dL3.pT[0] << dL2.lg_qT[0]);
  for (int i = 1; i < N; i++) {
    printf("UL , 0x%lx", dL3.pT[i] << dL2.lg_qT[i]);
  }
  printf("UL };\n");

  printf("p = { 0x%lx", dL2.pT[0]);
  for (int i = 1; i < N; i++) {
    printf("UL , 0x%lx", dL2.pT[i]);
  }
  printf("UL };\n");

  printf("q = { %d", dL2.lg_qT[0]);
  for (int i = 1; i < N; i++) {
    printf(", %d", dL2.lg_qT[i]);
  }
  printf(" };\n");
  */
  return 0;
}



