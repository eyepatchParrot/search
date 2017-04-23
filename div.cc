#define NDEBUG

#include <cstdio>
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
struct DivLut2 {
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut2() {
    memset(pT, 0, sizeof(pT[0]) * N);
    memset(lg_qT, 0, sizeof(lg_qT[0]) * N);

    // I've moved everything back by one to accomodate common case

    // multiplying by one is approximated by (2^64 - 1) / 2^64
    pT[1-1] = ~(T)0;
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
      const T p = d >> k; // +1 for math, -1 for offset
      return divFit(n, p, k + lg_qT[p]);
    } else {
      return divFit(n, d+(-1+1), lg_qT[d+(-1+1)]);
    }
  }

  T divFit(T n, T p, int lg_q) {
    assert(p <= N);
    // assuming that we're in the range
    T hi = ((DT)n * pT[p]) >> 64;
    return hi >> lg_q;
  }

  T pT[N];
  int lg_qT[N];
};

template <int LG_N>
struct DivLut3 {
  typedef __uint128_t DT; // double T
  typedef uint64_t T;

  const static int N = 1 << LG_N;
  DivLut3() {
    memset(pT, 0, sizeof(pT[0]) * N);

    // counting on lg_q being small enough that we can shift right

    // I've moved everything back by one to accomodate common case

    // multiplying by one is approximated by (2^64 - 1) / 2^64
    pT[1-1] = ~(T)0;

    // do powers of 2
    // start at zero since we multiply by 2^64 / 2
    for (int i = 2, lg_q = 0; i <= N; i*=2, lg_q++) {
      T p = (T)1 << 63;
      p = p >> lg_q;
      pT[-1+i] = p;
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
      if (pT[-1+d] != 0) {
        if (i == d) continue;
        pT[-1+i] = pT[-1+d] >> lg_q;
        continue;
      }

      // find a 2^n + 1 multiple
      int j = 64;
      for (; j < 64+lgl(d); j++) {
        DT x = ((DT)1 << j) + 1;
        int r = x % d;
        if (r == 0) {
          lg_q = lg_q + (j-64);
          pT[-1+i] = ((((DT)1 << j) + d) / d) >> lg_q;
          break;
        }
      }
      if (pT[-1+i] == 0) {
        pT[-1+i] = ((((T)1 << 63) / d) << 1) >> lg_q;
      }
    }
  }

  T div(T n, T d) {
    d--; // start the d-1 early
    if (d > N-1) {
      // note that this becomes ceiling log
      const int k = lgl_flr(d) - LG_N;
      const T p = d >> k; // +1 for math, -1 for offset
      return divFit(n, p) >> k;
    } else {
      // -1 for index adjustment + 1 for early adjustment
      return divFit(n, d); 
    }
  }

  T divFit(T n, T p) {
    assert(p <= N);
    // assuming that we're in the range
    return ((DT)n * pT[p]) >> 64;
  }

  T pT[N];
};

template <int LG_W, int LG_L>
struct DivLut {
  const static ll W = 1ULL << LG_W;
  const static ll L = 1 + (1ULL << LG_L);
  const static ll OVERFLOW_CHK = 1ULL << (64 - LG_W);

  DivLut() : rNT(), lg_rDT() {
    /*
     * 1/d ~= 1 / r / 2^k ~= rN / rD / 2^k
     */
    // skip zero since divide by zero makes no sense
    for (ll r = 1; r < L; r++) {
      uint64_t rN = (1ULL << 63) / r; // 1/r ~= rN / 2^63

      int lz = __builtin_clzll(rN);
      lg_rDT[r] = lz + LG_W-1;
      rNT[r] = rN >> (63 - lg_rDT[r]);
    }
  }
  ll rNT[L];
  int lg_rDT[L]; // 2^8 * lg(rD)
  
  // approximate division
  ll div(const ll n, const ll d) const {
    ll p;
    int lg_q;
    if (d >= L) {
      one_dScale(d, p, lg_q);
    } else {
      one_dBase(d, p, lg_q);
    }
    
    if (n <= OVERFLOW_CHK) {
      return quotFit(n, p, lg_q);
    } else {
      return quotOverflow(n, p, lg_q);
    }
  }

  void one_dScale(const ll d, ll &p, int &lg_q) const {
    const int LG_D = lgl_flr(d);
    const int k = LG_D - LG_L;
    const ll r = 1 + ((d-1) >> k);
    p = rNT[r];
    lg_q = lg_rDT[r] + k;
  }

  void one_dBase(const ll d, ll &p, int &lg_q) const {
    p = rNT[d];
    lg_q = lg_rDT[d];
  }

  static inline ll quotFit(const ll n, const ll p, const int lg_q) {
    if (lg_q > 63) return 0;
    return (n * p) >> lg_q;
  }

  static inline ll quotOverflow(const ll n, const ll p, const int lg_q) {
    return ((n >> LG_W) * p) >> (lg_q - LG_W);
  }
};

using namespace std;
int main() {
  const int LG_N = 8;
  const int N = 1 << LG_N;

  DivLut<22, 8> dL;
  DivLut2<8> dL2;
  DivLut3<8> dL3;

#ifndef NDEBUG
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
      l r2 = dL2.divFit(n, d);
      l r3 = dL.div(n, d);
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
    }
  }
  printf("max d2 %d sum d2 %d cnt d2 %d\n", max2, sum2, cnt2);
  printf("max d3 %d sum d3 %d cnt d3 %d\n", max3, sum3, cnt3);
#endif


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
  return 0;
}



