#ifndef DIV_H
#define DIV_H
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

#endif
