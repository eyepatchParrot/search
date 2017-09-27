#ifndef UTIL_H
#define UTIL_H

#include <cinttypes>

#ifdef NDEBUG
#define ERR(A, args...)
#else
#define ERR(A, args...) fprintf(stderr, A, args)
#endif

#ifndef N_RUNS
#ifndef NDEBUG
#define N_RUNS 1 << 1
#else
#define N_RUNS 5000
#endif
#endif

#ifndef N_SAMPLES
#define N_SAMPLES 10
#endif

#ifndef SUBSET_SIZE
#define SUBSET_SIZE -1
#endif

unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
}

inline unsigned lgl(uint64_t x) {
  assert(x >= 2); // subtracting and clz < 1 undefined
  return 64 - __builtin_clzll(x-1);
}

inline int lgl_flr(uint64_t x) {
  assert(x >= 1); // clz < 1 undefined
  return 64 - __builtin_clzll(x);
}

// https://locklessinc.com/articles/sat_arithmetic/
inline std::uint64_t sub_sat_u64(std::uint64_t x, std::uint64_t y) {
  std::uint64_t res = x-y;
  res &= -(res <= x);
  return res;
}

#endif
