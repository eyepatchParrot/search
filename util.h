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

#if IS == 1
#define T1(X,Y,Z,W) X
#elif BS == 1
#define T1(X,Y,Z,W) Y
#elif OS == 1
#define T1(X,Y,Z,W) Z
#elif IS2 == 1
#define T1(X,Y,Z,W) W
#else
#define T1(X,Y,Z,W)
#endif
#if IS == 2
#define T2(X,Y,Z,W) X
#elif BS == 2
#define T2(X,Y,Z,W) Y
#elif OS == 2
#define T2(X,Y,Z,W) Z
#elif IS2 == 2
#define T2(X,Y,Z,W) W
#else
#define T2(X,Y,Z,W)
#endif
#if IS == 3
#define T3(X,Y,Z,W) X
#elif BS == 3
#define T3(X,Y,Z,W) Y
#elif OS == 3
#define T3(X,Y,Z,W) Z
#elif IS2 == 3
#define T3(X,Y,Z,W) W
#else
#define T3(X,Y,Z,W)
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
