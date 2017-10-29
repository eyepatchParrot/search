#ifndef BIN_H
#define BIN_H

#include "lin.h"
#include "util.h"
#include "padded_vector.h"

#include <limits>

enum BsFn {
  BS_EQ,
  BS_LIN
};
template <BsFn f
         ,int MIN_EQ_SZ = 16
         ,bool useFor = true
         ,class Linear = LinearSIMD<>
         >
class Binary {
  PaddedVector<> A;
  int lg_v;

  Key bsEq(const Key x) {
    unsigned l = 0;
    unsigned r = A.size();

    while (r - l > 1) {
      assert(l < r);    // ordering check
      assert(l+r >= r); // overflow check
      unsigned m = l + (r-l) / 2;
      if (A[m] < x) {
        l = m + 1;
      } else if (A[m] > x) {
        r = m;
      } else {
        l = r = m;
      }
    }
    assert(A[l] == x);
    return A[l];
  }

  // https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
  auto bsLin(const Key x) {
    auto n = A.size();
    auto a = A.begin();
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
    if (MIN_EQ_SZ == 1) return a[leftIndex];

    auto guess = leftIndex + n/2;
    if (a[guess] < x) return a[Linear::forward(a,guess+1,x)];
    else return a[Linear::reverse(a,guess,x)];
  }

  public:
  Binary(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=0;
    for (auto n = A.size();n > MIN_EQ_SZ; n -= (n/2)) lg_v++;
  }
  Key operator()(const Key x) {
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
};

using BinaryNaive = Binary<BS_EQ> ;
using BinarySize = Binary<BS_LIN,1> ;
using BinaryLinear = Binary<BS_LIN> ;

#endif
