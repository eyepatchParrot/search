#ifndef BIN_H
#define BIN_H

#include "lin.h"
#include "util.h"

#include <limits>

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
  std::vector<Key> v;
  int lg_v;

  auto szA() { return v.size() - 32; };
  const Key* a() { return &v[16]; };

  Key bsEq(const Key x) {
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
  auto bsLin(const Key x) {
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
    if (MIN_EQ_SZ == 1) return a[leftIndex];

    auto guess = leftIndex + n/2;
    if (a[guess] < x) return a[baseForwardSearch(a,guess+1,x)];
    else return a[baseBackwardSearch(a,guess,x)];
  }

  public:
  Binary(const std::vector<Key>& _a) : v(_a.size() + 32,0) {
    // put barriers to allow fast linear search
    std::copy(_a.begin(), _a.end(), v.begin() + 16);
    std::fill(v.begin(), v.begin() + 16, std::numeric_limits<Key>::min());
    std::fill(v.end()-16, v.end(), std::numeric_limits<Key>::max());
    //auto n = szA();
    lg_v=0;
    for (auto n = szA();n > MIN_EQ_SZ; n -= (n/2)) lg_v++;

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
  const char* name() {
    switch (f) {
      case BsFn::BS_EQ:
        return "binary-naive";
      case BsFn::BS_LIN:
        if (MIN_EQ_SZ == 1) {
          return "binary-size";
        } else {
          return "binary-linear";
        }
        break;
      default:
        return "???";
    }
  }
};

#endif
