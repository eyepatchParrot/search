#ifndef BIN_H
#define BIN_H

#include <cinttypes>
#include <assert.h>

#include "lin.h"                                                                                          
#include "util.h"                                                                                         
#include "padded_vector.h"                                                                                

#include <limits>

template <int MIN_EQ_SZ=1
          ,bool TEST_EQ=true
          ,bool OVERFLOW_MATH=true
          ,bool FOR=false
          ,int BREAK=0>
class BinaryLR {
  PaddedVector<> A;
  int lg_v;

  public:
  BinaryLR(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=0;
    for (auto n = A.size();n > MIN_EQ_SZ; n -= (n/2)) lg_v++;
  }
    __attribute__((always_inline))
  Key operator()(const Key x) {
    unsigned l = 0;
    unsigned r = A.size();

    for (int i = BREAK==3? lg_v : 0; FOR ? (BREAK==3 ? 0 != i : i < lg_v) : r - l > 1; i = BREAK==3 ? i-1 : i+1) {
      assert(l <= r);    // ordering check
      assert(l+r >= r); // overflow check
      unsigned m = OVERFLOW_MATH ? l + (r-l) / 2 : (l+r)/2;
      if (TEST_EQ) {
        if (A[m] < x) {
          l = m + 1;
        } else if (A[m] > x) {
          r = m;
        } else {
          if (BREAK==4) return A[m];
          l = r = m;
          if (FOR && BREAK==1) i = lg_v;
          if (BREAK==2) break;
          if (FOR && BREAK==3) i = 1;
        }
      } else {
        if (A[m] <= x) l = m;
        else r = m;
      }
    }
    assert(A[l] == x);
    return A[l];
  }
};

template <int MIN_EQ_SZ = 16
,bool useFor = true
,class Linear = LinearSIMD<>
>
class BinarySize {
	PaddedVector<> A;
	int lg_v;

	public:
	BinarySize(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
		lg_v=0;
		for (auto n = A.size();n > MIN_EQ_SZ; n -= (n/2)) lg_v++;
	}

	__attribute__((always_inline))
		Key operator()(const Key x) {
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
};

using BinaryNaive = BinaryLR<>;
using BFor0 = BinaryLR<1, true, true, true, 0>;
using BFor1 = BinaryLR<1, true, true, true, 1>;
using BFor2 = BinaryLR<1, true, true, true, 2>;
using BFor3 = BinaryLR<1, true, true, true, 3>;
using BFor4 = BinaryLR<1, true, true, true, 4>;
using BinaryNaiveImm = BinaryLR<1,true,true>;
using BinarySizeRecurse = BinarySize<1> ;
using BinaryLinear = BinarySize<>;

#endif
