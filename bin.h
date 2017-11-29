#ifndef BIN_H
#define BIN_H

#include <cinttypes>
#include <assert.h>

#include "lin.h"                                                                                          
#include "util.h"                                                                                         
#include "padded_vector.h"                                                                                

#include <limits>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

template <int MIN_EQ_SZ=1
          ,bool TEST_EQ=true
          ,bool FOR=false
          ,bool OVERFLOW_MATH=false
          ,typename Index = unsigned
          ,class Linear = LinearUnroll<>
          >
class BinaryLR {
  PaddedVector<> A;
  int lg_v;

  public:
  BinaryLR(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=0;
    for (auto n = A.size(); n > MIN_EQ_SZ; n -= (n/2)) lg_v++;
  }
    __attribute__((always_inline))
  Key operator()(const Key x) {
    // use pointer
    Index l = 0;
    Index r = A.size();

    for (int i =  0; FOR ? i < lg_v : r - l > 1; i++) {
      //IACA_START
      assert(l <= r);    // ordering check
      assert(l+r >= r); // overflow check
      Index m = OVERFLOW_MATH ? l + (r-l) / 2 : (l+r)/2;
      if (TEST_EQ) {
        if (A[m] < x) {
          l = m + 1;
        } else if (A[m] > x) {
          r = m;
        } else {
          l = r = m;
        }
      } else {
        if (A[m] <= x) l = m;
        else r = m;
      }
    }
    if (MIN_EQ_SZ == 1) {
      //IACA_END
      return A[l];
    }

    Index guess = (l+r)/2;
    if (A[guess] < x) return A[Linear::forward(A.begin(), guess+1, x)];
    else return A[Linear::reverse(A.begin(), guess, x)];
  }
};

template <int MIN_EQ_SZ = 32
,bool useN = true
,bool useFor = true
,typename Index = unsigned long
,class Linear = LinearUnroll<>
>
class BinarySize {
  PaddedVector<> A;
  int lg_v, lg_min;

  public:
  BinarySize(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=lg_min=0;
    for (auto n = A.size();n > 1; n -= (n/2)) {
      lg_v++;
      if (n > MIN_EQ_SZ) lg_min++;
    }
  }

  __attribute__((always_inline))
    Key operator()(const Key x) {
      Index n = A.size();
      auto a = A.begin();
      Index leftIndex = 0L;
      if (useFor) {
        if (useN) {
          Index mid = n - (1UL << (lg_v-1));
          leftIndex = a[mid] <= x ? mid : leftIndex;
            n -= mid;
            // wrong # of iterations?
            IACA_START
              if (MIN_EQ_SZ == 1) {
#pragma unroll(8)
                for (int i = 1; i < lg_min; i++) {
                  n/=2;
                  assert(n == (1 << (lg_v-1)) >> i);
                  leftIndex = a[leftIndex + n] <= x ? leftIndex + n : leftIndex;
                }
              } else {
#pragma unroll(4)
                for (int i = 1; i < lg_min; i++) {
                  n/=2;
                  assert(n == (1 << (lg_v-1)) >> i);
                  leftIndex = a[leftIndex + n] <= x ? leftIndex + n : leftIndex;
                }
              }
            IACA_END
          assert(n == (1 << (lg_v-lg_min)));
          assert(n == MIN_EQ_SZ);
        } else {
          IACA_START
          for (int i = 0; i < lg_v; i++) {
            Index half = n / 2;
            n -= half;
            Index mid = leftIndex+half;
            leftIndex = a[mid] <= x ? mid : leftIndex;
          }
          IACA_END
        }
      } else {
				if (useN) {
          Index mid = n - (1UL << (lg_v-1));
					leftIndex = a[mid] <= x ? mid : leftIndex;
          n -= mid;
					//n = (1UL << (lg_v-1));
				}
				while (n > MIN_EQ_SZ) {
          IACA_START
          Index half = n / 2;
          leftIndex = a[leftIndex + half] <= x ? leftIndex + half : leftIndex;
          if (useN) n = half;
					else n -= half;
        }
      }
      if (MIN_EQ_SZ == 1) return a[leftIndex];

      Index guess = leftIndex + n/2;
      if (a[guess] < x) return a[Linear::forward(a,guess+1,x)];
      else return a[Linear::reverse(a,guess,x)];
    }
};

using BinaryNaive = BinaryLR<1, true, false, true>;
using BinaryOver = BinaryLR<1, true, false>;
using BinaryNoEq = BinaryLR<1, false, false>;
using BinaryFor = BinaryLR<1, true, true>;
using BinaryForNoEq = BinaryLR<1, false, true>;
using BinaryLinLR = BinaryLR<32, false, true>;

using BinarySizeRecurse = BinarySize<1, false, false> ;
using BinarySizePow = BinarySize<1, true, false>;
//using BinarySizeFor = BinarySize<1, false, true>;
using BinarySizeForPow = BinarySize<1, true, true>;
using BinaryLinSize = BinarySize<32, true, true> ;

using B0 = BinaryLinSize;
using B1 = BinarySize<32, true, true, unsigned long, LinearUnroll<>>;

#endif
