#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "div.h"
#include "lin.h"
#include "padded_vector.h"

#include <vector>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

class IBase {
public:
  using Index = int64_t;
  using PadVec = PaddedVector<>;

  static constexpr int Recurse = -1;
  static constexpr bool Precompute = true;

  template <bool saveFirst = false, bool sub = true>
  struct Lut {
    // maybe want a.size() since we truncate
    Lut(const PadVec& a) : A(a), lgScale(lg(A.size() - 1)), a0(A[0]) {
      if (sub) {
        d_range_width = (DivLut::Divisor((A.back() - A[0]) >>  lgScale) << lgScale) / (A.size() - 1);
      } else {
        d_range_width = (DivLut::Divisor(A.back() >> lgScale) << lgScale) / (A.size() - 1);
      }
    }

    const PadVec& A;
    int lgScale;
    DivLut::Divisor d_range_width;
    DivLut divisors;
    Key a0;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + (Key)(((x - A[left]) >> lgScale) * (right-left)) /
        divisors[(A[right] - A[left]) >> lgScale];
    }
    Index operator()(const Key x) {
      uint64_t d = sub? (x - (saveFirst? a0 : A[0])) : x;
      return d / d_range_width;
    }
  };
  template <bool precompute=false>
    struct Float {
      
      Float(const PadVec& a) : A(a), f_aL(A[0]),
      f_width_range( (double)(A.size() - 1) / (double)(A.back() - A[0])) {}

      const PadVec& A;
      const double f_aL;
      const double f_width_range;

      Index operator()(const Key x, const Index left, const Index right) {
        return left + ((double)x - (double)(A[left])) /
          (double)(A[right] - A[left]) * (double)(right-left);
      }

      Index operator()(const Key x) {
        return precompute? (Index)(((double)x - f_aL) * f_width_range) :
          (*this)(x, 0, A.size()-1);
      }
    };

  struct IntDiv {
    IntDiv(const PadVec& a) : A(a),
    i_range_width((A.back() - A[0]) / (A.size() - 1)) {}

    const PadVec& A;
    Key i_range_width;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + (x-A[left]) / ((A[right]-A[left]) / (right-left));
    }
    Index operator()(const Key x) {
      return (x - A[0]) / i_range_width;
    }
  };
protected:
  PadVec A;

  IBase(const std::vector<Key>& v) : A(v) {}
};


template <class Interpolate = IBase::Lut<>
         ,int nIter = 1
         ,class Linear = LinearUnroll<>
         ,int guardOff=0
         ,bool savePtr=false
         >
class Interpolation : public IBase {
  Interpolate interpolate;
  const Key* a0;

  __attribute__((always_inline))
  auto is(const Key x) {
    Index left = 0;
    Index right = A.size() - 1;
    assert(A.size() >= 1);
    auto a = savePtr? a0 : A.begin();

    Index mid = interpolate(x);
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      IACA_START
      if (a[mid] < x) left = mid+1;
      else if (a[mid] > x) right = mid-1;
      else return a[mid];
      if (left == right) return a[left];

      assert(left<right);
      mid = interpolate(x, left, right);
      if (nIter < 0) { 
        IACA_END
        if (mid+guardOff >= right) return a[Linear::reverse(a, right, x)];
        else if (mid-guardOff <= left) return a[Linear::forward(a, left, x)];
      }
      assert(mid >= left); assert(mid <= right);
    }

    if (a[mid] > x) {
    //IACA_END
      auto r = a[Linear::reverse(a, mid - 1, x)];
      return r;
    } else {
    //IACA_END
      auto r = a[Linear::forward(a, mid, x)];
      return r;
    }
  }

  public:
  Interpolation(const std::vector<Key>& v, const std::vector<int>& indexes) : IBase(v), interpolate(A), a0(A.cbegin()) { }

  __attribute__((always_inline))
  Key operator()(const Key x) { return is(x); }
};


using InterpolationNaive = Interpolation<IBase::Float<>,IBase::Recurse, LinearUnroll<>>;
using InterpolationPrecompute = Interpolation<IBase::Float<IBase::Precompute>,IBase::Recurse,LinearUnroll<>>;
using InterpolationLinearFp = Interpolation<IBase::Float<IBase::Precompute>>;
using InterpolationLinear = Interpolation<>;
using InterpolationLinearSave = Interpolation<IBase::Lut<true>>;
using InterpolationLinearSub = Interpolation<IBase::Lut<true, false>>;
//using InterpolationIDiv = Interpolation<IBase::IntDiv> ;
//using InterpolationLin_2 = Interpolation<IBase::Lut<>,2>;
//using B1 = Interpolation<IBase::Lut<true>, IBase::Recurse, LinearUnroll<>, 0, true> ;
//using B0 = Interpolation<IBase::Float<IBase::Precompute>, IBase::Recurse, LinearUnroll<>, 0, true> ;
#endif
