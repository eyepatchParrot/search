#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "div.h"
#include "lin.h"
#include "padded_vector.h"

#include <vector>
#include <sstream>

enum IsFn {
  IS_LIN
 ,IS_IDIV
 ,IS_FP
};

class IBase {
protected:
  static constexpr int pad = 32;
  int lgScale;
  DivLut::Divisor d_range_width;
  // describe all of the operations included
  DivLut tableDL;
  PaddedVector<> A;
  Key i_range_width;
  double f_aL, f_width_range;

public:
  // can put lgScale internally if change indexing logic, doesn't help because
  // you'd need a subraction instead.
  // can avoid it with hoisted because don't need to do lookup
  // lookup pays for additional shift
  // I have to use my buffer space since my division is off by one
  IBase(const std::vector<Key>& v);
};

template <IsFn f = IS_LIN
         ,int nIter = 1
         ,bool fastFirst = true
         ,SearchFn* baseForwardSearch = linSIMD
         ,SearchFn* baseBackwardSearch = linSIMD<true>
         >
class Interpolation : public IBase {
  using Index = int64_t;

  Index getMid(const Key x, const Index left, const Index right) {
    switch (f) {
      case IS_FP:
        return left + ((double)x - (double)(A[left])) * (double)(right-left) /
          (double)(A[right]-A[left]);
        break;
      case IS_IDIV:
        return left + (x-A[left]) / ((A[right]-A[left]) / (right-left));
        break;
      case IS_LIN:
        return left + (Key)(((x - A[left]) >> lgScale) * (right-left)) /
          tableDL[(A[right] - A[left]) >> lgScale];
        break;
    }
  }

  Index getMid(const Key x) {
    switch (f) {
      case IS_FP:
        return (Index)(((double)x - f_aL) * f_width_range);
        break;
      case IS_IDIV:
        return (x - A[0]) / i_range_width;
        break;
      case IS_LIN:
        return (uint64_t)(x - A[0]) / d_range_width;
        break;
    }
  }

  auto is(const Key x) {
    Index left = 0;
    Index right = A.size() - 1;
    assert(A.size() >= 1);

    Index mid = fastFirst? getMid(x) : getMid(x, left, right);
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      if (A[mid] < x) left = mid+1;
      else if (A[mid] > x) right = mid-1;
      else return A[mid];
      if (left==right) return A[left];
      
      assert(left<right);
      mid = getMid(x, left, right);
      if (nIter < 0) {
        if (mid >= right) return A[baseBackwardSearch(&A[0], right, x)];
        else if (mid <= left) return A[baseForwardSearch(&A[0], left, x)];
      }
      assert(mid >= left); assert(mid <= right);
    }

    if (A[mid] > x) return A[baseBackwardSearch(&A[0], mid - 1, x)];
    return A[baseForwardSearch(&A[0], mid, x)];
  }

  public:
  using IBase::IBase;
  Key operator()(const Key x) {
    switch (f) {
      case IsFn::IS_LIN:
      case IsFn::IS_FP:
      case IsFn::IS_IDIV:
        return is(x);
        break;
    };
    return 0;
  }
  std::string name() {
    std::stringstream ss;
    switch (f) {
      case IsFn::IS_LIN:
        if (nIter < 0) ss << "isRecurse";
        else ss << "interpolation-linear";
        //else ss << "isLin_" << nIter << '_' << fastFirst;
        break;
      case IsFn::IS_IDIV:
        ss << "isIDiv" << '_' << fastFirst; break;
      case IsFn::IS_FP:
        if (nIter < 0) {
          if (!fastFirst) ss << "interpolation-naive";
          else ss << "interpolation-recurse";
        } else {
          ss << "interpolation-linear-fp";
          //else ss << "isFp" << '_' << fastFirst; break;
        }
        break;
      default:
        ss << "is???"; break;
    }
    return ss.str();
  }
};

#endif
