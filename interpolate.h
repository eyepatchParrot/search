#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "div.h"
#include "lin.h"
#include "padded_vector.h"

#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <cmath>
#include <tuple>

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

  template <bool saveFirst = false, bool sub = true, bool fold=false>
  struct Lut {
    // maybe want a.size() since we truncate
    Lut(const PadVec& a) : A(a), lgScale(lg(A.size() - 1)), a0(A[0]) {
      if (fold) {
        divisors /= (A.size() - 1);
        //divisors <<= lgScale;
      }
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
      if (fold) return left + (Key)((x - A[left]) >> lgScale) / divisors[(A[right] - A[left])];
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
      assert(left >= 0); assert(right < A.size());
      mid = interpolate(x, left, right);
#ifndef NDEBUG
      auto fp = IBase::Float<IBase::Precompute>(A);
      auto d2 = (int)fp(x, left, right) - mid;
      if (d2*d2 > 1) {
        printf("%lu, %lu, %lu = %lu ~ %ld\n", x, left, right, mid, d2);
      }
#endif

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

// consider using AI techniques to find a good solution without trying to find
// an optimal solution. The connection from the points and slopes isn't clear
// with regards to how it impacts the solution. One thing to notice is that
// for each solution, we can consider how the distance from the interpolation
// changes over time and how it's different in the moment. However, since the
// goal is to hit a point or not at all, how to achieve that objective is less
// clear
class InterpolationSet : public IBase {
  static constexpr int vector_width = 8;
  using RiseRun = double;
  using KeyVector = std::array<Key, vector_width>;
  using RiseRunVector = std::array<RiseRun, vector_width>;

  struct Point {
    Key x;
    Index y;
  };

  struct Line {
    Index intercept;
    RiseRun slope;
    Line(Point a, Point b) {
      // want mx + b form
      // have (ax, ay), (bx, by)
      // m = (by - ay) / (bx - ax)
      // m * ax + b = ay
      // b = ay - m * ax
      // will never increase the floor function, but should help with float
      // arithmetic
      //auto epsilon = (b.x - a.x) / 2e19; 
      auto epsilon = 0.0;
      slope = (RiseRun)(b.y - a.y + epsilon) / (b.x - a.x);
      intercept = a.y - slope * a.x;
    }
    Line(Point p) : intercept(p.y), slope(0.0) { }
  };

  static Index interpolate(Line line, Key x) {
    return line.intercept + (Index)(line.slope * x);
  }

  static bool hit(Line line, Point point) {
    // note that once the interpolation exceeds the size of the set, since it
    // is non-decreasing, it will miss all remaining elements.
    return interpolate(line, point.x) == point.y;
  }

  KeyVector points;
  std::vector<RiseRunVector> slopes;


  int numMissed(const std::vector<Point>& set, int cap, Line line) {
    // return the number missed up to cap
    int missed = 0;
    for (int j = 0; missed < cap; j++) {
      if (j == set.size()) return std::min(cap, missed);
      missed += !hit(line, set[j]);
    }
    return missed;
  }

  Line removeBestCoveringLine(std::vector<Point>& set) {
    assert(!set.empty());

    Point left{.x=A[0], .y=0};

    // TODO refactor this so that making it a two dimensional loop is more
    // straightforward
    Line best_line(left);
    // evaluate the quality of each line
    for (auto [right, min_missed] = std::tuple{set.crbegin(), set.size()}; right != set.crend(); right++) {
      if (right->x == left.x) continue;
      Line line(left, *right); // should include epsilon
      int missed = numMissed(set, min_missed, line);
      if (missed < min_missed) {
        best_line = line;
        min_missed = missed;
      }
    }
    auto old_size = set.size();
    set.resize(std::distance(set.begin(),
        std::remove_if(set.begin(), set.end(), [best_line](auto point) {
        return hit(best_line, point); })));
    std::cout << (old_size - set.size()) << '\n';
    return best_line;
  }

  public:
    InterpolationSet(const std::vector<Key>& v, const std::vector<int>& indexes) : IBase(v) {
      for (auto& p : points) p = v.front();
      /* Want to build a minimal set of lines that covers the entire set.
       * First simplification is to always start from the leftmost point. Also
       * look into allowing for other left points, but this reduces the time
       * complexity.
       *
       * Greedily find each line by minimizing the number of missed points.
       * By making it a minimization problem rather than a maximization problem,
       * you can do early stopping.
       */
      std::vector<Point> uncovered;
      for (int i = 0; i < v.size(); i++) uncovered.emplace_back(Point{.x=v[i], .y=i});

      for (int vector_index = 0; !uncovered.empty();
          vector_index = (vector_index + 1) % vector_width) {
        if (vector_index == 0) slopes.emplace_back();
        Line l = removeBestCoveringLine(uncovered);
        slopes.back()[vector_index] = l.slope;
        // TODO add intercepts
      }
      std::cout << slopes.size() * vector_width << '\n';
      for (auto slope_v : slopes) for (auto slope : slope_v)
        std::cout << slope << '\n';
    }

    Key operator()(const Key x) { return x; }
};


using InterpolationNaive = Interpolation<IBase::Float<>,IBase::Recurse, LinearUnroll<>>;
using InterpolationPrecompute = Interpolation<IBase::Float<IBase::Precompute>,IBase::Recurse,LinearUnroll<>>;
using InterpolationRecurseGuard = Interpolation<IBase::Float<IBase::Precompute>, IBase::Recurse, LinearUnroll<>, 32, false>;
using InterpolationRecurse3 = Interpolation<IBase::Float<IBase::Precompute>, 3>;
using InterpolationRecurseLut = Interpolation<IBase::Lut<>, IBase::Recurse, LinearUnroll<>, 32, false>;
using InterpolationLinearFp = Interpolation<IBase::Float<IBase::Precompute>>;
using InterpolationLinear = Interpolation<>;
using i_simd = Interpolation<IBase::Lut<>, 1, LinearSIMD<>>;
using InterpolationLinearSave = Interpolation<IBase::Lut<true>>;
using InterpolationLinearSub = Interpolation<IBase::Lut<true, false>>;
//using InterpolationIDiv = Interpolation<IBase::IntDiv> ;
//using InterpolationLin_2 = Interpolation<IBase::Lut<>,2>;
using B1 = InterpolationRecurseLut;
using B0 = InterpolationRecurseGuard;
#endif
