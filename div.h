#ifndef DIV_H
#define DIV_H

class DivLut {
  using Numerator = unsigned long;
  using Denominator = int;
  constexpr static int lg_TableSize = 8;
  constexpr static int TableSize = 1 << lg_TableSize;

public:
  class Divisor {
    Numerator p;

    public:
    constexpr Divisor() : p(0) {}
    constexpr Divisor(Numerator d) : // ceiling 2^k/d
      p(1 == d ? (~0) : ((1UL << 63) - 1 + d) / d * 2) { } // k-64 = 0
    friend Numerator operator/(__uint128_t n, const Divisor& d) { return n * d.p >> 64; }
    Divisor operator<<(int n) const {
      Divisor d = *this;
      d.p >>= n;
      return d;
    }
    Divisor operator*(Numerator n) const {
      Divisor d = *this;
      d.p *= n;
      return d;
    }
  };

private:
  std::array<Divisor, TableSize> t;
public:

  constexpr DivLut() { for (auto i = 1; i < TableSize; i++) t[i] = i; }
  Divisor operator[](Numerator d) const {
    if (d > TableSize) {
      const Denominator k = lgl_flr(d) - lg_TableSize;
      return t[d >> k] << k;
    } else {
      return t[d]; 
    }
  }
  DivLut& operator<<=(int n) { for (auto& d : t) d = d << n; return *this;}
  DivLut& operator*=(Numerator n) { for (auto& d : t) d = d * n; return *this;}
};
#endif
