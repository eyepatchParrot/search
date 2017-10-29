#ifndef PADDED_VECTOR_H
#define PADDED_VECTOR_H

#include <vector>
#include "util.h"
#include <limits>

template <int pad=32>
class PaddedVector {
  std::vector<Key> v;

  public:
  PaddedVector(const std::vector<Key>& v) : v(v.size() + 2*pad) {
    std::copy(v.begin(), v.end(), this->v.begin() + pad);
    std::fill(this->v.begin(), this->v.begin() + pad,
        std::numeric_limits<Key>::min());
    std::fill(this->v.end() - pad, this->v.end(),
        std::numeric_limits<Key>::max());
  }
  Key& operator[](long ix) {
    // allow some inaccuracy to reduce needed precision
    assert(ix >= -pad); assert(ix <= size() + pad);
    return v[ix+pad]; 
  }
  size_t size() { return v.size() - 2*pad; }
  Key back() { return (*this)[size()-1]; }
};

#endif
