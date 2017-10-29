#include "interpolate.h"

IBase::IBase(const std::vector<Key>& v) : A(v) {
    // maybe want a.size() since we truncate
    lgScale = lg(A.size() - 1);
    //(tableDL <<= lgScale) *= (szA() - 1); 
    //tableDL <<= lgScale; 
    //auto d2 = tableDL[(_a.back() - _a.front()) >> lgScale];
    d_range_width = (DivLut::Divisor((A.back() - A[0]) >>  lgScale) << lgScale)
      / (A.size() - 1);
    i_range_width = (A.back() - A[0]) / (A.size() - 1);
    f_aL = A[0];
    f_width_range =  (double)(A.size() - 1) / (double)(A.back() - A[0]);
}
