//Search bs(const ll k, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  unsigned l = 0;
//  unsigned r = a.size();
//  int steps = 0;
//
//  while (r - l > 1) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    steps++;
//    if (a[m] < k) {
//      l = m + 1;
//    } else if (a[m] > k) {
//      r = m;
//    } else {
//      l = r = m;
//      //ret.ix = m;
//      //break;
//    }
//  }
//  assert(a[l] == k);
//  return Search{(int)l, steps};
//}

//Search bsNoEq(const ll k, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  unsigned l = 0;
//  unsigned r = a.size();
//  int steps = 0;
//
//  while (r - l > 1) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    steps++;
//    if (a[m] <= k) {
//      l = m;
//    } else if (a[m] > k) {
//      r = m;
//    }
//  }
//  assert(a[l] == k);
//  return Search{(int)l, steps};
//}

// PVK : https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
//Search bsPVK(const ll x, const std::vector<ll>& array, BinStruct& s) {
//  (void)s;
//  int leftIndex = 0;                                                               
//  int n = array.size();                                                            
//  int half;
//  int steps = 0;
//  if ((half = n) > 1) {
//    do {
//        half /= 2;
//        steps++;
//        n = array[leftIndex + half] == x ? 0 : n - half;
//        leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
//    } while ((half = n) > 1);
//  }                                                                                
//  assert(array[leftIndex] == x);  
//  return Search{leftIndex, steps};
//}

//template<int MIN_RANK_SZ>
//Search bsRank(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  int l = 0;
//  int n = a.size();
//  int steps = 0;
//  // try also doing it with half and quarter
//  int m;
//  while ((m = n / 4) >= MIN_RANK_SZ) {
//    int rank = (a[l+m] <= y) + (a[l+2*m] <= y) + (a[l+3*m] <= y);
//    steps += 3;
//    n = a[l+3*m] <= y ? n - 3 * m : m;
//    l += m * rank;
//  }
//  while ((m = n / 2) > 0) {
//    steps++;
//    n -= a[l+m] == y ? n : m;
//    l += a[l+m] <= y ? m : 0;
//  }
//  assert(a[l] == y);
//  return Search{l, steps};
//}

//template<int MIN_RANK_SZ, int K>
//Search bsRank2(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  (void)s;
//  const int MIN_M = MIN_RANK_SZ / K;
//  int l = 0;
//  int n = a.size();
//  int steps = 0;
//  // try also doing it with half and quarter
//  int m;
//  while ((m = n / K) >= MIN_M) {
//    int rank = a[l+m] <= y;
//    for (int i = 2; i < K; i++) {
//        rank += a[l+i*m] <= y;
//    }
//    steps += (K - 1);
//    n = a[l + (K - 1) * m] <= y ? n - (K-1) * m : m;
//    l += m * rank;
//  }
//  while ((m = n / 2) > 0) {
//    steps++;
//    n -= a[l+m] == y ? n : m;
//    l += a[l+m] <= y ? m : 0;
//  }
//  assert(a[l] == y);
//  return Search{l, steps};
//}

//Search bsArrBranchless(const ll y, const std::vector<ll>& a, BinStruct& s) {
//  // can we generate a branchless binary search?
//  // maybe use an integer to represent branching state. 
//  // could even put speculation in if there's enough depth
//  int steps = 0;
//  ll i[] = { 0, 0, a.size()};
//  while (i[2]-i[0] > 0) {
//    i[1] = (i[0] + i[2]) / 2;
//    ll x = a[i[1]] > i[1];
//    assert(x == 0 || x == 1);
//    i[0] = i[0+x];
//    i[2] = i[1+x];
//  }
//  if (a[i[0]] == y) {
//    return Search{(int)i[0], steps};
//  } else {
//    return Search{-1, steps};
//  }
//}

//Search leapSearch(const ll k, const std::vector<ll>& a, IntStruct& s) {
//  // 1. Interpolate
//  // 2. If the remaining bound is small, and not much progess is made, break;
//  // 3. Exponentially search for other bound.
//  // 4. Repeat
//  // 5. The bound is small and interpolation doesn't help much. Binary search.
//  //    No reason to go to interpolation search because bound is small.
//  // 
//  // can we change this to fix pt arithm?
//  // first priority is reducing number of lookups
//  // two ideas about systematically underestimating and overestimating
//  // 1. each repeated time you under-estimate, over-estimate by more the next time.
//  //    if you go on the other end, then reset and repeat.
//  // 2. Take the error by which you were off, and multiply your next estimate by that.
//  // 3. Look at euclid's alg and stuff like for for using error info to guide search
//  //
//  // have a notion of 'making progress'. If we interpolated and
//  // under-estimated, what is the distance from L to C? If we interpolated and
//  // over-estimated, what's the distance from R to C? If that distance is
//  // small, interpolation won't make good progress.
//  // 
//  // This can be extended to any cut. Interpolation is just one way of making
//  // a cut. We need to weigh the precision of attaining the other bound against
//  // the number of evaluations.
//  //
//  // If leap search hits half of the array, then we know that interpolation
//  // was extremely inaccurate. We should not repeat multiple bounded leap
//  // searches. 
//  //
//  // leaping expemplifies trust in the interpolation, so we can't measure
//  // progress.
//  //
//  // If we have to do log leaps, we shouldn't trust the interpolation, so do
//  // binary search.
//  int l = 0, r = a.size()-1;
//  Search ret = Search{-1, 0};
//  assert(r - l >= 0); // assume non-empty vector
//  const int MIN_SZ = 16;
//
//  while (r - l > 0) {
//    double p = 1.0 / (a[r] - a[l]) * (k - a[l]);
//    int m = l + (int)(p * (r-l));
//    ret.steps++;
//    // 4381798903807338309, 5921090195583789457, 6937039822968985763
//#ifndef NDEBUG
//    if (k == 0) {
//      printf("%d < %d < %d | %d i\n", l, m, r, r-l);
//    }
//#endif
//    if (k < a[m]) {
//      // over estimate
//      r = m - 1;
//      double progress = 1.0 - p; 
//      if (progress < 0.5) {
//        if ((r-l) <= MIN_SZ) break;
//        // potentially do x2
//        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)((1.0 - (1.0 / (a[r] - a[l]) * (k - a[l]))) * 2 * (r-l))));
//        int newR = r;
//        while (leap < (r-l)/2) {
//          m = r- leap;
//          ret.steps++;
//#ifndef NDEBUG
//          if (k == 0) printf("%d < %d < %d | %d li \n", l, m, newR, newR-l);
//#endif
//          if (k >= a[m]) {
//            break;
//          }
//          newR = m;
//          leap *=2;
//        }
//        if (leap >= (r-l)/2) {
//          r = newR;
//          break;
//        } else {
//          r = newR;
//          assert(m > l);
//          assert(m < r);
//          l = m;
//        }
//      }
//    } else if (k > a[m]) {
//      // under estimate
//      l = m + 1;
//      double progress = p;
//      if (progress < 0.5) {
//        if ((r-l) <= MIN_SZ) break;
//        //int leap = MIN_SZ;
//        int leap = std::min((r-l)/2, std::max(MIN_SZ, (int)(1.0 / (a[r] - a[l]) * (k - a[l]) * 2 * (r-l))));
//        int newL = l;
//        while (leap < (r-l)/2) {
//          m = l + leap;
//          ret.steps++;
//#ifndef NDEBUG
//          if (k == 0) printf("%d < %d < %d | %d ri \n", newL, m, r, r-newL);
//#endif
//          if (k <= a[m]) {
//            break;
//          }
//          newL = m;
//          leap *= 2;
//        }
//        if (leap >= (r-l)/2) {
//          l = newL;
//          break;
//        } else {
//          l = newL;
//          assert(m > l);
//          assert(m < r);
//          r = m;
//        }
//      }
//    } else {
//      ret.ix = m;
//      return ret;
//    }
//  }
//  while (r - l > 0) {
//    assert(l < r);    // ordering check
//    assert(l+r >= r); // overflow check
//    unsigned m = (l+r) / 2;
//    ret.steps++;
//#ifndef NDEBUG
//    if (k == 0) printf("%d < %d < %d | %d b \n", l, m, r, r-l);
//#endif
//    if (a[m] < k) {
//      l = m + 1;
//    } else if (a[m] > k) {
//      r = m;
//    } else {
//      ret.ix = m;
//      break;
//    }
//  }
//  if (k == a[l]) {
//    ret.ix = l;
//  }
//  return ret;
//}

// L + (R - L)(y - yL) / (yR - yL)
//Search isIntDiv(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  Search ret = Search{-1, 0};
//  ll l = 0, r = a.size() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yR = a[r], yL = a[l];
//  while (r - l > 0) {
//    assert(yR - yL > (1ULL << s.lgScale) && (y == yL || y - yL > (1ULL << s.lgScale)));
//    ll d = ((yR - yL) >> s.lgScale);
//    ll scOff = (r-l) * ((y - yL) >> s.lgScale) / d;
//    ll m = l + scOff;
//    assert(m <= r);
//    assert(m >= l);
//    ret.steps++;
//    if (y < a[m]) {
//      // over estimate
//      r = m - 1;
//      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//      yL = a[l];
//    } else {
//      ret.ix = m;
//      return ret;
//    }
//  }
//  if (y == a[l]) {
//    ret.ix = l;
//  }
//
//  return ret;
//}

// L + (R - L)(y - yL) / (yR - yL)
//Search is(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  ll l = 0, r = a.size() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yR = a[r], yL = a[l];
//  const unsigned lgScale = s.lgScale;
//  int steps = 0;
//  while (r - l > 0) {
//    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
//    ll n = (r-l)*((y-yL) >> lgScale);
//    ll d = ((yR - yL) >> lgScale);
//    ll scOff = dL.div(n,d);
//    ll m = l + scOff;
//    assert(m <= r);
//    assert(m >= l);
//    steps++;
//    if (y < a[m]) {
//      // over estimate
//      r = m - 1;
//      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//      yL = a[l];
//    } else {
//      return Search{(int)m, steps};
//    }
//  }
//
//  return Search{(int)l, steps};
//}

//Search is3(const ll y, const std::vector<ll>& a, IntStruct& s) {
//  ll l = 0, r = a.size() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yR = a[r], yL = a[l];
//  const unsigned lgScale = s.lgScale;
//  int steps = 0;
//  while (r - l > 0) {
//    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
//    ll n = (r-l)*((y-yL) >> lgScale);
//    ll d = (yR - yL) >> lgScale;
//    if (n < d) break; // if n < d, n / d == 0, so no progress made
//    //ll p;
//    //int lg_p, lg_q;
//    //dL2.one_d(d, p, lg_p, lg_q);
//
//    //ll n = ((y - yL) >> lgScale) * (r-l);
//    //if (n < d) break;
//    //ll scOff = dL2.timesFrac(n, p, lg_p, lg_q);
////
////    ll n = ((n >> lg_p) >> lgScale) * p * (r-l);
////    ll d = 
////
////    //ll n = (r-l)*((y-yL) >> lgScale);
////    ll n = y - yL;
////    if (n < 2) break; // could also work with n < d >> lgScale
////    int lg_n = lgl(n);
////
////    ll scOff = dL2.abc_d(n, lg_n, p, lg_p, r-l, lgScale, lg_q);
////    //if (n < d) break; // if n < d, n / d == 0, so no progress made
////
////    //steps++;
//    ll scOff = dL.div(n,d);
//    ll m = l + scOff;
//#ifndef NDEBUG
//    printf(" %llu", m);
//#endif
//    assert(m <= r);
//    assert(m >= l);
//    if (y < a[m]) {
//      // over estimate
//      r = m - 1;
//      yR = a[r];
//    } else if (y > a[m]) {
//      // under estimate
//      l = m + 1;
//      yL = a[l];
//    } else {
//      return Search{(int)m, steps};
//    }
//  }
//
//  while (a[l] < y && l < r) l++;
//  return Search{(int)l, steps};
//}

Search isDeriv(const ll y, const std::vector<ll>& a, IntStruct& s) {
  const int LG_W = 22;

  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  ll m2 = 0;
  const unsigned lgScale = s.lgScale;
  int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll d = ((yR - yL) >> lgScale);
    ll p;
    int lg_q;
    dL.one_d(d, p, lg_q);

    ll np = p*((y - yL) >> LG_W);
    ll n = (r-l)*(np >> lgScale);
    ll scOff = n >> (lg_q - LG_W);
    //ll n = (r-l)*(t1 >> lgScale);
    //if (n < d) break;
    //ll scOff = dL.div(n, d);
//
//
//    //steps++;
//    ll scOff = dL.quot(n,p, lg_q);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %llu,%llu", m, m2);
#endif
    assert(m <= r);
    assert(m >= l);
    if (y < a[m]) {
      // over estimate
      const int x = r - (m - 1);
      ll np2 = np * x;
      if (np2 < np) {
        printf(" OVERFLOW ");
      }
      np2 = np2 >> lg_q;

      r = m - 1;
      ll n2 = ((np + np2) >> lgScale)*(r-l);
      m2 = n2 >> (lg_q - LG_W);
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      return Search{(int)m, steps};
    }
  }

  while (a[l] < y && l < r) l++;
  return Search{(int)l, steps};
}

Search oracleLin(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  while (s.a[i] < y) i++;
  return Search{i, 0};
}

Search oracleLinRnd(const ll y, OracleStruct& s) {
  int i = s.i[s.j++];
  s.j = s.j >= s.i.size() ? 0 : s.j;
  // look at range of salaries compared to range of search as a reduction factor
  if (s.a[i] < y) {
    do { i++; } while (s.a[i] < y);
  } else {
    while (s.a[i] > y) { i--; }
  }
  return Search{i, 0};
}

Search isFull(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  //int steps = 0;
  while (r - l > 0) {
    assert(yR - yL > (1ULL << lgScale) && (y == yL || y - yL > (1ULL << lgScale)));
    ll n = (r-l)*((y-yL) >> lgScale);
    ll d = ((yR - yL) >> lgScale);
    if (n < d) break; // if n < d, n / d == 0, so no progress made

    //steps++;
    ll scOff = dL.div(n,d);
    ll m = l + scOff;
#ifndef NDEBUG
    printf(" %ld", m);
#endif
    assert(m <= r);
    assert(m >= l);
    if (y < a[m]) {
      // over estimate
      r = m - 1;
      yR = a[r];
    } else if (y > a[m]) {
      // under estimate
      l = m + 1;
      yL = a[l];
    } else {
      return Search{(int)m};
    }
  }

  while (a[l] < y && l < r) l++;
  return Search{(int)l};
}

Search bsUnroll(const ll x, const BinStruct& s) {
  const int MIN_EQ_SZ_LG = 9;
  const std::vector<ll>& array = s.a;
  int leftIndex = 0;                                                               
  int n = array.size();                                                            
  int half;
  if ((half = n) > (1 << MIN_EQ_SZ_LG)) {
    for (int i = 0; i < MIN_EQ_SZ_LG; i++) {
      half /= 2;
      n -= half;
      leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
    }
  }
  while ((half = n) > 2) {
    half /= 2;
    n -= half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  while ((half = n) > 1) {
    half /= 2;
    n = array[leftIndex + half] == x ? 0 : n - half;
    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
  }
  assert(array[leftIndex] == x);  
  return Search{leftIndex};
}

Search isSIMD(const ll y, const IntStruct& s) {
  const std::vector<ll>& a = s.a;
  ll l = 0, r = a.size() - 1;
  assert(r - l >= 0); // assume non-empty vector
  ll yR = a[r], yL = a[l];
  const unsigned lgScale = s.lgScale;
  ll n = (r-l)*((y-yL) >> lgScale);
  ll d = ((yR - yL) >> lgScale);
  ll m = l + dL.div(n,d);
  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]); // we know this because n would've been less than d
  if (a[m] > y) {
    do { m--; } while (a[m] > y);
    return Search{(int)m};
  }

  // n < d implies that we should start from the left
  // we know that l = m because we didn't go into the only path ewhere that's not true
  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
  typedef uint64_t v4u __attribute__ ((vector_size (32)));
//  int i = 0;
  while (m+8 < r) {
//    if (++i > 0) printf("\n*** %d ***\n", i);
    const ll* v = a.data() + m;
    // Compare a vector of the data with a vector of the searched for element
    // Gather together the results of each comparison into one byte per element
    // is mask ever not true?
    int mask1 = 
          _mm256_movemask_epi8(
            _mm256_cmpgt_epi64((v4u){y,y,y,y}, (v4u){v[0], v[1], v[2], v[3]}));
    //uint64_t mask2 = (uint64_t)
    //      _mm256_movemask_epi8(
    //        _mm256_cmpgt_epi64((v4u){y,y,y,y}, (v4u){v[4], v[5], v[6], v[7]})) << 32;
    //uint64_t mask3 = ~(mask2 | mask1);
 //   bool a = ~mask;
 //   bool b = !mask;
 //   assert(a==b);
//    printf("\n%x\n", mask);
//    mask = ~mask;
    int mask3 = mask1;
    if (mask3) return Search{(int)(m) + (__builtin_ctz(mask3) / 8)};
    m += 8;
  }
  if (y >= yR) return Search{(int)r};
  while (a[m] < y) m++;

  return Search{(int)m};
}

//auto bsLin(const ll x, BinStruct& s) {
//  const ll* array = s.a();
//  const int MIN_EQ_SZ = 32;
//  const int roll = 12;
//  auto leftIndex = 0;                                                               
//  auto n = s.szA();
//  while (n > MIN_EQ_SZ) {
//    auto half = n / 2;
//    n -= half;
//    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
//  }
//  auto guess = leftIndex + n/2;
//  if (array[guess] < x) {
//    guess++;
//    for (;;guess+=roll) for(int i=0;i<roll;i++)  if (array[guess +i] >= x)    return array[guess+i];
//  } else for (;;guess-=roll) for(int i=0;i<roll;i++) if (array[guess-i] <= x) return array[guess-i];
//}
//
//auto bsPVKEq2(const ll x, BinStruct& s) {
//  const ll* array = s.a();
//  const int MIN_EQ_SZ = 2;
//  long leftIndex = 0;                                                               
//  int n = s.szA();
//  int half;
//  while ((half = n) > MIN_EQ_SZ) {
//    half /= 2;
//    n -= half;
//    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
//  }
//  while ((half = n) > 1) {
//    half /= 2;
//    n = array[leftIndex + half] == x ? 0 : n - half;
//    leftIndex = array[leftIndex + half] <= x ? leftIndex + half : leftIndex;
//  }                                                                                
//  assert(array[leftIndex] == x);  
//  return array[leftIndex];
//}
//
//int is(const ll y, IntStruct& s) {
//  const ll* a = s.a();
//  ll l = 0, r = s.szA() - 1;
//  assert(r - l >= 0); // assume non-empty vector
//  ll yL = a[l];
//  const unsigned lgScale = s.lgScale;
//  ll n = (r-l)*((y-yL) >> lgScale);
//  ll m = l + dL.divFit(n,s.p, s.lg_q);
//
//  assert(m <= r);
//  assert(m >= l);
//  assert(a[m] >= a[l]); // we know this because n would've been less than d
//  if (a[m] > y) {
//    m--;
//    for (;; m -= 8)
//      for (ll i = 0; i < 8; i++)
//        if (a[m-i] <= y) return (int)(m - i);
//    return (int)m;
//  }
//
//  // n < d implies that we should start from the left
//  // we know that l = m because we didn't go into the only path ewhere that's not true
//  // note that (a[m] < y && a[m] < yR) was better than (a[m] < y && m < yR).
//  for (;; m += 8)
//    for (ll i = 0; i < 8; i++)
//      if (a[m+i] >= y) return (int)(m+i);
//  return (int)m;
//}

template <int roll>
int isLoop(const ll y, IntStruct& s) {
#define N (r-l)*((y-a[l]) >> lgScale)
  const ll* a = s.a();
  ll l = 0, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  const unsigned lgScale = s.lgScale;
  ll m = l + dL.divFit(N,s.p, s.lg_q);
  for (;;) {
    assert(m <= r); assert(m >= l); assert(a[m] >= a[l]); // we know this because n would've been less than d
    if (a[m] > y) {
      for (ll i = 0; i < roll; i++)
        if (a[m-i] <= y) return (int)(m - i);
      r = m;
    } else {
      for (ll i = 0; i < roll; i++)
        if (a[m+i] >= y) return (int)(m+i);
      l = m;
    }
    ll d = (a[r] - a[l]) >> lgScale;
    m = l + dL.div(N,d);
  }
#undef N
}

template <bool reverse=false,int n=12>
int64_t linUnroll2(const ll* a, int64_t m, ll k) {
  // does a check to see if it can skip the whole group
  for (;;m = (reverse?m-n:m+n)) {
    if (reverse?!(a[m-n+1]<=k):!(a[m+n-1]>=k)) continue;
    for (int i = 0; i < n; i++) {
      assert(m+i < 1032); assert((m-i) > -32);
      if (reverse?(a[m-i]<=k):(a[m+i]>=k)) return reverse?(m-i):(m+i);
    }
  }
}

// check for equality for the first k items before doing a LEQ check
template <int k, SearchFn* sF, SearchFn* sB, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps6(const ll y, PruneStruct& s) {
  auto a = s.a();
  ll m;
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  for (auto i=0;i<k;i++) if (a[m+i] == y) return a[m+i];
  m = sF(a,m+k,y);
  if (a[m] == y) return a[m];

  auto b = s.b();
  assert(s.szA()/s.szB() == 4);
  auto m2 = m >> 2;
  if (b[m2] > y) return b[baseBackwardSearch(b,m2-1,y)]; 
  return b[baseForwardSearch(b, m2, y)];
}

// use assignment to hoist out the checks
template <int k, SearchFn* sF, SearchFn* sB, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps7(const ll y, PruneStruct& s) {
  auto a = s.a();
  ll m;
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  for (auto i=0;i<k;i++) if (a[m+i] == y) return a[m+i];
  auto b = s.b();
  assert(s.szA()/s.szB() == 4);
  auto m2 = m >> 2;
  if (b[m2] > y) m2 = b[baseBackwardSearch(b,m2-1,y)]; 
  else m2 = b[baseForwardSearch(b, m2, y)];
  if (m2 == y) return m2;
  else return a[sF(a,m+k,y)];
}

auto ps(const ll y, PruneStruct& s) {
  auto a = s.a();
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  auto m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  if (a[m] > y) {
    while (m != 0 && a[m] > y) m--;
    if (a[m] == y) return a[m];
  } else {
    while (m < r && a[m] < y) m++;
    if (a[m] == y) return a[m];
  }
  auto b = s.b();
  auto i=0;
  for (;i<s.szB() && b[i] < y;i++) ;
  return b[i];
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps2(const ll y, PruneStruct& s) {
  auto a = s.a();
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  auto m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  if (a[m] > y) {
    m = baseBackwardSearch(a,m-1,y);
    if (a[m] == y) return a[m];
  } else {
    m = baseForwardSearch(a,m,y);
    if (a[m] == y) return a[m];
  }
  auto b = s.b();
  return b[baseForwardSearch(b, 0, y)];
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps3(const ll y, PruneStruct& s) {
  auto a = s.a();
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  auto m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  if (a[m] > y) {
    m = baseBackwardSearch(a,m-1,y);
    if (a[m] == y) return a[m];
  } else {
    m = baseForwardSearch(a,m,y);
    if (a[m] == y) return a[m];
  }
  auto b = s.b();
  auto m2 = s.szB() / 2;
  if (b[m2] > y) return b[baseBackwardSearch(b,m2-1,y)]; 
  return b[baseForwardSearch(b, m2, y)];
}

template <SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto ps4(const ll y, PruneStruct& s) {
  auto a = s.a();
  auto l = 0UL, r = s.szA() - 1;
  assert(r - l >= 0); // assume non-empty vector
  auto lgScale = s.lgScale;
  auto n = (r-l)*((y-a[l]) >> lgScale);
  auto m = l + dL.divFit(n,s.p, s.lg_q);

  assert(m <= r);
  assert(m >= l);
  assert(a[m] >= a[l]);

  if (a[m] > y) {
    m = linUnroll<true>(a,m-1,y);
    if (a[m] == y) return a[m];
  } else {
    m = linUnroll(a,m,y);
    if (a[m] == y) return a[m];
  }
  auto b = s.b();
  return b[binLin<512,baseForwardSearch,baseBackwardSearch>(b,0,s.szB(),y)];
}

template <int MIN_EQ_SZ, SearchFn* baseForwardSearch, SearchFn* baseBackwardSearch>
auto bsLinT(const ll x, BinStruct& s) {
  return s.a()[binLin<MIN_EQ_SZ, baseForwardSearch, baseBackwardSearch>(s.a(), 0, s.szA(), x)];
}

