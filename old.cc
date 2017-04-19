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

