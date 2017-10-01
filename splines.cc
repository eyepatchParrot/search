#include <functional>
#include <unistd.h>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "div.h"
using Ix = int;

template <typename T>
class Matrix {
  std::vector<T> v;
  const Ix W;

  public:
  Matrix(Ix h, Ix w) : v(w * h), W(w)  {}
  Matrix(Ix h, Ix w, T t) : v(w * h, t), W(w)  {}
  auto& at(Ix r, Ix c) { assert(r*W + c < v.size()); return v[r*W + c]; }
};

std::vector<long> readInts(std::string fileName) {
  std::ifstream f(fileName);
  std::vector<long> v;
  long n;
  f >> n;
  v.reserve(n);
  for (Ix i=0;i<n;i++) {
    long x;
    f >> x;
    v.push_back(x);
  }
  return v;
}

//Ix interpolate(long x, Ix n, long yL, long yR) {
//  return (x-yL) * n / (yR - yL);
//}

// add maxCost allowing you to reduce calculation time
float meanAbsErr(std::vector<long>& nums, Ix l, Ix r) {
  int lgScale = lgl(r-l);
  long absErr = 0;
  auto d = (DivLut::Divisor((nums[r] - nums[l]) >> lgScale) << lgScale) / (r-l);
  for (Ix y = 1+l; y < r; y++) {
    auto interpolate = l + (nums[y] - nums[l]) / d;
    auto diff = (long)y - (long)interpolate;
    absErr = diff < 0 ? absErr - diff : absErr + diff;
  }
  return 1.0 / (r - l) * absErr;
}

// use visited to record path
// consider floats and ints
std::vector<Ix> dijkstra(std::vector<long>& nums, long nLinks) {
  Ix s=0, t=nums.size()-1;
  Matrix<float> g(nums.size(), nums.size());
  for (Ix r = 0; r < nums.size(); r++)
    for (Ix c = r+2; c < nums.size(); c++)
      g.at(r,c) = meanAbsErr(nums, r, c);

  float maxCost = g.at(0, 999);
  std::vector<std::tuple<float, Ix,Ix,Ix>> q = {std::make_tuple(0.0, -1, s, -1) };
  Ix iter = 0,nDec=0,nSkip=0;
  Matrix<Ix> visited(nums.size(), nLinks, -1);
  while (q.size() > 0) {
    iter++;
    float dist;
    Ix len,cur,prev;
    std::tie(dist, len,cur,prev) = q[0];
    std::pop_heap(q.begin(), q.end(), std::greater());
    q.pop_back();

    if (cur == t) {
//      std::cout << "q=" << q.size() << " iter=" << iter << " cost=" << dist << " ndec=" << nDec << " nSkip=" << nSkip << '\n';
      std::vector<Ix> path = { cur };
      for (; len > 0; prev = visited.at(prev, --len)) {
        path.push_back(prev);
      }
      path.push_back(s);

      std::reverse(path.begin(), path.end());
      return path;
    }

    if (-1 != len) {
      if (-1 != visited.at(cur, len)) continue;
      visited.at(cur,len) = prev;
    }

    len++;

    assert(len <= nLinks);

    bool last = len == (nLinks - 1);
    for (Ix r = last? t : 1+cur; r < nums.size(); r++) {
      float cost = dist + g.at(cur,r);
      if (cost > maxCost) continue;
      if (r == t) {
        maxCost = cost;
        nDec++;
//        std::cout << std::count_if(q.begin(), q.end(), [cost](std::tuple<float,Ix,Ix,Ix> t){
//            return std::get<0>(t) > cost; }) << ' ' << q.size() << '\n';
//        std::cout << cost << ' ' << std::get<0>(q.back()) << '\n';
//        q.erase(std::remove_if(q.begin(), q.end(), [cost](std::tuple<float, Ix, Ix, Ix> t) {
//            return std::get<0>(t) > cost;}), q.end());
//        std::make_heap(q.begin(), q.end(), std::greater());
      }
      if (-1 != len && -1 != visited.at(r, len)) {
        nSkip++;
        continue;
      }
      q.emplace_back(dist + g.at(cur,r), len,r,cur);
      std::push_heap(q.begin(), q.end(), std::greater());
    }
  }

  return std::vector<Ix>();
}
int main(int argc, char* argv[]) {
  std::vector<int> nSplines;
  std::vector<std::string> fileNames;

  for (int opt = -1; -1 != (opt = getopt(argc, argv, "n:f:"));) {
    switch (opt) {
      case 'f':
        fileNames.push_back(optarg);
        break;
      case 'n':
        nSplines.push_back(std::stol(optarg));
        break;
      default:
        std::cerr << "./splines -f fileName [-n nSplines]" << std::endl;
        abort();
    }
  }

  std::cout << "file nSplines splines\n";
  for (auto fileName : fileNames) {
    auto nums = readInts(fileName);

    // dijkstras
    for (auto n : nSplines) {
      std::cout << fileName << ' ' << n << ' ';
      for (auto i : dijkstra(nums, n))
        std::cout << i << ',';
      std::cout << std::endl;
    }
  }

  return 0;
}






