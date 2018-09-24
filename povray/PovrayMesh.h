#ifndef POVRAYMESH_H
#define POVRAYMESH_H

#include <vector>
#include <tuple>
#include <unordered_map>

struct IntegerTriple {
  IntegerTriple(int i_, int j_, int k_) : i(i_), j(j_), k(k_) {}
  int i, j, k;
  friend std::ostream& operator<<(std::ostream& os, const IntegerTriple &t) { return os << t.i << ", " << t.j << ", " << t.k; }
};

template <typename T>
class PovrayMesh {
public:
  
  PovrayMesh() : _nextInt(0) {}
  void addPoint(const T point, const Vector3D& position) { _pointerToInt[point]=_nextInt; ++_nextInt; _positions.emplace_back(point, position); }
  void addTriangle(const T p1, const T p2, const T p3) { _triangles.emplace_back(p1, p2, p3); }
  
  std::vector<Vector3D> positions() const { 
    std::vector<Vector3D> pos;
    for(const auto &p : _positions) pos.emplace_back(p.second);
    return pos;
  }
  std::vector<IntegerTriple> triangles() const { 
    std::vector<IntegerTriple> tri;
    for(const auto &t : _triangles) tri.emplace_back(_pointerToInt.at(std::get<0>(t)), _pointerToInt.at(std::get<1>(t)), _pointerToInt.at(std::get<2>(t)));
    return tri;
  }

private:
  std::vector<std::pair<T, Vector3D> > _positions;
  std::vector< std::tuple<T, T, T> > _triangles;  
  std::unordered_map<T, int> _pointerToInt;
  int _nextInt;
};

#endif /* MESH3D_H */

