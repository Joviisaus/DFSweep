#ifndef __OCTTREE_H__
#define __OCTTREE_H__
#include <Eigen/Dense>
#include <algorithm>
#include <limits>
#include <memory>
#include <omp.h>
#include <unordered_map>

class OctreeNode {
public:
  Eigen::Vector3f center;
  float halfSize;
  std::vector<int> pointIndices;
  std::vector<std::shared_ptr<OctreeNode>> children;
  bool isLeaf;

  OctreeNode(const Eigen::Vector3f &center, float halfSize)
      : center(center), halfSize(halfSize), isLeaf(true) {}
};
#endif
