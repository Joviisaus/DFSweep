#ifndef __DFCONTAINER_H__
#define __DFCONTAINER_H__
#include "OctTree.h"
#include "SweepDirDetector.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <float.h>

inline double epsilon = 1e-2f;
inline double PI = 3.1415926;

class DistanceField {
public:
  DistanceField();
  DistanceField(MeshLib::CTMesh *mesh);
  void SetMesh(MeshLib::CTMesh *mesh);
  void GridScalar(int MinScatter);
  void ComputeDistanceField();
  void readPrime(string primefile);
  std::vector<std::vector<std::vector<float>>> getField() { return Field; };
  std::vector<std::vector<std::vector<Eigen::Vector3f>>> getCoord() {
    return Coord;
  };

  std::vector<std::vector<std::vector<std::vector<float>>>>
  GetSweepProjScalar() {
    return this->SweepProjScalar;
  };
  std::vector<std::vector<std::vector<int>>> getGradianceCount() {
    return this->GradianceCount;
  };
  std::vector<std::vector<std::vector<float>>> getGradianceDiff() {
    return this->GradianceDiff;
  }
  void SaveFieldToBinary(const std::string &filename);
  void SaveGradianceToBinary(const std::string &filename);

protected:
  MeshLib::CTMesh *mesh;
  float PatchSize;
  std::vector<std::vector<float>> PointList;
  std::vector<std::vector<std::vector<float>>> Field;
  std::vector<std::vector<std::vector<int>>> FieldLabel;
  std::vector<std::vector<std::vector<float>>> GradianceDiff;
  std::vector<std::vector<std::vector<int>>> GradianceCount;

  std::vector<std::vector<std::vector<Eigen::Vector3f>>> GradianceField;
  std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord;
  std::vector<Eigen::Vector3f> SweepDir;
  std::vector<PrimeData> primes;
  std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjScalar;

  int maxPointsPerNode = 32;
  int maxDepth = 8;
  void BuildOctree();
  void BuildOctreeRecursive(std::shared_ptr<OctreeNode> node,
                            const std::vector<int> &pointIndices, int depth);
  void SweepProjection_Regist();
  void SubdivideNode(std::shared_ptr<OctreeNode> node);

  Eigen::Vector4f DistanceToMesh(int x, int y, int z);
  void FindNearestPointsInOctree(const Eigen::Vector3f &point,
                                 std::shared_ptr<OctreeNode> node,
                                 std::vector<int> &candidateIndices);
  void ExtractSweepDir();
  float PointToTriangleDistance(const Eigen::Vector3f &point,
                                const Eigen::Vector3f &v0,
                                const Eigen::Vector3f &v1,
                                const Eigen::Vector3f &v2);
  double DisCompute(Eigen::Vector3f point, int label);
  Eigen::Vector3f ClosestPointOnTriangle(const Eigen::Vector3f &point,
                                         const Eigen::Vector3f &v0,
                                         const Eigen::Vector3f &v1,
                                         const Eigen::Vector3f &v2);
  std::shared_ptr<OctreeNode> octreeRoot;
};

#endif
