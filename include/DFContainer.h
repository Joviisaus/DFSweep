#ifndef __DFCONTAINER_H__
#define __DFCONTAINER_H__
#include "ColorImplementer.h"
#include "CuttingBox.h"
#include "CylinderCuttingBox.h"
#include "DevelopableReducer.h"
#include "MeshCutter.h"
#include "OctTree.h"
#include "SweepBlock.h"
#include "SweepHexMesher.h"
#include "SweepDirDetector.h"
#include "SweepDirFilter.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <float.h>
#include <unordered_map>

inline double epsilon = 1e-2f;
inline double PI = 3.1415926;
inline int SampleSize = 100;
inline float Alpha = 0.6;

void ComputeNearestPointsCPU(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<float>> &PointList,
    std::vector<std::vector<std::vector<int>>> &NearestIndex);

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
  std::vector<std::vector<std::vector<std::vector<float>>>>
  GetSweepProjEnergy() {
    return this->SweepProjEnergy;
  };
  const std::vector<std::string> &getSweepEnergyNames() const {
    return sweepEnergyNames;
  }
  std::vector<std::vector<std::vector<int>>> getGradianceCount() {
    return this->GradianceCount;
  };
  std::vector<std::vector<std::vector<bool>>> ForbiddenBoundaryPoints;
  std::vector<std::vector<std::vector<float>>> getGradianceDiff() {
    return this->GradianceDiff;
  }

  void exportPlanesToFile(const std::string &filename);
  std::vector<std::map<int, Eigen::Vector3f>> getCuttingHex() {
    return this->CuttingHexLists;
  };

  std::vector<Eigen::Vector3f> getSweepDir() { return this->SweepDir; }
  const std::vector<bool> &getSweepBlockNonPlanar() const {
    return sweepBlockNonPlanar;
  }
  const std::vector<Eigen::Vector3f> &getSweepBlockColors() const {
    return sweepBlockColors;
  }
  std::vector<int> getDisplayHexIndices() const;
  void SaveFieldToBinary(const std::string &filename);
  void SaveGradianceToBinary(const std::string &filename);

  /**
   * @brief Generalized sweep decomposition using iso-surface growth from
   * surface patches. Supports translational, rotational, and radial sweeps.
   * @param angularThreshold Angle threshold (radians) for direction constraint
   */
  void GeneralizedSweepDecomposition(float angularThreshold = 0.3f,
                                   bool cylinderPairsOnly = false);
  void AppendCylinderSweepDecomposition(float angularThreshold = 0.3f);
  /**
   * @brief 将 K≠0 的二次解析片替换为平面/柱面（K=0），并更新顶点法向。
   */
  void ReducePrimesToDevelopable(
      double curvatureThreshold = 1e-4,
      const std::string &exportPath = "");
  void RunCuttingBoxPipeline(bool cutMesh = false);
  /**
   * @brief 将模型分解为恰好两个扫掠体：
   *  - 柱面径向扫掠体（空心柱面，从内向外）
   *  - 垂直平移扫掠体（带柱孔的底座，沿柱轴方向）
   * 使用距离场体素归属来标记两个扫掠分块。
   */
  void DecomposeIntoTwoSweepBodies(float angularThreshold = 0.3f);
  void ApplySweepVisualization();
  bool HasNonPlanarPrimes() const;
  // Hex meshing is optional and must not change CuttingHex / mesh blocking.
  void GenerateSweepHexMeshes(int divisionsU = 8, int divisionsV = 8,
                              int divisionsW = 4, float targetCellSize = 0.0f);
  bool WriteSweepHexMeshesVTK(const std::string &path) const;
  const std::vector<SweepHexMesh> &GetSweepHexMeshes() const {
    return sweepHexMeshes;
  }
  bool PrimeLabelValid(int label) const;
  const PrimeData *GetPrimeByLabel(int label) const;
  std::vector<SweepBlockRegion> GetSweepBlocks() const {
    return this->sweepBlocks;
  }
  const std::vector<CylinderPairViz> &GetCylinderPairViz() const {
    return cylinderPairViz;
  }

protected:
  MeshLib::CTMesh *mesh;
  float PatchSize;
  std::vector<std::vector<float>> PointList;
  std::vector<int> PointIDList;
  std::vector<MeshLib::CToolVertex *> VertexPtrList;
  std::vector<std::vector<std::vector<float>>> Field;
  std::vector<std::vector<std::vector<int>>> FieldLabel;
  std::vector<std::vector<std::vector<float>>> GradianceDiff;
  std::vector<std::vector<std::vector<int>>> GradianceCount;
  std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists;

  std::vector<std::vector<std::vector<Eigen::Vector3f>>> GradianceField;
  std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord;
  std::vector<Eigen::Vector3f> SweepDir;
  std::vector<PrimeData> primes;
  std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjScalar;
  std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjEnergy;
  std::vector<std::string> sweepEnergyNames;

  int maxPointsPerNode = 32;
  int maxDepth = 8;
  void BuildOctree();
  void BuildOctreeRecursive(std::shared_ptr<OctreeNode> node,
                            const std::vector<int> &pointIndices, int depth);
  void SweepProjection_Regist(bool cutMesh = false);
  /**
   * @brief 计算各扫掠方向的投影标量与组合能量（Filter + Spliter + Alpha）。
   * @return 能量场数量（与 SweepDir 一一对应）
   */
  int ComputeSweepDirectionEnergies();
  /**
   * @brief 对每个已算好的扫掠能量，按对应方向各建一个 CuttingBox。
   * 规则：有几个能量就建几个框。
   */
  void BuildCuttingBoxesFromEnergies(bool cutMesh = false);
  void InitForbiddenBoundaryPoints();
  void ReindexPrimesById();
  void SubdivideNode(std::shared_ptr<OctreeNode> node);

  Eigen::Vector4f ComputeVertexDistance(const Eigen::Vector3f &point,
                                        MeshLib::CToolVertex *nearestVertex,
                                        int x, int y, int z);
  int FindNearestPointInOctree(const Eigen::Vector3f &point,
                               std::shared_ptr<OctreeNode> node,
                               float &bestDist);
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
  void DFS(MeshLib::CToolVertex *vert, int label);

  bool insideCuttingBox(Eigen::Vector3f point,
                        const std::map<int, Eigen::Vector3f> &verticesMap) const;
  std::shared_ptr<OctreeNode> octreeRoot;
  std::vector<SweepBlockRegion> sweepBlocks;
  std::vector<SweepHexMesh> sweepHexMeshes;
  std::vector<CylinderPairViz> cylinderPairViz;
  std::vector<bool> sweepBlockNonPlanar;
  std::vector<Eigen::Vector3f> sweepBlockColors;
  std::map<VoxelIndex, int> voxelToBlock;

  VoxelIndex WorldToVoxel(const Eigen::Vector3f &p) const;
  std::map<int, Eigen::Vector3f>
  BuildOrientedHex(const std::vector<Eigen::Vector3f> &pts,
                   const Eigen::Vector3f &axis,
                   const Eigen::Vector3f &origin) const;

  static Eigen::Vector3f RandomSweepColor(int seed);
  void EnsureSweepBlockColors();
  void AppendSweepBlocks(const std::vector<SweepBlockRegion> &blocks,
                         const std::vector<std::map<int, Eigen::Vector3f>> &hexes,
                         bool nonPlanar);
  int FindSweepBlockForPoint(const Eigen::Vector3f &position) const;
  int FindSweepBlockForPrimeLabel(int primeLabel) const;
  int SweepBlockToHexIndex(int sweepBlockIdx) const;
  int HexIndexToSweepBlockIndex(int hexIdx) const;
  int FindBasePlanarHexIndex() const;
  bool IsPointInCylinderSweepBlock(int hexIdx,
                                   const Eigen::Vector3f &position) const;
  static float CuttingHexVolume(const std::map<int, Eigen::Vector3f> &verticesMap);
  int FindSweepBlockForFace(const Eigen::Vector3f &position,
                            const std::unordered_map<int, int> &labelVotes) const;
  Eigen::Vector3f SweepDirectionAt(int blockIdx,
                                 const Eigen::Vector3f &position) const;
  static Eigen::Vector3f EncodeSweepDirColor(const Eigen::Vector3f &dir);
};

#endif
