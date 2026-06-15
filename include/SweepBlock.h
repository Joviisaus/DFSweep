#ifndef __SWEEP_BLOCK_H__
#define __SWEEP_BLOCK_H__

#include "CTMesh.h"
#include "PrimeData.h"
#include <Eigen/Eigen>
#include <limits>
#include <map>
#include <set>
#include <vector>

struct VoxelIndex {
  int x, y, z;
  bool operator==(const VoxelIndex &o) const {
    return x == o.x && y == o.y && z == o.z;
  }
  bool operator<(const VoxelIndex &o) const {
    if (x != o.x)
      return x < o.x;
    if (y != o.y)
      return y < o.y;
    return z < o.z;
  }
};

struct VoxelIndexHash {
  std::size_t operator()(const VoxelIndex &v) const {
    std::size_t h = std::hash<int>()(v.x);
    h ^= std::hash<int>()(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>()(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

enum class SweepKind {
  Translational,  // 平移扫掠：沿固定方向
  Rotational,     // 旋转扫掠：绕轴的切向
  Radial,         // 中心向外：从中心指向外的径向
  CylindricalBase // 柱面：从柱轴向内向外径向扫掠，配对柱面为外底面
};

/**
 * @brief 从解析表面片向内生长得到的扫掠体分块
 * 拓扑：外底面（表面片）+ 内底面（最深等值面）+ 侧面（满足约束的扫掠边界）
 */
struct SweepBlockRegion {
  int patchId;
  int pairedPatchId = -1; // 柱面配对时的另一底面 patch 索引，-1 表示无
  int primeId;
  int pairedPrimeId = -1;
  SweepKind kind;
  float minIsoValue; // 最深等值面（最负）
  float maxIsoValue; // 表面附近（≈0）
  float axialLower = 0.0f; // 柱轴方向范围下界（相对 sweepOrigin）
  float axialUpper = 0.0f; // 柱轴方向范围上界
  float radialInner = 0.0f; // 径向扫掠内界（近柱轴）
  float radialOuter = 0.0f; // 径向扫掠外界（柱面）
  int bottomPatchId = -1;   // 上下底面中的底面 patch
  std::vector<int> memberPrimeIds; // 该扫掠区域包含的全部柱面 prime
  std::set<VoxelIndex> coveredVoxels;
  Eigen::Vector3f sweepAxis;   // 平移方向 / 旋转轴 / 径向参考轴
  Eigen::Vector3f sweepOrigin; // 旋转中心 / 径向中心
  Eigen::Vector3f crossDirY;   // 垂直于扫掠轴的截面基向量
  Eigen::Vector3f crossDirZ;
  bool isValid;
};

/** 柱面配对可视化：两片柱面 + 径向扫掠内外界 */
struct CylinderPairViz {
  int patchIdA = -1;
  int patchIdB = -1;
  int primeIdA = -1;
  int primeIdB = -1;
  Eigen::Vector3f sweepAxis = Eigen::Vector3f::UnitY();
  Eigen::Vector3f sweepOrigin = Eigen::Vector3f::Zero();
  Eigen::Vector3f crossDirY = Eigen::Vector3f::UnitX();
  Eigen::Vector3f crossDirZ = Eigen::Vector3f::UnitZ();
  float axialMid = 0.0f;
  float radialInner = 0.0f;
  float radialOuter = 0.0f;
  float innerCrossMinY = 0.0f;
  float innerCrossMaxY = 0.0f;
  float innerCrossMinZ = 0.0f;
  float innerCrossMaxZ = 0.0f;
  float outerCrossMinY = 0.0f;
  float outerCrossMaxY = 0.0f;
  float outerCrossMinZ = 0.0f;
  float outerCrossMaxZ = 0.0f;
};

/**
 * @brief 广义扫掠体分块（与 CuttingBox 正交）
 *
 * 1. 按解析式表面片（PrimeData）分块种子
 * 2. 从表面沿距离场等值面向内生长，区域尽可能大
 * 3. 中间体素：场方向与局部等值面法向垂直或相切（用户阈值）
 * 4. 贪心集合覆盖，用最少的分块覆盖全部内部体素
 * 5. 校验两底面 + 一侧面 的扫掠体拓扑
 */
class SweepDecomposer {
public:
  SweepDecomposer(
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
      const std::vector<std::vector<std::vector<float>>> &Field,
      const std::vector<std::vector<std::vector<int>>> &FieldLabel,
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField,
      const std::vector<PrimeData> &primes, float angularThreshold,
      bool cylinderPairsOnly = false,
      MeshLib::CTMesh *mesh = nullptr);

  std::vector<SweepBlockRegion> GetBlocks() const { return selectedBlocks; }
  std::vector<std::map<int, Eigen::Vector3f>> GetBlockHexVertices() const;
  std::vector<std::pair<int, int>> GetCylinderPairs() const {
    return cylinderPairs;
  }
  const std::vector<CylinderPairViz> &GetCylinderPairViz() const {
    return cylinderPairViz;
  }

protected:
  const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord;
  const std::vector<std::vector<std::vector<float>>> &Field;
  const std::vector<std::vector<std::vector<int>>> &FieldLabel;
  const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField;
  const std::vector<PrimeData> &primes;
  MeshLib::CTMesh *mesh;
  float threshold;
  float stepSize;
  int D1, D2, D3;

  bool cylinderPairsOnly;
  std::vector<std::pair<int, int>> cylinderPairs;
  std::vector<CylinderPairViz> cylinderPairViz;
  std::vector<SweepBlockRegion> candidateBlocks;
  std::vector<SweepBlockRegion> selectedBlocks;
  std::vector<std::map<int, Eigen::Vector3f>> blockHexVertices;
  std::set<VoxelIndex> allInteriorVoxels;

  void CollectInteriorVoxels();
  void GrowBlockFromPatch(int patchId);
  void GrowBlockBetweenCylinderPair(int patchA, int patchB);
  void GrowCylinderSweepRegion(const std::vector<int> &groupPatchIds);
  std::vector<std::vector<int>> BuildCylinderSweepGroups() const;
  float CylinderSweepEnergy(int x, int y, int z,
                            const Eigen::Vector3f &sweepDir) const;
  Eigen::Vector3f LocalCylinderSweepDirection(
      int x, int y, int z, const Eigen::Vector3f &axis,
      const Eigen::Vector3f &origin, int bottomPrimeId, float axBottom,
      const Eigen::Vector3f &bottomSweepDir) const;
  void BuildCylinderPairBlockFromMesh(int patchA, int patchB);
  void BuildCylinderPatchBlockFromMesh(int patchId);
  void BuildCylinderSinglePatchBlockFromMesh(
      int patchId, int pairedPatchId, const Eigen::Vector3f &axis,
      const Eigen::Vector3f &origin, const std::vector<Eigen::Vector3f> &patchPts,
      const CylinderPairViz &pairViz);
  std::vector<std::pair<int, int>> FindMatchingCylinderPairs() const;
  Eigen::Vector3f ComputePatchCentroidFromMesh(int patchId) const;
  Eigen::Vector3f ComputePatchAvgGradientFromMesh(int patchId,
                                                  int *sampleCount) const;
  Eigen::Vector3f EvalPrimeGradient(int patchId,
                                    const Eigen::Vector3f &p) const;
  Eigen::Vector3f
  ComputePatchAvgAnalyticGradientFromMesh(int patchId,
                                          int *sampleCount) const;
  std::map<int, Eigen::Vector3f>
  BuildHexFromPoints(const std::vector<Eigen::Vector3f> &points,
                     const Eigen::Vector3f &axis,
                     const Eigen::Vector3f &origin) const;
  std::map<int, Eigen::Vector3f>
  BuildCylinderRadialHex(const std::vector<Eigen::Vector3f> &points,
                         const Eigen::Vector3f &axisDir,
                         const Eigen::Vector3f &axisOrigin, float rInner,
                         float rOuter, float axMin, float axMax,
                         const Eigen::Vector3f &crossY,
                         const Eigen::Vector3f &crossZ) const;
  CylinderPairViz
  BuildCylinderPairViz(int patchA, int patchB,
                       const Eigen::Vector3f &axis,
                       const Eigen::Vector3f &origin,
                       const std::vector<Eigen::Vector3f> &ptsA,
                       const std::vector<Eigen::Vector3f> &ptsB) const;
  bool IsCylinderPatch(int patchId) const;
  bool GradientsMatch(int patchA, int patchB) const;
  Eigen::Vector3f ExtractCylinderAxis(int patchId) const;
  Eigen::Vector3f ComputePatchCentroid(int patchId) const;
  Eigen::Vector3f ComputePatchAvgGradient(int patchId, int *sampleCount) const;
  SweepKind ClassifyPatch(int patchId) const;
  Eigen::Vector3f ComputePatchFrame(int patchId, SweepKind kind,
                                    Eigen::Vector3f &outOrigin) const;
  Eigen::Vector3f LocalSweepDirection(SweepKind kind, int x, int y, int z,
                                      const Eigen::Vector3f &axis,
                                      const Eigen::Vector3f &origin) const;
  Eigen::Vector3f ComputeIsoSurfaceNormal(int x, int y, int z) const;
  bool CheckDirectionConstraint(int x, int y, int z,
                                const Eigen::Vector3f &isoNormal) const;
  float SweepParameter(SweepKind kind, const Eigen::Vector3f &p,
                       const Eigen::Vector3f &axis,
                       const Eigen::Vector3f &origin) const;
  bool ValidateSweepTopology(const SweepBlockRegion &block) const;
  void GreedySetCover();
};

#endif
