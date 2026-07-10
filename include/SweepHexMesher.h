#ifndef __SWEEP_HEX_MESHER_H__
#define __SWEEP_HEX_MESHER_H__

#include "SweepBlock.h"
#include <Eigen/Eigen>
#include <array>
#include <map>
#include <string>
#include <vector>

/// 单个扫掠分块生成的体网格（HEX8 + WEDGE6）
struct SweepHexMesh {
  int blockIndex = -1;
  std::string name;
  SweepKind kind = SweepKind::Translational;
  std::vector<Eigen::Vector3f> nodes;
  std::vector<std::array<int, 8>> hexes;
  std::vector<std::array<int, 6>> wedges; // 三角形截面挤出的三棱柱
  /// 压印到底面的特征线（可视化）
  std::vector<std::array<Eigen::Vector3f, 2>> imprintEdges;
};

struct SweepHexMesherConfig {
  int divisionsU = 8; // 底面四边形网格 U 方向
  int divisionsV = 8; // 底面四边形网格 V 方向
  int divisionsW = 4; // 沿扫掠方向的层数
  float targetCellSize = 0.0f; // >0 时按块尺寸自动估计层数
  bool useGmsh = true;         // 压印后用 Gmsh 做 quasi-structured 四边形剖分
  std::string gmshWorkDir = "gmsh_cap_work";
};

/**
 * @brief 从扫掠分块在源/目标面生成四边形网格，并沿扫掠方向拉出六面体。
 *
 * - Translational：底面/顶面为截面四边形，沿 sweepAxis 平移扫掠
 * - CylindricalBase：内柱面/外柱面为源/目标，参数 (theta, axial) 四边形，沿径向扫掠
 */
class SweepHexMesher {
public:
  static std::vector<SweepHexMesh>
  Generate(const std::vector<SweepBlockRegion> &blocks,
           const std::vector<std::map<int, Eigen::Vector3f>> &cuttingHexes,
           const SweepHexMesherConfig &cfg = {},
           MeshLib::CTMesh *surfaceMesh = nullptr);

  static bool WriteVTK(const std::vector<SweepHexMesh> &meshes,
                       const std::string &path);

  /// 提取六面体线框，供 Polyscope 预览
  static void BuildWireframe(
      const SweepHexMesh &mesh,
      std::vector<Eigen::Vector3f> &vertices,
      std::vector<std::array<size_t, 2>> &edges);
};

#endif
