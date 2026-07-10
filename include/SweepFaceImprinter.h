#ifndef __SWEEP_FACE_IMPRINTER_H__
#define __SWEEP_FACE_IMPRINTER_H__

#include "CTMesh.h"
#include "SweepBlock.h"
#include <Eigen/Eigen>
#include <array>
#include <string>
#include <vector>

struct SweepHexMesh;

/// 扫掠块在截面平面上的定向包围盒（底面参数域）
struct SweepCapFrame {
  Eigen::Vector3f origin = Eigen::Vector3f::Zero();
  Eigen::Vector3f axis = Eigen::Vector3f::UnitY();
  Eigen::Vector3f crossY = Eigen::Vector3f::UnitX();
  Eigen::Vector3f crossZ = Eigen::Vector3f::UnitZ();
  float axMin = 0.0f;
  float axMax = 0.0f;
  float crossMinY = 0.0f;
  float crossMaxY = 0.0f;
  float crossMinZ = 0.0f;
  float crossMaxZ = 0.0f;
};

/// 压印后的底面四边形网格 + 沿扫掠方向的层数
struct ImprintedCapMesh {
  std::vector<Eigen::Vector3f> capNodes; // 底面节点（3D）
  std::vector<std::array<int, 4>> quads;
  std::vector<float> uSplits;           // 截面 crossY 方向分割
  std::vector<float> vSplits;           // 截面 crossZ 方向分割
  std::vector<std::array<Eigen::Vector3f, 2>> imprintEdges; // 压印特征线（3D，在底面上）
  int sweepLayers = 1;
};

class SweepFaceImprinter {
public:
  static bool BuildFrameFromHex(const std::map<int, Eigen::Vector3f> &hex,
                                const Eigen::Vector3f &axisHint,
                                const Eigen::Vector3f &originHint,
                                SweepCapFrame &frame);

  /// 提取底/顶面特征边，将顶面特征压印到底面参数域。
  /// @param buildLocalQuads false 时仅压印分割，四边形由 Gmsh 生成
  static bool ImprintTranslationalCap(
      MeshLib::CTMesh *mesh, const SweepBlockRegion &block,
      const SweepCapFrame &frame, int minDivU, int minDivV, int sweepLayers,
      float holeRadius, ImprintedCapMesh &out, bool buildLocalQuads = false);

  static void BuildLocalQuadsFromSplits(const SweepCapFrame &frame,
                                        float holeRadius,
                                        ImprintedCapMesh &cap);

  static void SweepCapToHex(const ImprintedCapMesh &cap,
                            const SweepCapFrame &frame, int blockIndex,
                            const std::string &name, SweepHexMesh &out);

  static void SweepCylindricalCapToHex(const ImprintedCapMesh &innerCap,
                                       const SweepBlockRegion &block,
                                       int radialLayers, int blockIndex,
                                       const std::string &name,
                                       SweepHexMesh &out);
};

#endif
