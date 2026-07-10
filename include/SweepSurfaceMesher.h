#ifndef __SWEEP_SURFACE_MESHER_H__
#define __SWEEP_SURFACE_MESHER_H__

#include "CTMesh.h"
#include "SweepBlock.h"
#include "SweepFaceImprinter.h"
#include "SweepHexMesher.h"
#include <Eigen/Eigen>
#include <array>
#include <vector>

/// 原表面源面三角片（柱面侧壁 / 带孔底面），紧贴输入网格
struct SurfaceCapPatch {
  int blockIndex = -1;
  std::vector<int> faceIds;
  std::vector<Eigen::Vector3f> nodes;
  std::vector<std::array<int, 3>> tris;
  std::vector<std::array<Eigen::Vector3f, 2>> boundaryEdges;
  std::vector<std::array<Eigen::Vector3f, 2>> holeEdges; // 镂空内边界
  float holeRadius = 0.0f;
  float outerRadius = 0.0f;
};

/// 四边形占优截面：quad + 剩余 triangle（Q-Morph 风格）
struct QuadDominantCap {
  std::vector<Eigen::Vector3f> nodes;
  std::vector<std::array<int, 4>> quads;
  std::vector<std::array<int, 3>> tris;
  std::vector<std::array<Eigen::Vector3f, 2>> imprintEdges;
  int sweepLayers = 1;
};

/**
 * @brief 从原表面提取扫掠源面，做四边形占优剖分，再沿扫掠方向挤出。
 *
 * sweepFaceType:
 *  - 5 = 平移块底面（带孔圆环，原表面）
 *  - 6 = 柱面块源侧壁（内柱面，原表面）
 *  - 7 = 柱面块目标侧壁（外柱面，原表面）
 *
 * 扫掠深度/半径优先用原 mesh 面片几何，不依赖 cutting box。
 */
class SweepSurfaceMesher {
public:
  /// 标注平移块底面(5) + 柱面内外侧壁(6/7)
  static int MarkSourceFacesOnSurface(MeshLib::CTMesh *mesh,
                                      const std::vector<SweepBlockRegion> &blocks,
                                      const std::vector<bool> &blockNonPlanar);

  /// 兼容旧接口
  static int MarkBottomFacesOnSurface(MeshLib::CTMesh *mesh,
                                      const std::vector<SweepBlockRegion> &blocks,
                                      const std::vector<bool> &blockNonPlanar) {
    return MarkSourceFacesOnSurface(mesh, blocks, blockNonPlanar);
  }

  /// 从原 mesh 面片统计扫掠框（不依赖 cutting hex）
  static bool BuildFrameFromMeshFaces(MeshLib::CTMesh *mesh, int blockIndex,
                                      const SweepBlockRegion &block,
                                      SweepCapFrame &frame);

  /// 提取平移块带孔底面（type 5），镂空边界写入 holeEdges
  static bool ExtractBottomCapPatch(MeshLib::CTMesh *mesh, int blockIndex,
                                    const SweepCapFrame &frame,
                                    float holeRadius, SurfaceCapPatch &out);

  /// 提取柱面内壁源面（type 6），作为径向扫掠源
  static bool ExtractCylindricalWallPatch(MeshLib::CTMesh *mesh, int blockIndex,
                                          const SweepBlockRegion &block,
                                          bool innerWall, SurfaceCapPatch &out);

  /// 简易 Q-Morph；若有 holeRadius>0，剔除孔内单元
  static bool BuildQuadDominantFromTris(const SurfaceCapPatch &patch,
                                        float minQuadQuality,
                                        QuadDominantCap &out);

  /// 平移扫掠：quad→hex，tri→wedge；顶层贴合柱面下沿接口（几何拼合）
  static bool SweepQuadDominantToVolume(
      const QuadDominantCap &cap, const SweepCapFrame &frame,
      MeshLib::CTMesh *mesh, int blockIndex, const std::string &name,
      int layers, SweepHexMesh &out,
      const SweepBlockRegion *interfaceTube = nullptr);

  /// 径向扫掠：内壁源面 → 外壁；节点贴合原表面内外柱面
  static bool SweepCylindricalWallToVolume(
      const QuadDominantCap &cap, const SweepBlockRegion &block,
      MeshLib::CTMesh *mesh, int blockIndex, const std::string &name,
      int radialLayers, SweepHexMesh &out);
};

#endif
