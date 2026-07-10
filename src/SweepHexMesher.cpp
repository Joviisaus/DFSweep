#include "SweepHexMesher.h"
#include "GmshCapMesher.h"
#include "Mesh/iterators.h"
#include "SweepFaceImprinter.h"
#include "SweepSurfaceMesher.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {
constexpr float kPi = 3.14159265358979323846f;

void BuildCrossFrame(const Eigen::Vector3f &axisDir, Eigen::Vector3f &crossY,
                     Eigen::Vector3f &crossZ) {
  Eigen::Vector3f arbitrary =
      (std::abs(axisDir.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axisDir.dot(arbitrary) * axisDir).normalized();
  crossZ = axisDir.cross(crossY);
  if (crossZ.norm() < 1e-8f) {
    crossZ = Eigen::Vector3f::UnitZ();
  } else {
    crossZ.normalize();
  }
}

struct OrientedBox {
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

bool ParseOrientedBoxFromHex(const std::map<int, Eigen::Vector3f> &hex,
                             const Eigen::Vector3f &axisHint,
                             const Eigen::Vector3f &originHint,
                             OrientedBox &box) {
  if (hex.size() != 8) {
    return false;
  }
  box.origin = originHint;
  box.axis = axisHint.normalized();
  BuildCrossFrame(box.axis, box.crossY, box.crossZ);

  box.axMin = std::numeric_limits<float>::max();
  box.axMax = std::numeric_limits<float>::lowest();
  box.crossMinY = std::numeric_limits<float>::max();
  box.crossMaxY = std::numeric_limits<float>::lowest();
  box.crossMinZ = std::numeric_limits<float>::max();
  box.crossMaxZ = std::numeric_limits<float>::lowest();

  for (const auto &[_, p] : hex) {
    Eigen::Vector3f rel = p - box.origin;
    float a = rel.dot(box.axis);
    float cy = rel.dot(box.crossY);
    float cz = rel.dot(box.crossZ);
    box.axMin = std::min(box.axMin, a);
    box.axMax = std::max(box.axMax, a);
    box.crossMinY = std::min(box.crossMinY, cy);
    box.crossMaxY = std::max(box.crossMaxY, cy);
    box.crossMinZ = std::min(box.crossMinZ, cz);
    box.crossMaxZ = std::max(box.crossMaxZ, cz);
  }
  return box.axMax > box.axMin + 1e-6f;
}

int ClampDiv(int value, int minV = 2, int maxV = 48) {
  return std::max(minV, std::min(maxV, value));
}

int AutoDiv(float extent, float targetSize, int fallback) {
  if (targetSize <= 1e-8f || extent <= 1e-8f) {
    return fallback;
  }
  return ClampDiv(static_cast<int>(std::lround(extent / targetSize)));
}

template <typename NodeFn>
void AddStructuredHexes(int nu, int nv, int nw, NodeFn &&nodeAt,
                        SweepHexMesh &mesh) {
  const int base = static_cast<int>(mesh.nodes.size());
  mesh.nodes.reserve(base + (nu + 1) * (nv + 1) * (nw + 1));
  for (int iw = 0; iw <= nw; ++iw) {
    for (int iv = 0; iv <= nv; ++iv) {
      for (int iu = 0; iu <= nu; ++iu) {
        mesh.nodes.push_back(nodeAt(iu, iv, iw));
      }
    }
  }

  auto nid = [&](int iu, int iv, int iw) {
    return base + iw * (nu + 1) * (nv + 1) + iv * (nu + 1) + iu;
  };

  for (int iw = 0; iw < nw; ++iw) {
    for (int iv = 0; iv < nv; ++iv) {
      for (int iu = 0; iu < nu; ++iu) {
        std::array<int, 8> hex = {
            nid(iu, iv, iw),         nid(iu + 1, iv, iw),
            nid(iu + 1, iv + 1, iw), nid(iu, iv + 1, iw),
            nid(iu, iv, iw + 1),     nid(iu + 1, iv, iw + 1),
            nid(iu + 1, iv + 1, iw + 1), nid(iu, iv + 1, iw + 1)};
        mesh.hexes.push_back(hex);
      }
    }
  }
}

SweepHexMesh GenerateTranslationalMesh(const SweepBlockRegion &block,
                                       const OrientedBox &box, int nu, int nv,
                                       int nw) {
  SweepHexMesh mesh;
  mesh.kind = SweepKind::Translational;
  const float du = 1.0f / static_cast<float>(nu);
  const float dv = 1.0f / static_cast<float>(nv);
  const float dw = 1.0f / static_cast<float>(nw);

  AddStructuredHexes(
      nu, nv, nw,
      [&](int iu, int iv, int iw) {
        float u = static_cast<float>(iu) * du;
        float v = static_cast<float>(iv) * dv;
        float w = static_cast<float>(iw) * dw;
        float a = box.axMin + w * (box.axMax - box.axMin);
        float cy = box.crossMinY + u * (box.crossMaxY - box.crossMinY);
        float cz = box.crossMinZ + v * (box.crossMaxZ - box.crossMinZ);
        return box.origin + a * box.axis + cy * box.crossY + cz * box.crossZ;
      },
      mesh);
  return mesh;
}

SweepHexMesh GenerateCylindricalRadialMesh(const SweepBlockRegion &block,
                                           int nTheta, int nAxial, int nRadial) {
  SweepHexMesh mesh;
  mesh.kind = SweepKind::CylindricalBase;

  Eigen::Vector3f axis = block.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f origin = block.sweepOrigin;
  Eigen::Vector3f crossY = block.crossDirY;
  Eigen::Vector3f crossZ = block.crossDirZ;
  if (crossY.norm() < 1e-6f || crossZ.norm() < 1e-6f ||
      std::abs(crossY.normalized().dot(axis)) > 0.99f) {
    BuildCrossFrame(axis, crossY, crossZ);
  } else {
    crossY.normalize();
    crossZ = axis.cross(crossY);
    if (crossZ.norm() < 1e-8f) {
      BuildCrossFrame(axis, crossY, crossZ);
    } else {
      crossZ.normalize();
    }
  }

  const float r0 = std::max(0.0f, block.radialInner);
  const float r1 = std::max(r0 + 1e-4f, block.radialOuter);
  const float a0 = block.axialLower;
  const float a1 = block.axialUpper;
  const float dTheta = 2.0f * kPi / static_cast<float>(nTheta);
  const float dAxial = 1.0f / static_cast<float>(nAxial);
  const float dRadial = 1.0f / static_cast<float>(nRadial);

  const int base = 0;
  mesh.nodes.reserve((nTheta + 1) * (nAxial + 1) * (nRadial + 1));
  for (int ir = 0; ir <= nRadial; ++ir) {
    for (int ia = 0; ia <= nAxial; ++ia) {
      for (int it = 0; it <= nTheta; ++it) {
        float theta = static_cast<float>(it) * dTheta;
        float ax = a0 + static_cast<float>(ia) * dAxial * (a1 - a0);
        float r = r0 + static_cast<float>(ir) * dRadial * (r1 - r0);
        Eigen::Vector3f radialDir =
            std::cos(theta) * crossY + std::sin(theta) * crossZ;
        mesh.nodes.push_back(origin + ax * axis + r * radialDir);
      }
    }
  }

  auto nid = [&](int it, int ia, int ir) {
    return base + ir * (nTheta + 1) * (nAxial + 1) + ia * (nTheta + 1) + it;
  };
  for (int ir = 0; ir < nRadial; ++ir) {
    for (int ia = 0; ia < nAxial; ++ia) {
      for (int it = 0; it < nTheta; ++it) {
        std::array<int, 8> hex = {
            nid(it, ia, ir),         nid(it + 1, ia, ir),
            nid(it + 1, ia + 1, ir), nid(it, ia + 1, ir),
            nid(it, ia, ir + 1),     nid(it + 1, ia, ir + 1),
            nid(it + 1, ia + 1, ir + 1), nid(it, ia + 1, ir + 1)};
        mesh.hexes.push_back(hex);
      }
    }
  }
  return mesh;
}
} // namespace

float FindHoleRadiusForBlock(
    const std::vector<SweepBlockRegion> &blocks, int blockIndex) {
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    if (i == blockIndex) {
      continue;
    }
    const SweepBlockRegion &other = blocks[static_cast<size_t>(i)];
    if (other.kind == SweepKind::CylindricalBase && other.radialInner > 0.0f) {
      return other.radialInner;
    }
  }
  return 0.0f;
}

const SweepBlockRegion *
FindTubeInterfaceBlock(const std::vector<SweepBlockRegion> &blocks,
                       int blockIndex, int *outTubeIndex = nullptr) {
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    if (i == blockIndex) {
      continue;
    }
    const SweepBlockRegion &other = blocks[static_cast<size_t>(i)];
    if (other.kind == SweepKind::CylindricalBase &&
        other.radialOuter > other.radialInner + 1e-6f) {
      if (outTubeIndex) {
        *outTubeIndex = i;
      }
      return &other;
    }
  }
  if (outTubeIndex) {
    *outTubeIndex = -1;
  }
  return nullptr;
}

/// 从已生成的柱面体网格取下沿轴向（只用单元实际引用的节点）
float TubeMeshLowerAxial(const SweepHexMesh &tubeMesh,
                         const SweepBlockRegion &tube) {
  if (tubeMesh.nodes.empty()) {
    return tube.axialLower;
  }
  Eigen::Vector3f axis = tube.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  const Eigen::Vector3f &origin = tube.sweepOrigin;

  std::vector<char> used(tubeMesh.nodes.size(), 0);
  for (const auto &h : tubeMesh.hexes) {
    for (int nid : h) {
      if (nid >= 0 && nid < static_cast<int>(used.size())) {
        used[static_cast<size_t>(nid)] = 1;
      }
    }
  }
  for (const auto &w : tubeMesh.wedges) {
    for (int nid : w) {
      if (nid >= 0 && nid < static_cast<int>(used.size())) {
        used[static_cast<size_t>(nid)] = 1;
      }
    }
  }

  float axMin = std::numeric_limits<float>::max();
  int cnt = 0;
  for (size_t i = 0; i < tubeMesh.nodes.size(); ++i) {
    if (!used[i]) {
      continue;
    }
    axMin = std::min(axMin, (tubeMesh.nodes[i] - origin).dot(axis));
    ++cnt;
  }
  if (cnt == 0) {
    for (const auto &p : tubeMesh.nodes) {
      axMin = std::min(axMin, (p - origin).dot(axis));
    }
  }
  return axMin;
}

const SweepHexMesh *FindGeneratedTubeMesh(const std::vector<SweepHexMesh> &meshes,
                                          int tubeBlockIndex) {
  for (const auto &m : meshes) {
    if (m.blockIndex == tubeBlockIndex &&
        m.kind == SweepKind::CylindricalBase && !m.nodes.empty()) {
      return &m;
    }
  }
  return nullptr;
}

std::vector<SweepHexMesh>
SweepHexMesher::Generate(const std::vector<SweepBlockRegion> &blocks,
                         const std::vector<std::map<int, Eigen::Vector3f>> &cuttingHexes,
                         const SweepHexMesherConfig &cfg,
                         MeshLib::CTMesh *surfaceMesh) {
  std::vector<SweepHexMesh> meshes;
  meshes.reserve(blocks.size());

  for (int bi = 0; bi < static_cast<int>(blocks.size()); ++bi) {
    const SweepBlockRegion &block = blocks[static_cast<size_t>(bi)];
    if (!block.isValid) {
      continue;
    }

    SweepHexMesh hexMesh;
    hexMesh.blockIndex = bi;
    hexMesh.name = (block.kind == SweepKind::CylindricalBase)
                       ? "CylinderSweepHex"
                       : "TranslationalSweepHex";
    hexMesh.name += std::to_string(bi);

    if (block.kind == SweepKind::CylindricalBase) {
      int nTheta = cfg.divisionsU;
      int nAxial = cfg.divisionsV;
      int nRadial = cfg.divisionsW;
      if (cfg.targetCellSize > 0.0f) {
        nTheta = AutoDiv(2.0f * kPi * block.radialOuter, cfg.targetCellSize,
                         nTheta);
        nAxial = AutoDiv(block.axialUpper - block.axialLower,
                         cfg.targetCellSize, nAxial);
        nRadial = AutoDiv(block.radialOuter - block.radialInner,
                          cfg.targetCellSize, nRadial);
      }
      nTheta = ClampDiv(nTheta);
      nAxial = ClampDiv(nAxial);
      nRadial = ClampDiv(nRadial, 2, 24);

      bool usedSurface = false;
      // 优先：原 mesh 内柱面侧壁 → Q-Morph → 径向扫掠到外壁（几何拼合）
      if (surfaceMesh) {
        SurfaceCapPatch wall;
        if (SweepSurfaceMesher::ExtractCylindricalWallPatch(
                surfaceMesh, bi, block, /*innerWall=*/true, wall)) {
          QuadDominantCap qd;
          if (SweepSurfaceMesher::BuildQuadDominantFromTris(wall, 0.05f, qd)) {
            if (SweepSurfaceMesher::SweepCylindricalWallToVolume(
                    qd, block, surfaceMesh, bi,
                    "CylinderSweepHex" + std::to_string(bi), nRadial, hexMesh)) {
              usedSurface = true;
              std::cout << "[SweepHexMesher] block " << bi
                        << " surface wall radial sweep: quads=" << qd.quads.size()
                        << " tris=" << qd.tris.size() << " -> hex="
                        << hexMesh.hexes.size()
                        << " wedge=" << hexMesh.wedges.size() << "\n";
            }
          }
        }
      }

      bool usedGmsh = false;
      if (!usedSurface && cfg.useGmsh && GmshCapMesher::IsAvailable()) {
        SweepBlockRegion cylBlock = block;
        Eigen::Vector3f cylAxis = cylBlock.sweepAxis;
        if (cylAxis.norm() < 1e-8f) {
          cylAxis = Eigen::Vector3f::UnitY();
        } else {
          cylAxis.normalize();
        }
        cylBlock.sweepAxis = cylAxis;
        BuildCrossFrame(cylAxis, cylBlock.crossDirY, cylBlock.crossDirZ);

        ImprintedCapMesh innerCap;
        const std::string workDir =
            cfg.gmshWorkDir + "/block_" + std::to_string(bi);
        if (GmshCapMesher::MeshCylindricalInnerCap(cylBlock, nTheta, nAxial,
                                                   innerCap, workDir)) {
          const size_t minQuads =
              static_cast<size_t>(std::max(1, nTheta * nAxial / 2));
          if (innerCap.quads.size() < minQuads) {
            std::cerr << "[SweepHexMesher] block " << bi
                      << " gmsh cylindrical cap too few quads ("
                      << innerCap.quads.size() << "), fallback\n";
          } else {
            SweepFaceImprinter::SweepCylindricalCapToHex(
                innerCap, cylBlock, nRadial, bi,
                "CylinderSweepHex" + std::to_string(bi), hexMesh);
            usedGmsh = !hexMesh.hexes.empty();
            if (usedGmsh) {
              std::cout << "[SweepHexMesher] block " << bi
                        << " gmsh inner cap + radial sweep: "
                        << innerCap.quads.size() << " quads x " << nRadial
                        << " radial -> " << hexMesh.hexes.size() << " elements\n";
            }
          }
        }
      }
      if (!usedSurface && !usedGmsh) {
        hexMesh = GenerateCylindricalRadialMesh(block, nTheta, nAxial, nRadial);
        std::cout << "[SweepHexMesher] block " << bi
                  << " cylindrical radial hex (fallback): " << nTheta << "x"
                  << nAxial << "x" << nRadial << " -> " << hexMesh.hexes.size()
                  << " elements\n";
      }
    } else {
      // 扫掠框优先从原 mesh 源/目标面统计，不依赖 cutting box
      SweepCapFrame frame;
      bool haveFrame = false;
      if (surfaceMesh) {
        haveFrame = SweepSurfaceMesher::BuildFrameFromMeshFaces(
            surfaceMesh, bi, block, frame);
      }
      if (!haveFrame && bi >= 0 && bi < static_cast<int>(cuttingHexes.size())) {
        haveFrame = SweepFaceImprinter::BuildFrameFromHex(
            cuttingHexes[static_cast<size_t>(bi)], block.sweepAxis,
            block.sweepOrigin, frame);
      }
      if (!haveFrame) {
        std::cerr << "[SweepHexMesher] block " << bi
                  << " missing mesh frame / cutting hex, skipped\n";
        continue;
      }

      int nu = cfg.divisionsU;
      int nv = cfg.divisionsV;
      int nw = cfg.divisionsW;
      bool usedSurface = false;
      float holeR = FindHoleRadiusForBlock(blocks, bi);

      // 优先：原表面带孔底面 → Q-Morph → 轴向扫掠到已生成柱面 hex 下沿
      int tubeIdx = -1;
      SweepBlockRegion tubeIfaceLocal;
      const SweepBlockRegion *tubeIface = nullptr;
      if (const SweepBlockRegion *found =
              FindTubeInterfaceBlock(blocks, bi, &tubeIdx)) {
        tubeIfaceLocal = *found;
        float ifaceAx = found->axialLower;
        if (const SweepHexMesh *tubeMesh =
                FindGeneratedTubeMesh(meshes, tubeIdx)) {
          ifaceAx = TubeMeshLowerAxial(*tubeMesh, tubeIfaceLocal);
        }
        tubeIfaceLocal.axialLower = ifaceAx;
        tubeIface = &tubeIfaceLocal;
        std::cout << "[SweepHexMesher] block " << bi
                  << " tube interface: ifaceAx=" << ifaceAx
                  << " (blockAx=" << found->axialLower << ") r=["
                  << found->radialInner << "," << found->radialOuter << "]\n";
      }
      if (surfaceMesh) {
        SurfaceCapPatch patch;
        if (SweepSurfaceMesher::ExtractBottomCapPatch(surfaceMesh, bi, frame,
                                                      holeR, patch)) {
          QuadDominantCap qd;
          if (SweepSurfaceMesher::BuildQuadDominantFromTris(patch, 0.05f, qd)) {
            qd.sweepLayers = nw;
            if (SweepSurfaceMesher::SweepQuadDominantToVolume(
                    qd, frame, surfaceMesh, bi,
                    "TranslationalSweepHex" + std::to_string(bi), nw, hexMesh,
                    tubeIface)) {
              // 顶层环带：用柱面 hex 节点的绝对坐标拼合（不依赖 axial 标架一致性）
              if (tubeIface) {
                if (const SweepHexMesh *tubeMesh =
                        FindGeneratedTubeMesh(meshes, tubeIdx)) {
                  Eigen::Vector3f axis = tubeIface->sweepAxis;
                  if (axis.norm() < 1e-8f) {
                    axis = Eigen::Vector3f::UnitY();
                  } else {
                    axis.normalize();
                  }
                  const Eigen::Vector3f &origin = tubeIface->sweepOrigin;
                  const float rIn = tubeIface->radialInner;
                  const float rOut = tubeIface->radialOuter;

                  // 柱面体网格实际下沿：只用单元引用节点
                  std::vector<char> used(tubeMesh->nodes.size(), 0);
                  for (const auto &h : tubeMesh->hexes) {
                    for (int nid : h) {
                      if (nid >= 0 &&
                          nid < static_cast<int>(used.size())) {
                        used[static_cast<size_t>(nid)] = 1;
                      }
                    }
                  }
                  for (const auto &w : tubeMesh->wedges) {
                    for (int nid : w) {
                      if (nid >= 0 &&
                          nid < static_cast<int>(used.size())) {
                        used[static_cast<size_t>(nid)] = 1;
                      }
                    }
                  }
                  float tubeAxMin = std::numeric_limits<float>::max();
                  float tubeYMin = std::numeric_limits<float>::max();
                  float tubeYMax = std::numeric_limits<float>::lowest();
                  for (size_t i = 0; i < tubeMesh->nodes.size(); ++i) {
                    if (!used[i]) {
                      continue;
                    }
                    const auto &q = tubeMesh->nodes[i];
                    tubeAxMin =
                        std::min(tubeAxMin, (q - origin).dot(axis));
                    tubeYMin = std::min(tubeYMin, q.y());
                    tubeYMax = std::max(tubeYMax, q.y());
                  }
                  const float axTol = 0.02f * std::max(
                      1e-3f, tubeIface->axialUpper - tubeIface->axialLower);

                  std::vector<Eigen::Vector3f> rim;
                  rim.reserve(tubeMesh->nodes.size() / 8 + 8);
                  for (size_t i = 0; i < tubeMesh->nodes.size(); ++i) {
                    if (!used[i]) {
                      continue;
                    }
                    const auto &q = tubeMesh->nodes[i];
                    float a = (q - origin).dot(axis);
                    if (a > tubeAxMin + axTol) {
                      continue;
                    }
                    float r = ((q - origin) - a * axis).norm();
                    if (r >= rIn * 0.9f && r <= rOut * 1.1f) {
                      rim.push_back(q);
                    }
                  }

                  const int nCap = static_cast<int>(qd.nodes.size());
                  const int topBase = nw * nCap;
                  float topAxBeforeMin = std::numeric_limits<float>::max();
                  float topAxBeforeMax = std::numeric_limits<float>::lowest();
                  for (int i = 0; i < nCap; ++i) {
                    float a =
                        (hexMesh.nodes[static_cast<size_t>(topBase + i)] -
                         origin)
                            .dot(axis);
                    topAxBeforeMin = std::min(topAxBeforeMin, a);
                    topAxBeforeMax = std::max(topAxBeforeMax, a);
                  }

                  int snapped = 0;
                  int lifted = 0;
                  int leftOuter = 0;
                  float topYMin = std::numeric_limits<float>::max();
                  float topYMax = std::numeric_limits<float>::lowest();
                  for (int i = 0; i < nCap; ++i) {
                    Eigen::Vector3f &p =
                        hexMesh.nodes[static_cast<size_t>(topBase + i)];
                    Eigen::Vector3f rel = p - origin;
                    float a = rel.dot(axis);
                    Eigen::Vector3f rad = rel - a * axis;
                    float r = rad.norm();

                    // 只拼接管壁环带内的顶层节点；外缘底板保持原径向，
                    // 避免把底板外缘“粘”到柱壁上形成尖刺
                    if (r < rIn - 0.05f || r > rOut + 0.08f) {
                      // 外缘/孔内：仅抬到接口平面，不吸附 rim
                      if (rad.norm() > 1e-8f) {
                        p = origin + tubeAxMin * axis + rad;
                      } else {
                        p = origin + tubeAxMin * axis;
                      }
                      ++leftOuter;
                      ++lifted;
                    } else {
                      p = origin + tubeAxMin * axis + rad;
                      ++lifted;
                      if (!rim.empty()) {
                        float best = std::numeric_limits<float>::max();
                        Eigen::Vector3f bestP = p;
                        for (const auto &q : rim) {
                          float d = (q - p).squaredNorm();
                          if (d < best) {
                            best = d;
                            bestP = q;
                          }
                        }
                        // 只允许很近的吸附，防止跨到错误方位
                        if (best < 0.08f * 0.08f) {
                          p = bestP;
                          ++snapped;
                        }
                      }
                    }
                    topYMin = std::min(topYMin, p.y());
                    topYMax = std::max(topYMax, p.y());
                  }
                  std::cout << "[SweepHexMesher] block " << bi
                            << " join tube: tubeAxMin=" << tubeAxMin
                            << " topAxBefore=[" << topAxBeforeMin << ","
                            << topAxBeforeMax << "] lifted=" << lifted
                            << " snapped=" << snapped << "/" << nCap
                            << " leftOuter=" << leftOuter
                            << " rim=" << rim.size()
                            << " topY=[" << topYMin << "," << topYMax << "]"
                            << " tubeY=[" << tubeYMin << "," << tubeYMax << "]"
                            << "\n";
                }
              }
              usedSurface = true;
              std::cout << "[SweepHexMesher] block " << bi
                        << " surface annulus Q-Morph sweep: quads="
                        << qd.quads.size() << " tris=" << qd.tris.size()
                        << " holeEdges=" << patch.holeEdges.size() << " -> hex="
                        << hexMesh.hexes.size()
                        << " wedge=" << hexMesh.wedges.size()
                        << " ifaceAx=" << (tubeIface ? tubeIface->axialLower : 0)
                        << "\n";
            }
          }
        }
      }

      bool usedImprint = false;
      if (!usedSurface && surfaceMesh) {
        ImprintedCapMesh cap;
        if (SweepFaceImprinter::ImprintTranslationalCap(
                surfaceMesh, block, frame, nu, nv, nw, holeR, cap, false)) {
          bool meshed = false;
          if (cfg.useGmsh && GmshCapMesher::IsAvailable()) {
            const std::string workDir =
                cfg.gmshWorkDir + "/block_" + std::to_string(bi);
            meshed =
                GmshCapMesher::MeshTranslationalCap(cap, frame, holeR, workDir);
          }
          if (!meshed) {
            SweepFaceImprinter::BuildLocalQuadsFromSplits(frame, holeR, cap);
          }
          if (!cap.quads.empty()) {
            SweepFaceImprinter::SweepCapToHex(
                cap, frame, bi,
                "TranslationalSweepHex" + std::to_string(bi), hexMesh);
            usedImprint = true;
            std::cout << "[SweepHexMesher] block " << bi
                      << (meshed ? " gmsh" : " local")
                      << "+imprint+sweep hex: " << cap.quads.size()
                      << " quads x " << cap.sweepLayers << " layers -> "
                      << hexMesh.hexes.size() << " elements\n";
          }
        }
      }
      if (!usedSurface && !usedImprint) {
        OrientedBox box;
        box.origin = frame.origin;
        box.axis = frame.axis;
        box.crossY = frame.crossY;
        box.crossZ = frame.crossZ;
        box.axMin = frame.axMin;
        box.axMax = frame.axMax;
        box.crossMinY = frame.crossMinY;
        box.crossMaxY = frame.crossMaxY;
        box.crossMinZ = frame.crossMinZ;
        box.crossMaxZ = frame.crossMaxZ;
        if (cfg.targetCellSize > 0.0f) {
          nu = AutoDiv(box.crossMaxY - box.crossMinY, cfg.targetCellSize, nu);
          nv = AutoDiv(box.crossMaxZ - box.crossMinZ, cfg.targetCellSize, nv);
          nw = AutoDiv(box.axMax - box.axMin, cfg.targetCellSize, nw);
        }
        nu = ClampDiv(nu);
        nv = ClampDiv(nv);
        nw = ClampDiv(nw);
        hexMesh = GenerateTranslationalMesh(block, box, nu, nv, nw);
        std::cout << "[SweepHexMesher] block " << bi
                  << " translational hex (no imprint): " << nu << "x" << nv
                  << "x" << nw << " -> " << hexMesh.hexes.size() << " elements\n";
      }
    }

    hexMesh.blockIndex = bi;
    if (!hexMesh.hexes.empty() || !hexMesh.wedges.empty()) {
      meshes.push_back(std::move(hexMesh));
    }
  }
  return meshes;
}

bool SweepHexMesher::WriteVTK(const std::vector<SweepHexMesh> &meshes,
                              const std::string &path) {
  size_t totalNodes = 0;
  size_t totalHexes = 0;
  size_t totalWedges = 0;
  for (const auto &mesh : meshes) {
    totalNodes += mesh.nodes.size();
    totalHexes += mesh.hexes.size();
    totalWedges += mesh.wedges.size();
  }
  const size_t totalCells = totalHexes + totalWedges;
  if (totalCells == 0) {
    std::cerr << "[SweepHexMesher] no volume elements to write\n";
    return false;
  }

  std::ofstream out(path);
  if (!out.is_open()) {
    std::cerr << "[SweepHexMesher] cannot open " << path << "\n";
    return false;
  }

  out << "# vtk DataFile Version 3.0\n";
  out << "DFSweep swept hex+wedge mesh\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << totalNodes << " float\n";
  out << std::setprecision(9);
  for (const auto &mesh : meshes) {
    for (const auto &p : mesh.nodes) {
      out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
  }

  // HEX8: 9 ints/cell, WEDGE6: 7 ints/cell
  size_t listSize = totalHexes * 9 + totalWedges * 7;
  out << "CELLS " << totalCells << " " << listSize << "\n";
  int nodeOffset = 0;
  for (const auto &mesh : meshes) {
    for (const auto &hex : mesh.hexes) {
      out << "8";
      for (int nid : hex) {
        out << " " << (nodeOffset + nid);
      }
      out << "\n";
    }
    for (const auto &w : mesh.wedges) {
      out << "6";
      for (int nid : w) {
        out << " " << (nodeOffset + nid);
      }
      out << "\n";
    }
    nodeOffset += static_cast<int>(mesh.nodes.size());
  }

  out << "CELL_TYPES " << totalCells << "\n";
  for (const auto &mesh : meshes) {
    for (size_t i = 0; i < mesh.hexes.size(); ++i) {
      out << "12\n"; // VTK_HEXAHEDRON
    }
    for (size_t i = 0; i < mesh.wedges.size(); ++i) {
      out << "13\n"; // VTK_WEDGE
    }
  }

  out << "CELL_DATA " << totalCells << "\n";
  out << "SCALARS block_id int 1\n";
  out << "LOOKUP_TABLE default\n";
  for (const auto &mesh : meshes) {
    for (size_t i = 0; i < mesh.hexes.size(); ++i) {
      out << mesh.blockIndex << "\n";
    }
    for (size_t i = 0; i < mesh.wedges.size(); ++i) {
      out << mesh.blockIndex << "\n";
    }
  }

  std::cout << "[SweepHexMesher] wrote hex=" << totalHexes
            << " wedge=" << totalWedges << " to " << path << "\n";
  return true;
}

void SweepHexMesher::BuildWireframe(
    const SweepHexMesh &mesh, std::vector<Eigen::Vector3f> &vertices,
    std::vector<std::array<size_t, 2>> &edges) {
  static const std::array<std::array<int, 2>, 12> kHexEdges = {
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
       {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}}};
  static const std::array<std::array<int, 2>, 9> kWedgeEdges = {
      {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}}};

  vertices.clear();
  edges.clear();
  for (const auto &hex : mesh.hexes) {
    size_t base = vertices.size();
    for (int local : hex) {
      vertices.push_back(mesh.nodes[static_cast<size_t>(local)]);
    }
    for (const auto &e : kHexEdges) {
      edges.push_back({base + static_cast<size_t>(e[0]),
                       base + static_cast<size_t>(e[1])});
    }
  }
  for (const auto &w : mesh.wedges) {
    size_t base = vertices.size();
    for (int local : w) {
      vertices.push_back(mesh.nodes[static_cast<size_t>(local)]);
    }
    for (const auto &e : kWedgeEdges) {
      edges.push_back({base + static_cast<size_t>(e[0]),
                       base + static_cast<size_t>(e[1])});
    }
  }
}
