#include "SweepFaceImprinter.h"
#include "SweepHexMesher.h"
#include "Mesh/iterators.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>

namespace {
constexpr float kPi = 3.14159265358979323846f;

void BuildCrossFrame(const Eigen::Vector3f &axisDir, Eigen::Vector3f &crossY,
                     Eigen::Vector3f &crossZ) {
  Eigen::Vector3f arbitrary =
      (std::abs(axisDir.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axisDir.dot(arbitrary) * axisDir).normalized();
  crossZ = axisDir.cross(crossY).normalized();
}

Eigen::Vector2f ToCapUV(const Eigen::Vector3f &p, const SweepCapFrame &frame) {
  Eigen::Vector3f rel = p - frame.origin;
  return {rel.dot(frame.crossY), rel.dot(frame.crossZ)};
}

Eigen::Vector3f FromCapUV(float u, float v, float ax,
                          const SweepCapFrame &frame) {
  return frame.origin + ax * frame.axis + u * frame.crossY + v * frame.crossZ;
}

float AxialOf(const Eigen::Vector3f &p, const SweepCapFrame &frame) {
  return (p - frame.origin).dot(frame.axis);
}

void MergeSplit(std::vector<float> &splits, float value, float tol) {
  for (float existing : splits) {
    if (std::abs(existing - value) <= tol) {
      return;
    }
  }
  splits.push_back(value);
}

void SortUniqueSplits(std::vector<float> &splits, float tol) {
  std::sort(splits.begin(), splits.end());
  std::vector<float> merged;
  for (float v : splits) {
    if (merged.empty() || std::abs(v - merged.back()) > tol) {
      merged.push_back(v);
    }
  }
  splits.swap(merged);
}

bool OnCap(const Eigen::Vector3f &p, const SweepCapFrame &frame, float capAx,
           float tol) {
  return std::abs(AxialOf(p, frame) - capAx) <= tol;
}

Eigen::Vector3f ProjectToBottomCap(const Eigen::Vector3f &p,
                                   const SweepCapFrame &frame) {
  float ax = AxialOf(p, frame);
  return p - (ax - frame.axMin) * frame.axis;
}

bool IsFeatureEdgeOnCap(MeshLib::CToolEdge *edge, float capTol,
                        const SweepCapFrame &frame, float capAx) {
  if (!edge || !edge->halfedge(0)) {
    return false;
  }
  auto *he0 = static_cast<MeshLib::CToolHalfEdge *>(edge->halfedge(0));
  auto *he1 = static_cast<MeshLib::CToolHalfEdge *>(edge->halfedge(1));
  auto *v0 = static_cast<MeshLib::CToolVertex *>(he0->source());
  auto *v1 = static_cast<MeshLib::CToolVertex *>(he0->target());
  Eigen::Vector3f p0(v0->point()[0], v0->point()[1], v0->point()[2]);
  Eigen::Vector3f p1(v1->point()[0], v1->point()[1], v1->point()[2]);
  if (!OnCap(p0, frame, capAx, capTol) || !OnCap(p1, frame, capAx, capTol)) {
    return false;
  }
  if (edge->boundary() || edge->sharp()) {
    return true;
  }
  if (v0->FeaturePoint() && v1->FeaturePoint()) {
    return true;
  }
  if (he0->face() && he1 && he1->face()) {
    auto *f0 = static_cast<MeshLib::CToolFace *>(he0->face());
    auto *f1 = static_cast<MeshLib::CToolFace *>(he1->face());
    if (f0->label() != f1->label() &&
        (v0->FeaturePoint() || v1->FeaturePoint())) {
      return true;
    }
  }
  return false;
}

void ExtractCapFeatureSegments(MeshLib::CTMesh *mesh, const SweepCapFrame &frame,
                               float capAx, float capTol, float minSegLen,
                               std::vector<std::array<Eigen::Vector2f, 2>>
                                   &segments2d,
                               std::vector<std::array<Eigen::Vector3f, 2>>
                                   &segments3d) {
  if (!mesh) {
    return;
  }
  for (MeshLib::MeshEdgeIterator eiter(mesh); !eiter.end(); ++eiter) {
    auto *edge = static_cast<MeshLib::CToolEdge *>(eiter.value());
    if (!IsFeatureEdgeOnCap(edge, capTol, frame, capAx)) {
      continue;
    }
    auto *he = static_cast<MeshLib::CToolHalfEdge *>(edge->halfedge(0));
    auto *v0 = static_cast<MeshLib::CToolVertex *>(he->source());
    auto *v1 = static_cast<MeshLib::CToolVertex *>(he->target());
    Eigen::Vector3f p0(v0->point()[0], v0->point()[1], v0->point()[2]);
    Eigen::Vector3f p1(v1->point()[0], v1->point()[1], v1->point()[2]);
    float segLen = (p1 - p0).norm();
    bool cornerEdge = v0->FeaturePoint() && v1->FeaturePoint();
    if (!cornerEdge && segLen < minSegLen) {
      continue;
    }
    Eigen::Vector2f uv0 = ToCapUV(p0, frame);
    Eigen::Vector2f uv1 = ToCapUV(p1, frame);
    segments2d.push_back({uv0, uv1});
    segments3d.push_back({p0, p1});
  }
}

void LimitSplits(std::vector<float> &splits, int maxCells) {
  if (static_cast<int>(splits.size()) - 1 <= maxCells) {
    return;
  }
  float lo = splits.front();
  float hi = splits.back();
  splits.clear();
  splits.push_back(lo);
  for (int i = 1; i < maxCells; ++i) {
    float t = static_cast<float>(i) / static_cast<float>(maxCells);
    splits.push_back(lo + t * (hi - lo));
  }
  splits.push_back(hi);
}

void AddUniformInteriorSplits(std::vector<float> &splits, float lo, float hi,
                              int targetCount) {
  if (targetCount <= 1 || hi <= lo + 1e-6f) {
    return;
  }
  for (int i = 1; i < targetCount; ++i) {
    float t = static_cast<float>(i) / static_cast<float>(targetCount);
    splits.push_back(lo + t * (hi - lo));
  }
}

void ImprintSegmentsToSplits(
    const std::vector<std::array<Eigen::Vector2f, 2>> &segments,
    float tol, std::vector<float> &uSplits, std::vector<float> &vSplits) {
  for (const auto &seg : segments) {
    MergeSplit(uSplits, seg[0].x(), tol);
    MergeSplit(vSplits, seg[0].y(), tol);
    MergeSplit(uSplits, seg[1].x(), tol);
    MergeSplit(vSplits, seg[1].y(), tol);
  }
}

bool QuadCenterValid(float u0, float u1, float v0, float v1,
                     float holeRadius, const SweepCapFrame &frame) {
  if (holeRadius <= 0.0f) {
    return true;
  }
  float cu = 0.5f * (u0 + u1);
  float cv = 0.5f * (v0 + v1);
  Eigen::Vector3f c = FromCapUV(cu, cv, frame.axMin, frame);
  Eigen::Vector3f rel = c - frame.origin;
  float radial = (rel - rel.dot(frame.axis) * frame.axis).norm();
  return radial > holeRadius;
}
} // namespace

bool SweepFaceImprinter::BuildFrameFromHex(
    const std::map<int, Eigen::Vector3f> &hex,
    const Eigen::Vector3f &axisHint, const Eigen::Vector3f &originHint,
    SweepCapFrame &frame) {
  if (hex.size() != 8) {
    return false;
  }
  frame.origin = originHint;
  frame.axis = axisHint.normalized();
  BuildCrossFrame(frame.axis, frame.crossY, frame.crossZ);
  frame.axMin = std::numeric_limits<float>::max();
  frame.axMax = std::numeric_limits<float>::lowest();
  frame.crossMinY = std::numeric_limits<float>::max();
  frame.crossMaxY = std::numeric_limits<float>::lowest();
  frame.crossMinZ = std::numeric_limits<float>::max();
  frame.crossMaxZ = std::numeric_limits<float>::lowest();
  for (const auto &[_, p] : hex) {
    Eigen::Vector2f uv = ToCapUV(p, frame);
    float ax = AxialOf(p, frame);
    frame.axMin = std::min(frame.axMin, ax);
    frame.axMax = std::max(frame.axMax, ax);
    frame.crossMinY = std::min(frame.crossMinY, uv.x());
    frame.crossMaxY = std::max(frame.crossMaxY, uv.x());
    frame.crossMinZ = std::min(frame.crossMinZ, uv.y());
    frame.crossMaxZ = std::max(frame.crossMaxZ, uv.y());
  }
  return frame.axMax > frame.axMin + 1e-6f;
}

bool SweepFaceImprinter::ImprintTranslationalCap(
    MeshLib::CTMesh *mesh, const SweepBlockRegion &block,
    const SweepCapFrame &frame, int minDivU, int minDivV, int sweepLayers,
    float holeRadius, ImprintedCapMesh &out, bool buildLocalQuads) {
  out = ImprintedCapMesh{};
  out.sweepLayers = std::max(1, sweepLayers);

  const float capTol =
      std::max(0.01f, 0.05f * (frame.axMax - frame.axMin));
  const float splitTol =
      std::max(1e-3f, 0.02f * std::min(frame.crossMaxY - frame.crossMinY,
                                       frame.crossMaxZ - frame.crossMinZ));

  const float crossExtent = std::min(frame.crossMaxY - frame.crossMinY,
                                     frame.crossMaxZ - frame.crossMinZ);
  const float minSegLen = std::max(0.05f, 0.12f * crossExtent);

  std::vector<std::array<Eigen::Vector2f, 2>> bottomSeg2d;
  std::vector<std::array<Eigen::Vector3f, 2>> bottomSeg3d;
  ExtractCapFeatureSegments(mesh, frame, frame.axMin, capTol, minSegLen,
                            bottomSeg2d, bottomSeg3d);

  std::vector<std::array<Eigen::Vector2f, 2>> topSeg2d;
  std::vector<std::array<Eigen::Vector3f, 2>> topSeg3d;
  ExtractCapFeatureSegments(mesh, frame, frame.axMax, capTol, minSegLen,
                            topSeg2d, topSeg3d);

  std::vector<std::array<Eigen::Vector2f, 2>> imprintSeg2d = bottomSeg2d;
  std::vector<std::array<Eigen::Vector3f, 2>> imprintSeg3d = bottomSeg3d;
  for (size_t i = 0; i < topSeg3d.size(); ++i) {
    Eigen::Vector3f q0 = ProjectToBottomCap(topSeg3d[i][0], frame);
    Eigen::Vector3f q1 = ProjectToBottomCap(topSeg3d[i][1], frame);
    imprintSeg2d.push_back({ToCapUV(q0, frame), ToCapUV(q1, frame)});
    imprintSeg3d.push_back({q0, q1});
  }

  std::vector<float> uSplits = {frame.crossMinY, frame.crossMaxY};
  std::vector<float> vSplits = {frame.crossMinZ, frame.crossMaxZ};
  ImprintSegmentsToSplits(imprintSeg2d, splitTol, uSplits, vSplits);
  SortUniqueSplits(uSplits, splitTol);
  SortUniqueSplits(vSplits, splitTol);

  while (static_cast<int>(uSplits.size()) - 1 < minDivU) {
    AddUniformInteriorSplits(uSplits, frame.crossMinY, frame.crossMaxY,
                             minDivU + 1);
    SortUniqueSplits(uSplits, splitTol);
  }
  while (static_cast<int>(vSplits.size()) - 1 < minDivV) {
    AddUniformInteriorSplits(vSplits, frame.crossMinZ, frame.crossMaxZ,
                             minDivV + 1);
    SortUniqueSplits(vSplits, splitTol);
  }

  const int maxCells = std::max(minDivU, minDivV) + 8;
  LimitSplits(uSplits, maxCells);
  LimitSplits(vSplits, maxCells);
  SortUniqueSplits(uSplits, splitTol);
  SortUniqueSplits(vSplits, splitTol);

  const int nu = static_cast<int>(uSplits.size()) - 1;
  const int nv = static_cast<int>(vSplits.size()) - 1;
  out.uSplits = uSplits;
  out.vSplits = vSplits;
  out.imprintEdges = imprintSeg3d;

  if (buildLocalQuads) {
    BuildLocalQuadsFromSplits(frame, holeRadius, out);
  }

  std::cout << "[SweepFaceImprinter] imprint: bottomFeat=" << bottomSeg3d.size()
            << " topFeat=" << topSeg3d.size()
            << " imprinted=" << imprintSeg3d.size() << " grid=" << nu << "x"
            << nv << " quads=" << out.quads.size() << " layers="
            << out.sweepLayers << "\n";
  return !out.uSplits.empty() && !out.vSplits.empty();
}

void SweepFaceImprinter::BuildLocalQuadsFromSplits(const SweepCapFrame &frame,
                                                   float holeRadius,
                                                   ImprintedCapMesh &cap) {
  const auto &uSplits = cap.uSplits;
  const auto &vSplits = cap.vSplits;
  const int nu = static_cast<int>(uSplits.size()) - 1;
  const int nv = static_cast<int>(vSplits.size()) - 1;
  cap.quads.clear();
  cap.capNodes.clear();
  cap.capNodes.reserve(static_cast<size_t>((nu + 1) * (nv + 1)));
  for (int iv = 0; iv <= nv; ++iv) {
    for (int iu = 0; iu <= nu; ++iu) {
      cap.capNodes.push_back(
          FromCapUV(uSplits[static_cast<size_t>(iu)],
                    vSplits[static_cast<size_t>(iv)], frame.axMin, frame));
    }
  }
  auto nid = [&](int iu, int iv) { return iv * (nu + 1) + iu; };
  for (int iv = 0; iv < nv; ++iv) {
    for (int iu = 0; iu < nu; ++iu) {
      float u0 = uSplits[static_cast<size_t>(iu)];
      float u1 = uSplits[static_cast<size_t>(iu + 1)];
      float v0 = vSplits[static_cast<size_t>(iv)];
      float v1 = vSplits[static_cast<size_t>(iv + 1)];
      if (!QuadCenterValid(u0, u1, v0, v1, holeRadius, frame)) {
        continue;
      }
      cap.quads.push_back({nid(iu, iv), nid(iu + 1, iv), nid(iu + 1, iv + 1),
                           nid(iu, iv + 1)});
    }
  }
}

void SweepFaceImprinter::SweepCapToHex(const ImprintedCapMesh &cap,
                                       const SweepCapFrame &frame,
                                       int blockIndex, const std::string &name,
                                       SweepHexMesh &out) {
  out = SweepHexMesh{};
  out.blockIndex = blockIndex;
  out.name = name;
  out.kind = SweepKind::Translational;
  out.imprintEdges = cap.imprintEdges;

  const int nw = cap.sweepLayers;
  const int nCapNodes = static_cast<int>(cap.capNodes.size());
  if (nCapNodes <= 0 || cap.quads.empty()) {
    return;
  }

  out.nodes.reserve(static_cast<size_t>(nCapNodes * (nw + 1)));
  for (int iw = 0; iw <= nw; ++iw) {
    float t = static_cast<float>(iw) / static_cast<float>(nw);
    float ax = frame.axMin + t * (frame.axMax - frame.axMin);
    for (const auto &bottomNode : cap.capNodes) {
      Eigen::Vector2f uv = ToCapUV(bottomNode, frame);
      out.nodes.push_back(FromCapUV(uv.x(), uv.y(), ax, frame));
    }
  }

  auto layerBase = [&](int iw) { return iw * nCapNodes; };
  for (int iw = 0; iw < nw; ++iw) {
    int b0 = layerBase(iw);
    int b1 = layerBase(iw + 1);
    for (const auto &q : cap.quads) {
      std::array<int, 8> hex = {
          b0 + q[0], b0 + q[1], b0 + q[2], b0 + q[3],
          b1 + q[0], b1 + q[1], b1 + q[2], b1 + q[3]};
      out.hexes.push_back(hex);
    }
  }
}

void SweepFaceImprinter::SweepCylindricalCapToHex(
    const ImprintedCapMesh &innerCap, const SweepBlockRegion &block,
    int radialLayers, int blockIndex, const std::string &name,
    SweepHexMesh &out) {
  out = SweepHexMesh{};
  out.blockIndex = blockIndex;
  out.name = name;
  out.kind = SweepKind::CylindricalBase;
  out.imprintEdges = innerCap.imprintEdges;

  const int nCapNodes = static_cast<int>(innerCap.capNodes.size());
  const int nr = std::max(1, radialLayers);
  if (nCapNodes <= 0 || innerCap.quads.empty()) {
    return;
  }

  Eigen::Vector3f axis = block.sweepAxis.normalized();
  Eigen::Vector3f origin = block.sweepOrigin;
  const float r0 = std::max(0.0f, block.radialInner);
  const float r1 = std::max(r0 + 1e-4f, block.radialOuter);

  out.nodes.reserve(static_cast<size_t>(nCapNodes * (nr + 1)));
  for (int ir = 0; ir <= nr; ++ir) {
    float t = static_cast<float>(ir) / static_cast<float>(nr);
    float r = r0 + t * (r1 - r0);
    for (const auto &pInner : innerCap.capNodes) {
      Eigen::Vector3f rel = pInner - origin;
      float ax = rel.dot(axis);
      Eigen::Vector3f radial = rel - ax * axis;
      if (r0 > 1e-6f && radial.norm() > 1e-8f) {
        radial = radial.normalized() * r;
      } else {
        Eigen::Vector3f crossY = block.crossDirY;
        Eigen::Vector3f crossZ = block.crossDirZ;
        if (crossY.norm() < 1e-6f) {
          BuildCrossFrame(axis, crossY, crossZ);
        } else {
          crossY.normalize();
          crossZ = axis.cross(crossY).normalized();
        }
        radial = r * crossY;
      }
      out.nodes.push_back(origin + ax * axis + radial);
    }
  }

  auto layerBase = [&](int ir) { return ir * nCapNodes; };
  for (int ir = 0; ir < nr; ++ir) {
    int b0 = layerBase(ir);
    int b1 = layerBase(ir + 1);
    for (const auto &q : innerCap.quads) {
      std::array<int, 8> hex = {
          b0 + q[0], b0 + q[1], b0 + q[2], b0 + q[3],
          b1 + q[0], b1 + q[1], b1 + q[2], b1 + q[3]};
      out.hexes.push_back(hex);
    }
  }
}
