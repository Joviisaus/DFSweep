#include "SweepSurfaceMesher.h"
#include "Mesh/iterators.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>

namespace {
constexpr float kPi = 3.14159265358979323846f;

Eigen::Vector3f FaceCentroid(MeshLib::CToolFace *f) {
  Eigen::Vector3f c = Eigen::Vector3f::Zero();
  int n = 0;
  for (MeshLib::CTMesh::FaceVertexIterator fv(f); !fv.end(); ++fv) {
    auto *v = static_cast<MeshLib::CToolVertex *>(fv.value());
    c += Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
    ++n;
  }
  if (n > 0) {
    c /= static_cast<float>(n);
  }
  return c;
}

Eigen::Vector3f FaceNormal(MeshLib::CToolFace *f) {
  auto *he = f->halfedge();
  if (!he) {
    return Eigen::Vector3f::UnitY();
  }
  CPoint p1 = he->target()->point() - he->source()->point();
  CPoint p2 =
      he->he_next()->target()->point() - he->he_next()->source()->point();
  CPoint n = p1 ^ p2;
  if (n.norm() < 1e-12) {
    return Eigen::Vector3f::UnitY();
  }
  n /= n.norm();
  return Eigen::Vector3f(static_cast<float>(n[0]), static_cast<float>(n[1]),
                         static_cast<float>(n[2]));
}

void BuildCrossFrame(const Eigen::Vector3f &axisDir, Eigen::Vector3f &crossY,
                     Eigen::Vector3f &crossZ) {
  Eigen::Vector3f axis = axisDir;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f arbitrary =
      (std::abs(axis.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axis.dot(arbitrary) * axis).normalized();
  crossZ = axis.cross(crossY);
  if (crossZ.norm() < 1e-8f) {
    crossZ = Eigen::Vector3f::UnitZ();
  } else {
    crossZ.normalize();
  }
}

float RadialOf(const Eigen::Vector3f &p, const Eigen::Vector3f &origin,
               const Eigen::Vector3f &axis) {
  Eigen::Vector3f rel = p - origin;
  return (rel - rel.dot(axis) * axis).norm();
}

float AxialOf(const Eigen::Vector3f &p, const Eigen::Vector3f &origin,
              const Eigen::Vector3f &axis) {
  return (p - origin).dot(axis);
}

float AxialOfFrame(const Eigen::Vector3f &p, const SweepCapFrame &frame) {
  return (p - frame.origin).dot(frame.axis);
}

float QuadQuality(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                  const Eigen::Vector3f &c, const Eigen::Vector3f &d) {
  Eigen::Vector3f e0 = b - a;
  Eigen::Vector3f e1 = c - b;
  Eigen::Vector3f e2 = d - c;
  Eigen::Vector3f e3 = a - d;
  float l0 = e0.norm(), l1 = e1.norm(), l2 = e2.norm(), l3 = e3.norm();
  if (l0 < 1e-12f || l1 < 1e-12f || l2 < 1e-12f || l3 < 1e-12f) {
    return 0.0f;
  }
  float lMin = std::min(std::min(l0, l1), std::min(l2, l3));
  float lMax = std::max(std::max(l0, l1), std::max(l2, l3));
  float aspect = lMin / lMax;

  // 两对角剖分的法向必须同向，否则是蝴蝶结 / 自交
  Eigen::Vector3f nAbc = e0.cross(c - a);
  Eigen::Vector3f nAcd = (c - a).cross(d - a);
  Eigen::Vector3f nAbd = e0.cross(d - a);
  Eigen::Vector3f nBcd = (c - b).cross(d - b);
  if (nAbc.norm() < 1e-14f || nAcd.norm() < 1e-14f || nAbd.norm() < 1e-14f ||
      nBcd.norm() < 1e-14f) {
    return 0.0f;
  }
  if (nAbc.dot(nAcd) <= 0.0f || nAbd.dot(nBcd) <= 0.0f) {
    return 0.0f; // bowtie 或严重非凸
  }

  Eigen::Vector3f n1 = nAbc.normalized();
  Eigen::Vector3f n2 = nAcd.normalized();
  float planar = std::max(0.0f, n1.dot(n2));
  if (planar < 0.05f) {
    return 0.0f;
  }

  // 连续边转角：相对平均法向，四次转弯应同号（凸四边形）
  Eigen::Vector3f nAvg = (n1 + n2).normalized();
  auto turn = [&](const Eigen::Vector3f &u, const Eigen::Vector3f &v) {
    return nAvg.dot(u.cross(v));
  };
  float t0 = turn(e0, e1);
  float t1 = turn(e1, e2);
  float t2 = turn(e2, e3);
  float t3 = turn(e3, e0);
  if (t0 * t1 <= 0.0f || t1 * t2 <= 0.0f || t2 * t3 <= 0.0f ||
      t3 * t0 <= 0.0f) {
    return 0.0f;
  }

  Eigen::Vector3f mAc = 0.5f * (a + c);
  Eigen::Vector3f mBd = 0.5f * (b + d);
  float diagSpan = 0.5f * ((a - c).norm() + (b - d).norm());
  float midGap = (mAc - mBd).norm() / std::max(1e-6f, diagSpan);
  // 凸四边形对角线中点应接近；过大则非平面/扭曲
  float convex = std::max(0.0f, 1.0f - 2.0f * midGap);
  return aspect * planar * convex;
}

/// 将四边形整理为沿边界连续的凸绕序；失败返回 false
bool MakeConvexQuadOrder(int vOppA, int e0, int vOppB, int e1,
                         const std::vector<Eigen::Vector3f> &nodes,
                         std::array<int, 4> &outOrder, float &outQuality) {
  std::array<std::array<int, 4>, 2> orders = {
      {{{vOppA, e0, vOppB, e1}}, {{vOppA, e1, vOppB, e0}}}};
  outQuality = 0.0f;
  bool ok = false;
  for (const auto &ord : orders) {
    float q = QuadQuality(nodes[static_cast<size_t>(ord[0])],
                          nodes[static_cast<size_t>(ord[1])],
                          nodes[static_cast<size_t>(ord[2])],
                          nodes[static_cast<size_t>(ord[3])]);
    if (q > outQuality) {
      outQuality = q;
      outOrder = ord;
      ok = true;
    }
  }
  return ok && outQuality > 0.0f;
}

/// 使四边形绕序与参考方向一致（右手系，法向与 refDir 同向）
void OrientQuadCCW(std::array<int, 4> &q,
                   const std::vector<Eigen::Vector3f> &nodes,
                   const Eigen::Vector3f &refDir) {
  const Eigen::Vector3f &a = nodes[static_cast<size_t>(q[0])];
  const Eigen::Vector3f &b = nodes[static_cast<size_t>(q[1])];
  const Eigen::Vector3f &c = nodes[static_cast<size_t>(q[2])];
  Eigen::Vector3f n = (b - a).cross(c - a);
  if (n.dot(refDir) < 0.0f) {
    std::swap(q[1], q[3]);
  }
}

void OrientTriCCW(std::array<int, 3> &t,
                  const std::vector<Eigen::Vector3f> &nodes,
                  const Eigen::Vector3f &refDir) {
  const Eigen::Vector3f &a = nodes[static_cast<size_t>(t[0])];
  const Eigen::Vector3f &b = nodes[static_cast<size_t>(t[1])];
  const Eigen::Vector3f &c = nodes[static_cast<size_t>(t[2])];
  Eigen::Vector3f n = (b - a).cross(c - a);
  if (n.dot(refDir) < 0.0f) {
    std::swap(t[1], t[2]);
  }
}

/// HEX8 体积符号：底面 0-1-2-3，顶面 4-5-6-7；要求正体积
float HexSignedVolume(const std::array<Eigen::Vector3f, 8> &p) {
  Eigen::Vector3f c = Eigen::Vector3f::Zero();
  for (const auto &v : p) {
    c += v;
  }
  c /= 8.0f;
  auto tetVol = [&](const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                    const Eigen::Vector3f &c3, const Eigen::Vector3f &d) {
    return (b - a).dot((c3 - a).cross(d - a));
  };
  // 分解为 5/6 个四面体近似
  return tetVol(p[0], p[1], p[3], p[4]) + tetVol(p[1], p[2], p[3], p[6]) +
         tetVol(p[1], p[3], p[4], p[6]) + tetVol(p[1], p[4], p[5], p[6]) +
         tetVol(p[3], p[4], p[6], p[7]) + tetVol(p[4], p[5], p[6], p[7]);
}

std::array<int, 8> MakeHexFromSweep(int b0a, int b0b, int b0c, int b0d, int b1a,
                                    int b1b, int b1c, int b1d,
                                    const std::vector<Eigen::Vector3f> &nodes) {
  std::array<int, 8> h = {b0a, b0b, b0c, b0d, b1a, b1b, b1c, b1d};
  std::array<Eigen::Vector3f, 8> p;
  for (int k = 0; k < 8; ++k) {
    p[static_cast<size_t>(k)] = nodes[static_cast<size_t>(h[static_cast<size_t>(k)])];
  }
  if (HexSignedVolume(p) < 0.0f) {
    // 翻转底/顶面绕序
    h = {b0a, b0d, b0c, b0b, b1a, b1d, b1c, b1b};
  }
  return h;
}

std::array<int, 6> MakeWedgeFromSweep(int b0a, int b0b, int b0c, int b1a,
                                      int b1b, int b1c,
                                      const std::vector<Eigen::Vector3f> &nodes) {
  // VTK_WEDGE: 底 0,1,2 顶 3,4,5
  std::array<int, 6> w = {b0a, b0b, b0c, b1a, b1b, b1c};
  Eigen::Vector3f n0 =
      (nodes[static_cast<size_t>(b0b)] - nodes[static_cast<size_t>(b0a)])
          .cross(nodes[static_cast<size_t>(b0c)] -
                 nodes[static_cast<size_t>(b0a)]);
  Eigen::Vector3f up =
      nodes[static_cast<size_t>(b1a)] - nodes[static_cast<size_t>(b0a)];
  if (n0.dot(up) < 0.0f) {
    w = {b0a, b0c, b0b, b1a, b1c, b1b};
  }
  return w;
}

bool ProjectAlongAxisToSurface(const Eigen::Vector3f &p,
                               const Eigen::Vector3f &axis,
                               MeshLib::CTMesh *mesh, int targetBlock,
                               int targetType, Eigen::Vector3f &out) {
  if (!mesh) {
    return false;
  }
  Eigen::Vector3f best = p;
  float bestDist = std::numeric_limits<float>::max();
  bool found = false;
  for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
    auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
    if (targetBlock >= 0 && f->sweeplabel() != targetBlock) {
      continue;
    }
    if (f->sweepFaceType() != targetType) {
      continue;
    }
    Eigen::Vector3f c = FaceCentroid(f);
    Eigen::Vector3f delta = c - p;
    float dAx = std::abs(delta.dot(axis));
    Eigen::Vector3f radial = delta - delta.dot(axis) * axis;
    float dR = radial.norm();
    float score = dR + 0.15f * dAx;
    if (score < bestDist) {
      bestDist = score;
      best = p + delta.dot(axis) * axis;
      if (dR < 0.75f) {
        best = 0.75f * best + 0.25f * c;
      }
      found = true;
    }
  }
  if (found) {
    out = best;
  }
  return found;
}

bool ProjectRadiallyToSurface(const Eigen::Vector3f &p,
                              const Eigen::Vector3f &origin,
                              const Eigen::Vector3f &axis, float targetR,
                              MeshLib::CTMesh *mesh, int targetBlock,
                              int targetType, Eigen::Vector3f &out) {
  Eigen::Vector3f rel = p - origin;
  float ax = rel.dot(axis);
  Eigen::Vector3f rad = rel - ax * axis;
  if (rad.norm() < 1e-8f) {
    return false;
  }
  // 严格沿原射线改半径，绝不改角度（否则外层节点会串到错误方位）
  Eigen::Vector3f dir = rad.normalized();
  Eigen::Vector3f guess = origin + ax * axis + dir * targetR;

  if (!mesh) {
    out = guess;
    return true;
  }

  // 仅用同方位角附近的外壁面片微调半径，轴向/角度保持
  float bestR = targetR;
  float bestScore = std::numeric_limits<float>::max();
  bool found = false;
  for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
    auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
    if (targetBlock >= 0 && f->sweeplabel() != targetBlock) {
      continue;
    }
    if (targetType >= 0 && f->sweepFaceType() != targetType) {
      continue;
    }
    Eigen::Vector3f c = FaceCentroid(f);
    Eigen::Vector3f cRel = c - origin;
    float cAx = cRel.dot(axis);
    Eigen::Vector3f cRad = cRel - cAx * axis;
    if (cRad.norm() < 1e-8f) {
      continue;
    }
    float dAx = std::abs(cAx - ax);
    float ang = 1.0f - std::max(-1.0f, std::min(1.0f, dir.dot(cRad.normalized())));
    // 角度差过大则忽略（ang=0 同向，ang=2 反向）
    if (ang > 0.15f) { // ~cos^{-1}(0.85) ≈ 32°
      continue;
    }
    float score = dAx + 2.0f * ang;
    if (score < bestScore) {
      bestScore = score;
      bestR = cRad.norm();
      found = true;
    }
  }
  out = origin + ax * axis + dir * (found ? bestR : targetR);
  // 半径钳制，防止吸附到底面外缘等离群面片
  float rOut = (out - origin - ax * axis).norm();
  float rMax = std::max(targetR * 1.05f, targetR + 0.05f);
  float rMin = std::max(0.0f, targetR * 0.85f);
  if (rOut > rMax || rOut < rMin) {
    out = guess;
  }
  return true;
}

/// 在柱面参数平面 (θ, z) 内按绕质心逆时针排序，得到连续凸绕序
void OrderQuadByCylindricalParam(std::array<int, 4> &q,
                                 const std::vector<Eigen::Vector3f> &nodes,
                                 const Eigen::Vector3f &origin,
                                 const Eigen::Vector3f &axis) {
  struct PV {
    int id;
    float theta;
    float z;
    float ang;
  };
  std::array<PV, 4> pv{};
  float th0 = 0.0f;
  for (int k = 0; k < 4; ++k) {
    const Eigen::Vector3f &p = nodes[static_cast<size_t>(q[static_cast<size_t>(k)])];
    Eigen::Vector3f rel = p - origin;
    float z = rel.dot(axis);
    Eigen::Vector3f rad = rel - z * axis;
    // 用与轴垂直的任意正交基展开 θ
    Eigen::Vector3f e1 =
        (std::abs(axis.dot(Eigen::Vector3f::UnitX())) < 0.9f)
            ? Eigen::Vector3f::UnitX()
            : Eigen::Vector3f::UnitY();
    e1 = (e1 - e1.dot(axis) * axis).normalized();
    Eigen::Vector3f e2 = axis.cross(e1);
    float theta = std::atan2(rad.dot(e2), rad.dot(e1));
    pv[static_cast<size_t>(k)] = {q[static_cast<size_t>(k)], theta, z, 0.0f};
    if (k == 0) {
      th0 = theta;
    }
  }
  // 展开 θ 避免跨 -π/π 跳变
  for (int k = 1; k < 4; ++k) {
    float d = pv[static_cast<size_t>(k)].theta - th0;
    while (d > kPi) {
      d -= 2.0f * kPi;
    }
    while (d < -kPi) {
      d += 2.0f * kPi;
    }
    pv[static_cast<size_t>(k)].theta = th0 + d;
  }
  float cTh = 0.0f, cZ = 0.0f;
  for (int k = 0; k < 4; ++k) {
    cTh += pv[static_cast<size_t>(k)].theta;
    cZ += pv[static_cast<size_t>(k)].z;
  }
  cTh *= 0.25f;
  cZ *= 0.25f;
  for (int k = 0; k < 4; ++k) {
    pv[static_cast<size_t>(k)].ang = std::atan2(
        pv[static_cast<size_t>(k)].z - cZ, pv[static_cast<size_t>(k)].theta - cTh);
  }
  std::sort(pv.begin(), pv.end(),
            [](const PV &a, const PV &b) { return a.ang < b.ang; });
  for (int k = 0; k < 4; ++k) {
    q[static_cast<size_t>(k)] = pv[static_cast<size_t>(k)].id;
  }

  // 3D 法向应与局部外径向同向
  const Eigen::Vector3f &a = nodes[static_cast<size_t>(q[0])];
  const Eigen::Vector3f &b = nodes[static_cast<size_t>(q[1])];
  const Eigen::Vector3f &c = nodes[static_cast<size_t>(q[2])];
  Eigen::Vector3f n = (b - a).cross(c - a);
  Eigen::Vector3f cen = 0.25f * (nodes[static_cast<size_t>(q[0])] +
                                 nodes[static_cast<size_t>(q[1])] +
                                 nodes[static_cast<size_t>(q[2])] +
                                 nodes[static_cast<size_t>(q[3])]);
  Eigen::Vector3f rad = (cen - origin) - (cen - origin).dot(axis) * axis;
  if (rad.norm() > 1e-8f && n.dot(rad) < 0.0f) {
    std::swap(q[1], q[3]);
  }
}

void OrderTriByCylindricalParam(std::array<int, 3> &t,
                                const std::vector<Eigen::Vector3f> &nodes,
                                const Eigen::Vector3f &origin,
                                const Eigen::Vector3f &axis) {
  const Eigen::Vector3f &a = nodes[static_cast<size_t>(t[0])];
  const Eigen::Vector3f &b = nodes[static_cast<size_t>(t[1])];
  const Eigen::Vector3f &c = nodes[static_cast<size_t>(t[2])];
  Eigen::Vector3f n = (b - a).cross(c - a);
  Eigen::Vector3f cen = (a + b + c) / 3.0f;
  Eigen::Vector3f rad = (cen - origin) - (cen - origin).dot(axis) * axis;
  if (rad.norm() > 1e-8f && n.dot(rad) < 0.0f) {
    std::swap(t[1], t[2]);
  }
}

void CollectBoundaryAndHoleEdges(SurfaceCapPatch &out,
                                 const Eigen::Vector3f &origin,
                                 const Eigen::Vector3f &axis) {
  std::map<std::pair<int, int>, int> edgeCount;
  auto addEdge = [&](int a, int b) {
    if (a > b) {
      std::swap(a, b);
    }
    edgeCount[{a, b}]++;
  };
  for (const auto &t : out.tris) {
    addEdge(t[0], t[1]);
    addEdge(t[1], t[2]);
    addEdge(t[2], t[0]);
  }

  std::vector<std::pair<int, int>> boundary;
  for (const auto &[e, cnt] : edgeCount) {
    if (cnt == 1) {
      boundary.push_back(e);
    }
  }
  if (boundary.empty()) {
    return;
  }

  // 按边界边中点半径聚类：较小半径 = 镂空内环
  std::vector<float> midR;
  midR.reserve(boundary.size());
  float rMin = std::numeric_limits<float>::max();
  float rMax = 0.0f;
  for (const auto &e : boundary) {
    Eigen::Vector3f mid =
        0.5f * (out.nodes[static_cast<size_t>(e.first)] +
                out.nodes[static_cast<size_t>(e.second)]);
    float r = RadialOf(mid, origin, axis);
    midR.push_back(r);
    rMin = std::min(rMin, r);
    rMax = std::max(rMax, r);
  }
  float split = 0.5f * (rMin + rMax);
  // 若内外半径差很小，则全部当外边界
  bool haveHole = (rMax - rMin) > 0.15f * std::max(rMax, 1e-3f);
  out.holeRadius = haveHole ? rMin : 0.0f;
  out.outerRadius = rMax;

  for (size_t i = 0; i < boundary.size(); ++i) {
    const auto &e = boundary[i];
    std::array<Eigen::Vector3f, 2> seg = {
        out.nodes[static_cast<size_t>(e.first)],
        out.nodes[static_cast<size_t>(e.second)]};
    if (haveHole && midR[i] < split) {
      out.holeEdges.push_back(seg);
    } else {
      out.boundaryEdges.push_back(seg);
    }
  }
}

bool AppendFaceTris(MeshLib::CToolFace *f, SurfaceCapPatch &out,
                    std::map<MeshLib::CToolVertex *, int> &vertIndex) {
  auto getOrAdd = [&](MeshLib::CToolVertex *v) {
    auto it = vertIndex.find(v);
    if (it != vertIndex.end()) {
      return it->second;
    }
    int id = static_cast<int>(out.nodes.size());
    out.nodes.emplace_back(v->point()[0], v->point()[1], v->point()[2]);
    vertIndex[v] = id;
    return id;
  };
  std::vector<int> ids;
  for (MeshLib::CTMesh::FaceVertexIterator fv(f); !fv.end(); ++fv) {
    ids.push_back(getOrAdd(static_cast<MeshLib::CToolVertex *>(fv.value())));
  }
  if (ids.size() == 3) {
    out.tris.push_back({ids[0], ids[1], ids[2]});
    return true;
  }
  if (ids.size() == 4) {
    out.tris.push_back({ids[0], ids[1], ids[2]});
    out.tris.push_back({ids[0], ids[2], ids[3]});
    return true;
  }
  return false;
}
} // namespace

int SweepSurfaceMesher::MarkSourceFacesOnSurface(
    MeshLib::CTMesh *mesh, const std::vector<SweepBlockRegion> &blocks,
    const std::vector<bool> &blockNonPlanar) {
  if (!mesh) {
    return 0;
  }

  int markedBottom = 0;
  int markedInner = 0;
  int markedOuter = 0;

  for (int bi = 0; bi < static_cast<int>(blocks.size()); ++bi) {
    const auto &b = blocks[static_cast<size_t>(bi)];
    Eigen::Vector3f axis = b.sweepAxis;
    if (axis.norm() < 1e-8f) {
      axis = Eigen::Vector3f::UnitY();
    } else {
      axis.normalize();
    }
    Eigen::Vector3f origin = b.sweepOrigin;
    bool isCyl = (bi < static_cast<int>(blockNonPlanar.size()) &&
                  blockNonPlanar[static_cast<size_t>(bi)]) ||
                 b.kind == SweepKind::CylindricalBase;

    if (isCyl) {
      float rInner = std::max(1e-4f, b.radialInner);
      float rOuter = std::max(rInner + 1e-4f, b.radialOuter);
      float split = 0.5f * (rInner + rOuter);

      for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
        auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
        if (f->sweeplabel() != bi) {
          continue;
        }
        Eigen::Vector3f c = FaceCentroid(f);
        Eigen::Vector3f n = FaceNormal(f);
        if (std::abs(n.dot(axis)) > 0.4f) {
          continue; // 只要侧壁
        }
        float r = RadialOf(c, origin, axis);
        float a = AxialOf(c, origin, axis);
        if (b.axialUpper > b.axialLower + 1e-4f) {
          if (a < b.axialLower - 0.25f || a > b.axialUpper + 0.25f) {
            continue;
          }
        }
        // 按相对内外半径分类；靠近内半径 → 源侧壁(6)
        float dIn = std::abs(r - rInner);
        float dOut = std::abs(r - rOuter);
        if (dIn <= dOut || r <= split) {
          f->sweepFaceType() = 6;
          ++markedInner;
        } else {
          f->sweepFaceType() = 7;
          ++markedOuter;
        }
      }
    } else {
      // 平移块：底面圆环（法向对齐轴，靠近轴向最底）
      float axMin = std::numeric_limits<float>::max();
      float axMax = std::numeric_limits<float>::lowest();
      for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
        auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
        if (f->sweeplabel() != bi) {
          continue;
        }
        float a = AxialOf(FaceCentroid(f), origin, axis);
        axMin = std::min(axMin, a);
        axMax = std::max(axMax, a);
      }
      float span = std::max(1e-4f, axMax - axMin);
      float bottomBand = std::max(0.12f * span, 0.08f);
      float topBand = std::max(0.12f * span, 0.08f);

      for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
        auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
        if (f->sweeplabel() != bi) {
          continue;
        }
        Eigen::Vector3f n = FaceNormal(f);
        Eigen::Vector3f c = FaceCentroid(f);
        float align = std::abs(n.dot(axis));
        float a = AxialOf(c, origin, axis);
        bool nearBottom = a <= axMin + bottomBand;
        bool nearTop = a >= axMax - topBand;
        bool planarLike = (align > 0.7f) || (f->sweepFaceType() == 1) ||
                          (f->sweepFaceType() == 3);
        if (nearBottom && planarLike) {
          f->sweepFaceType() = 5;
          ++markedBottom;
        } else if (nearTop && planarLike) {
          f->sweepFaceType() = 3; // 目标顶面（含镂空上沿）
        }
      }
    }
  }

  std::cout << "[SweepSurfaceMesher] marked source faces: bottom(5)="
            << markedBottom << " innerWall(6)=" << markedInner
            << " outerWall(7)=" << markedOuter << "\n";
  return markedBottom + markedInner + markedOuter;
}

bool SweepSurfaceMesher::BuildFrameFromMeshFaces(MeshLib::CTMesh *mesh,
                                                 int blockIndex,
                                                 const SweepBlockRegion &block,
                                                 SweepCapFrame &frame) {
  if (!mesh) {
    return false;
  }
  Eigen::Vector3f axis = block.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f origin = block.sweepOrigin;
  Eigen::Vector3f crossY, crossZ;
  BuildCrossFrame(axis, crossY, crossZ);

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float minY = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::lowest();
  float minZ = std::numeric_limits<float>::max();
  float maxZ = std::numeric_limits<float>::lowest();
  int count = 0;

  for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
    auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
    if (f->sweeplabel() != blockIndex) {
      continue;
    }
    int ft = f->sweepFaceType();
    // 用源面 + 目标面共同定框
    if (ft != 5 && ft != 6 && ft != 7 && ft != 1 && ft != 3) {
      continue;
    }
    Eigen::Vector3f c = FaceCentroid(f);
    float a = AxialOf(c, origin, axis);
    float cy = (c - origin).dot(crossY);
    float cz = (c - origin).dot(crossZ);
    axMin = std::min(axMin, a);
    axMax = std::max(axMax, a);
    minY = std::min(minY, cy);
    maxY = std::max(maxY, cy);
    minZ = std::min(minZ, cz);
    maxZ = std::max(maxZ, cz);
    ++count;
  }

  if (count == 0) {
    // 回退：该块全部面片
    for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
      auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
      if (f->sweeplabel() != blockIndex) {
        continue;
      }
      Eigen::Vector3f c = FaceCentroid(f);
      float a = AxialOf(c, origin, axis);
      float cy = (c - origin).dot(crossY);
      float cz = (c - origin).dot(crossZ);
      axMin = std::min(axMin, a);
      axMax = std::max(axMax, a);
      minY = std::min(minY, cy);
      maxY = std::max(maxY, cy);
      minZ = std::min(minZ, cz);
      maxZ = std::max(maxZ, cz);
      ++count;
    }
  }
  if (count == 0 || axMax <= axMin + 1e-6f) {
    return false;
  }

  frame.origin = origin;
  frame.axis = axis;
  frame.crossY = crossY;
  frame.crossZ = crossZ;
  frame.axMin = axMin;
  frame.axMax = axMax;
  frame.crossMinY = minY;
  frame.crossMaxY = maxY;
  frame.crossMinZ = minZ;
  frame.crossMaxZ = maxZ;

  // 平移块底面很薄：扫掠高度优先用块 axial 范围（已含柱面下沿接口）
  if (block.kind == SweepKind::Translational) {
    if (block.axialUpper > block.axialLower + 1e-4f) {
      frame.axMin = std::min(frame.axMin, block.axialLower);
      frame.axMax = std::max(frame.axMax, block.axialUpper);
    }
    // 源面在底：axMin 取 type-5 面片轴向均值附近
    float srcAxSum = 0.0f;
    int srcCnt = 0;
    for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
      auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
      if (f->sweeplabel() != blockIndex || f->sweepFaceType() != 5) {
        continue;
      }
      srcAxSum += AxialOf(FaceCentroid(f), origin, axis);
      ++srcCnt;
    }
    if (srcCnt > 0) {
      float srcAx = srcAxSum / static_cast<float>(srcCnt);
      frame.axMin = srcAx;
      if (frame.axMax < frame.axMin + 1e-4f) {
        frame.axMax = frame.axMin + 0.2f;
      }
    } else if (frame.axMax - frame.axMin < 0.15f) {
      float rMax = 0.0f;
      for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
        auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
        if (f->sweeplabel() != blockIndex || f->sweepFaceType() != 5) {
          continue;
        }
        rMax = std::max(rMax, RadialOf(FaceCentroid(f), origin, axis));
      }
      float thick = std::max(0.2f, 0.25f * rMax);
      frame.axMax = frame.axMin + thick;
    }
  }

  std::cout << "[SweepSurfaceMesher] mesh frame block " << blockIndex
            << " ax=[" << frame.axMin << "," << frame.axMax << "] faces=" << count
            << "\n";
  return true;
}

bool SweepSurfaceMesher::ExtractBottomCapPatch(MeshLib::CTMesh *mesh,
                                               int blockIndex,
                                               const SweepCapFrame &frame,
                                               float holeRadius,
                                               SurfaceCapPatch &out) {
  out = SurfaceCapPatch{};
  out.blockIndex = blockIndex;
  out.holeRadius = holeRadius;
  if (!mesh) {
    return false;
  }

  std::map<MeshLib::CToolVertex *, int> vertIndex;
  int faceCount = 0;
  const float axTol =
      0.4f * std::max(1e-3f, frame.axMax - frame.axMin) + 0.2f;

  for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
    auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
    if (f->sweeplabel() != blockIndex) {
      continue;
    }
    if (f->sweepFaceType() != 5 && f->sweepFaceType() != 1) {
      continue;
    }
    Eigen::Vector3f c = FaceCentroid(f);
    if (std::abs(AxialOfFrame(c, frame) - frame.axMin) > axTol &&
        f->sweepFaceType() != 5) {
      continue;
    }
    // 镂空内部：中心半径小于 holeRadius 的面不要
    float r = RadialOf(c, frame.origin, frame.axis);
    if (holeRadius > 1e-4f && r < holeRadius * 0.85f) {
      continue;
    }
    if (AppendFaceTris(f, out, vertIndex)) {
      out.faceIds.push_back(faceCount++);
    }
  }

  CollectBoundaryAndHoleEdges(out, frame.origin, frame.axis);
  if (out.holeRadius <= 1e-4f && holeRadius > 1e-4f) {
    out.holeRadius = holeRadius;
  }

  std::cout << "[SweepSurfaceMesher] block " << blockIndex
            << " bottom annulus: nodes=" << out.nodes.size()
            << " tris=" << out.tris.size()
            << " outerEdges=" << out.boundaryEdges.size()
            << " holeEdges=" << out.holeEdges.size()
            << " holeR=" << out.holeRadius << "\n";
  return !out.tris.empty();
}

bool SweepSurfaceMesher::ExtractCylindricalWallPatch(
    MeshLib::CTMesh *mesh, int blockIndex, const SweepBlockRegion &block,
    bool innerWall, SurfaceCapPatch &out) {
  out = SurfaceCapPatch{};
  out.blockIndex = blockIndex;
  if (!mesh) {
    return false;
  }

  Eigen::Vector3f axis = block.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f origin = block.sweepOrigin;
  const int wantType = innerWall ? 6 : 7;

  // 严格限制在管壁环带，避免底面外缘/法兰顶点混入造成“尖刺”
  const float rLo = block.radialInner - 0.08f * std::max(1e-3f, block.radialOuter);
  const float rHi = block.radialOuter + 0.08f * std::max(1e-3f, block.radialOuter);
  const float axLo = block.axialLower - 0.05f;
  const float axHi = block.axialUpper + 0.05f;
  const float midR = 0.5f * (block.radialInner + block.radialOuter);

  std::map<MeshLib::CToolVertex *, int> vertIndex;
  int faceCount = 0;
  float rSum = 0.0f;
  int rCnt = 0;
  int skipped = 0;

  for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
    auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
    if (f->sweeplabel() != blockIndex) {
      continue;
    }
    Eigen::Vector3f c = FaceCentroid(f);
    float r = RadialOf(c, origin, axis);
    float a = AxialOf(c, origin, axis);
    if (r < rLo || r > rHi || a < axLo || a > axHi) {
      ++skipped;
      continue;
    }

    bool accept = (f->sweepFaceType() == wantType);
    if (!accept) {
      // 回退：侧壁几何（法向⊥轴）且半径靠近内/外
      Eigen::Vector3f n = FaceNormal(f);
      if (std::abs(n.dot(axis)) > 0.35f) {
        continue;
      }
      if (innerWall && r > midR) {
        continue;
      }
      if (!innerWall && r < midR) {
        continue;
      }
      accept = true;
    }
    if (!accept) {
      continue;
    }

    // 再检查面片顶点：任一顶点远超环带则丢弃（防止共点把底面外缘拉进来）
    bool vertOk = true;
    for (MeshLib::CTMesh::FaceVertexIterator fv(f); !fv.end(); ++fv) {
      auto *v = static_cast<MeshLib::CToolVertex *>(fv.value());
      Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
      float vr = RadialOf(p, origin, axis);
      if (vr > rHi * 1.15f || vr < rLo * 0.75f) {
        vertOk = false;
        break;
      }
    }
    if (!vertOk) {
      ++skipped;
      continue;
    }

    rSum += r;
    ++rCnt;
    if (AppendFaceTris(f, out, vertIndex)) {
      out.faceIds.push_back(faceCount++);
    }
  }

  CollectBoundaryAndHoleEdges(out, origin, axis);
  if (rCnt > 0) {
    float avgR = rSum / static_cast<float>(rCnt);
    if (innerWall) {
      out.holeRadius = avgR;
    } else {
      out.outerRadius = avgR;
    }
  }

  std::cout << "[SweepSurfaceMesher] block " << blockIndex
            << (innerWall ? " inner" : " outer")
            << " wall patch: nodes=" << out.nodes.size()
            << " tris=" << out.tris.size() << " skipped=" << skipped
            << " rBand=[" << rLo << "," << rHi << "]\n";
  return !out.tris.empty();
}

bool SweepSurfaceMesher::BuildQuadDominantFromTris(const SurfaceCapPatch &patch,
                                                   float minQuadQuality,
                                                   QuadDominantCap &out) {
  out = QuadDominantCap{};
  out.nodes = patch.nodes;
  out.imprintEdges = patch.boundaryEdges;
  out.imprintEdges.insert(out.imprintEdges.end(), patch.holeEdges.begin(),
                          patch.holeEdges.end());
  if (patch.tris.empty()) {
    return false;
  }

  // 参考法向：用前几个三角形平均法向，保证后续绕序一致
  Eigen::Vector3f refN = Eigen::Vector3f::Zero();
  int nRef = 0;
  for (size_t i = 0; i < std::min<size_t>(patch.tris.size(), 32); ++i) {
    const auto &t = patch.tris[i];
    Eigen::Vector3f a = out.nodes[static_cast<size_t>(t[0])];
    Eigen::Vector3f b = out.nodes[static_cast<size_t>(t[1])];
    Eigen::Vector3f c = out.nodes[static_cast<size_t>(t[2])];
    Eigen::Vector3f n = (b - a).cross(c - a);
    if (n.norm() > 1e-12f) {
      refN += n.normalized();
      ++nRef;
    }
  }
  if (nRef > 0) {
    refN.normalize();
  } else {
    refN = Eigen::Vector3f::UnitY();
  }

  const int nT = static_cast<int>(patch.tris.size());
  std::map<std::pair<int, int>, std::vector<int>> edgeToTris;
  auto key = [](int a, int b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
  };
  for (int ti = 0; ti < nT; ++ti) {
    const auto &t = patch.tris[static_cast<size_t>(ti)];
    edgeToTris[key(t[0], t[1])].push_back(ti);
    edgeToTris[key(t[1], t[2])].push_back(ti);
    edgeToTris[key(t[2], t[0])].push_back(ti);
  }

  struct Candidate {
    int t0, t1;
    int a, b, c, d;
    float quality;
  };
  std::vector<Candidate> cands;
  std::set<std::pair<int, int>> seenPair;

  for (const auto &[e, tris] : edgeToTris) {
    if (tris.size() != 2) {
      continue;
    }
    int t0 = tris[0], t1 = tris[1];
    auto pkey = key(t0, t1);
    if (seenPair.count(pkey)) {
      continue;
    }
    seenPair.insert(pkey);

    const auto &A = patch.tris[static_cast<size_t>(t0)];
    const auto &B = patch.tris[static_cast<size_t>(t1)];
    int vOppA = -1, vOppB = -1;
    for (int k = 0; k < 3; ++k) {
      if (A[static_cast<size_t>(k)] != e.first &&
          A[static_cast<size_t>(k)] != e.second) {
        vOppA = A[static_cast<size_t>(k)];
      }
      if (B[static_cast<size_t>(k)] != e.first &&
          B[static_cast<size_t>(k)] != e.second) {
        vOppB = B[static_cast<size_t>(k)];
      }
    }
    if (vOppA < 0 || vOppB < 0) {
      continue;
    }

    std::array<int, 4> bestOrder{};
    float bestQ = 0.0f;
    if (!MakeConvexQuadOrder(vOppA, e.first, vOppB, e.second, out.nodes,
                             bestOrder, bestQ) ||
        bestQ < minQuadQuality) {
      continue;
    }
    float sharedLen =
        (out.nodes[static_cast<size_t>(e.first)] -
         out.nodes[static_cast<size_t>(e.second)])
            .norm();
    float score = bestQ / std::max(1e-4f, sharedLen);
    cands.push_back({t0, t1, bestOrder[0], bestOrder[1], bestOrder[2],
                     bestOrder[3], score});
  }

  std::sort(cands.begin(), cands.end(),
            [](const Candidate &x, const Candidate &y) {
              return x.quality > y.quality;
            });

  std::vector<char> used(static_cast<size_t>(nT), 0);
  for (const auto &c : cands) {
    if (used[static_cast<size_t>(c.t0)] || used[static_cast<size_t>(c.t1)]) {
      continue;
    }
    used[static_cast<size_t>(c.t0)] = 1;
    used[static_cast<size_t>(c.t1)] = 1;
    out.quads.push_back({c.a, c.b, c.c, c.d});
  }

  // 不再做无质量检查的强制合并，避免蝴蝶结四边形
  for (int ti = 0; ti < nT; ++ti) {
    if (!used[static_cast<size_t>(ti)]) {
      out.tris.push_back(patch.tris[static_cast<size_t>(ti)]);
    }
  }

  // 镂空：去掉中心落在孔内的单元（用相对扫掠轴的径向距离）
  if (patch.holeRadius > 1e-4f && !patch.holeEdges.empty() &&
      !out.nodes.empty()) {
    // 用孔边界边估计轴原点：取孔边中点平均作为轴上参考点附近的圆心投影
    Eigen::Vector3f holeCen = Eigen::Vector3f::Zero();
    for (const auto &seg : patch.holeEdges) {
      holeCen += 0.5f * (seg[0] + seg[1]);
    }
    holeCen /= static_cast<float>(patch.holeEdges.size());

    auto radialFromHole = [&](const Eigen::Vector3f &p) {
      // 在底面近似平面内：到孔心的距离
      return (p - holeCen).norm();
    };
    auto keepQuad = [&](const std::array<int, 4> &q) {
      Eigen::Vector3f m = 0.25f * (out.nodes[static_cast<size_t>(q[0])] +
                                   out.nodes[static_cast<size_t>(q[1])] +
                                   out.nodes[static_cast<size_t>(q[2])] +
                                   out.nodes[static_cast<size_t>(q[3])]);
      return radialFromHole(m) >= patch.holeRadius * 0.92f;
    };
    auto keepTri = [&](const std::array<int, 3> &t) {
      Eigen::Vector3f m = (out.nodes[static_cast<size_t>(t[0])] +
                           out.nodes[static_cast<size_t>(t[1])] +
                           out.nodes[static_cast<size_t>(t[2])]) /
                          3.0f;
      return radialFromHole(m) >= patch.holeRadius * 0.92f;
    };
    std::vector<std::array<int, 4>> q2;
    for (const auto &q : out.quads) {
      if (keepQuad(q)) {
        q2.push_back(q);
      }
    }
    out.quads.swap(q2);
    std::vector<std::array<int, 3>> t2;
    for (const auto &t : out.tris) {
      if (keepTri(t)) {
        t2.push_back(t);
      }
    }
    out.tris.swap(t2);
  }

  const float quadRatio =
      nT > 0 ? (2.0f * static_cast<float>(out.quads.size()) /
                static_cast<float>(nT))
             : 0.0f;

  // 不在这里做全局绕序：柱面整圈平均法向会抵消，导致大量翻转。
  // 绕序由后续 Sweep*ToVolume 按局部扫掠方向（轴向/局部径向）统一处理。
  (void)refN;

  std::cout << "[SweepSurfaceMesher] Q-Morph-like: quads=" << out.quads.size()
            << " tris=" << out.tris.size() << " (from " << nT
            << " tris, quad-cover≈" << quadRatio << ")\n";
  return !out.quads.empty() || !out.tris.empty();
}

bool SweepSurfaceMesher::SweepQuadDominantToVolume(
    const QuadDominantCap &cap, const SweepCapFrame &frame,
    MeshLib::CTMesh *mesh, int blockIndex, const std::string &name, int layers,
    SweepHexMesh &out, const SweepBlockRegion *interfaceTube) {
  out = SweepHexMesh{};
  out.blockIndex = blockIndex;
  out.name = name;
  out.kind = SweepKind::Translational;
  out.imprintEdges = cap.imprintEdges;

  const int nw = std::max(1, layers);
  const int nCap = static_cast<int>(cap.nodes.size());
  if (nCap <= 0 || (cap.quads.empty() && cap.tris.empty())) {
    return false;
  }

  Eigen::Vector3f axis = frame.axis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }

  // 目标轴向：优先贴到柱面下沿接口，实现几何拼合
  float axTarget = frame.axMax;
  Eigen::Vector3f tubeOrigin = frame.origin;
  float tubeRInner = 0.0f;
  float tubeROuter = 0.0f;
  if (interfaceTube) {
    axTarget = interfaceTube->axialLower;
    tubeOrigin = interfaceTube->sweepOrigin;
    tubeRInner = interfaceTube->radialInner;
    tubeROuter = interfaceTube->radialOuter;
    // 保证扫掠方向从底面指向接口
    if ((axTarget - frame.axMin) * (frame.axMax - frame.axMin) < 0.0f) {
      // 若轴向约定相反，取管体上沿
      axTarget = interfaceTube->axialUpper;
    }
  }

  out.nodes.reserve(static_cast<size_t>(nCap * (nw + 1)));
  for (int iw = 0; iw <= nw; ++iw) {
    float t = static_cast<float>(iw) / static_cast<float>(nw);
    float ax = frame.axMin + t * (axTarget - frame.axMin);
    for (int i = 0; i < nCap; ++i) {
      const Eigen::Vector3f &p0 = cap.nodes[static_cast<size_t>(i)];
      float a0 = AxialOfFrame(p0, frame);
      Eigen::Vector3f p = p0 + (ax - a0) * axis;

      if (iw == nw) {
        // 顶层：轴向贴到柱面下沿；孔内节点贴内半径；外环保持原径向（不夹到 rOuter）
        if (interfaceTube && tubeROuter > tubeRInner + 1e-6f) {
          Eigen::Vector3f rel = p0 - tubeOrigin;
          float aRel = rel.dot(axis);
          Eigen::Vector3f rad = rel - aRel * axis;
          float r = rad.norm();
          if (r > 1e-8f) {
            Eigen::Vector3f dir = rad.normalized();
            float rKeep = r;
            // 孔内：贴到管体内半径，与柱壁下口内缘拼合
            if (r < tubeRInner) {
              rKeep = tubeRInner;
            }
            p = tubeOrigin + axTarget * axis + dir * rKeep;
          } else {
            p = tubeOrigin + axTarget * axis +
                Eigen::Vector3f::UnitX() * tubeRInner;
          }
        } else if (mesh) {
          Eigen::Vector3f snapped;
          if (ProjectAlongAxisToSurface(p, axis, mesh, blockIndex, 3,
                                        snapped)) {
            p = snapped;
          }
        }
      }
      out.nodes.push_back(p);
    }
  }

  // 截面绕序相对扫掠方向
  Eigen::Vector3f sweepDir = ((axTarget - frame.axMin) >= 0.0f) ? axis : -axis;
  auto orientCap = [&](std::array<int, 4> q) {
    OrientQuadCCW(q, cap.nodes, sweepDir);
    return q;
  };
  auto orientTri = [&](std::array<int, 3> t) {
    OrientTriCCW(t, cap.nodes, sweepDir);
    return t;
  };

  auto base = [&](int iw) { return iw * nCap; };
  int flipped = 0;
  for (int iw = 0; iw < nw; ++iw) {
    int b0 = base(iw);
    int b1 = base(iw + 1);
    for (const auto &qin : cap.quads) {
      auto q = orientCap(qin);
      auto hex = MakeHexFromSweep(b0 + q[0], b0 + q[1], b0 + q[2], b0 + q[3],
                                  b1 + q[0], b1 + q[1], b1 + q[2], b1 + q[3],
                                  out.nodes);
      // 统计是否发生了翻转修正
      if (hex[1] == b0 + q[3]) {
        ++flipped;
      }
      out.hexes.push_back(hex);
    }
    for (const auto &trin : cap.tris) {
      auto tr = orientTri(trin);
      out.wedges.push_back(MakeWedgeFromSweep(
          b0 + tr[0], b0 + tr[1], b0 + tr[2], b1 + tr[0], b1 + tr[1],
          b1 + tr[2], out.nodes));
    }
  }

  std::cout << "[SweepSurfaceMesher] translational sweep: hex="
            << out.hexes.size() << " wedge=" << out.wedges.size()
            << " ax=[" << frame.axMin << "," << axTarget << "]"
            << " flipped=" << flipped << "\n";
  return !out.hexes.empty() || !out.wedges.empty();
}

bool SweepSurfaceMesher::SweepCylindricalWallToVolume(
    const QuadDominantCap &cap, const SweepBlockRegion &block,
    MeshLib::CTMesh *mesh, int blockIndex, const std::string &name,
    int radialLayers, SweepHexMesh &out) {
  out = SweepHexMesh{};
  out.blockIndex = blockIndex;
  out.name = name;
  out.kind = SweepKind::CylindricalBase;
  out.imprintEdges = cap.imprintEdges;

  const int nr = std::max(1, radialLayers);
  const int nCap = static_cast<int>(cap.nodes.size());
  if (nCap <= 0 || (cap.quads.empty() && cap.tris.empty())) {
    return false;
  }

  Eigen::Vector3f axis = block.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f origin = block.sweepOrigin;

  // 每个源节点的局部半径与单位径向（避免全局平均径向在整圈上抵消）
  std::vector<float> nodeR(static_cast<size_t>(nCap), 0.0f);
  std::vector<Eigen::Vector3f> nodeRadDir(static_cast<size_t>(nCap),
                                          Eigen::Vector3f::UnitX());
  std::vector<float> nodeAx(static_cast<size_t>(nCap), 0.0f);
  float r0Sum = 0.0f;
  int r0Cnt = 0;
  const float rBandLo = block.radialInner * 0.85f;
  const float rBandHi = block.radialOuter * 1.08f;
  int clampedSrc = 0;
  for (int i = 0; i < nCap; ++i) {
    const Eigen::Vector3f &p0 = cap.nodes[static_cast<size_t>(i)];
    float ax = AxialOf(p0, origin, axis);
    Eigen::Vector3f rad = (p0 - origin) - ax * axis;
    float r = rad.norm();
    nodeAx[static_cast<size_t>(i)] = ax;
    if (r > 1e-8f) {
      nodeRadDir[static_cast<size_t>(i)] = rad / r;
    }
    // 钳制离群源点（底面外缘误入），避免径向扫掠长出尖刺
    if (r > rBandHi) {
      r = block.radialInner;
      ++clampedSrc;
    } else if (r < rBandLo && block.radialInner > 1e-4f) {
      r = block.radialInner;
      ++clampedSrc;
    }
    nodeR[static_cast<size_t>(i)] = r;
    r0Sum += r;
    ++r0Cnt;
  }
  float r0 = (r0Cnt > 0) ? (r0Sum / static_cast<float>(r0Cnt)) : block.radialInner;
  float r1 = std::max(r0 + 1e-4f, block.radialOuter);
  if (mesh) {
    float rOuterSum = 0.0f;
    int rOuterCnt = 0;
    for (MeshLib::MeshFaceIterator mf(mesh); !mf.end(); ++mf) {
      auto *f = static_cast<MeshLib::CToolFace *>(mf.value());
      if (f->sweeplabel() != blockIndex || f->sweepFaceType() != 7) {
        continue;
      }
      float rr = RadialOf(FaceCentroid(f), origin, axis);
      if (rr < rBandLo || rr > rBandHi) {
        continue;
      }
      rOuterSum += rr;
      ++rOuterCnt;
    }
    if (rOuterCnt > 0) {
      r1 = std::max(r0 + 1e-4f, rOuterSum / static_cast<float>(rOuterCnt));
    }
  }
  r1 = std::min(r1, block.radialOuter * 1.05f);

  auto localOutward = [&](const std::array<int, 4> &q) -> Eigen::Vector3f {
    Eigen::Vector3f c = Eigen::Vector3f::Zero();
    for (int k = 0; k < 4; ++k) {
      c += cap.nodes[static_cast<size_t>(q[static_cast<size_t>(k)])];
    }
    c *= 0.25f;
    Eigen::Vector3f rad = (c - origin) - AxialOf(c, origin, axis) * axis;
    if (rad.norm() > 1e-8f) {
      return rad.normalized();
    }
    return Eigen::Vector3f::UnitX();
  };
  (void)localOutward;

  // 预定向：在 (θ,z) 参数平面排序，法向对齐局部外径向
  std::vector<std::array<int, 4>> quads;
  std::vector<std::array<int, 3>> tris = cap.tris;
  int splitBad = 0;
  quads.reserve(cap.quads.size());
  for (auto q : cap.quads) {
    OrderQuadByCylindricalParam(q, cap.nodes, origin, axis);
    float qQual =
        QuadQuality(cap.nodes[static_cast<size_t>(q[0])],
                    cap.nodes[static_cast<size_t>(q[1])],
                    cap.nodes[static_cast<size_t>(q[2])],
                    cap.nodes[static_cast<size_t>(q[3])]);
    if (qQual < 0.05f) {
      // 不合格四边形拆成两个三角形，避免蝴蝶结 hex
      tris.push_back({q[0], q[1], q[2]});
      tris.push_back({q[0], q[2], q[3]});
      ++splitBad;
      continue;
    }
    quads.push_back(q);
  }
  for (auto &t : tris) {
    OrderTriByCylindricalParam(t, cap.nodes, origin, axis);
  }

  out.nodes.reserve(static_cast<size_t>(nCap * (nr + 1)));
  for (int ir = 0; ir <= nr; ++ir) {
    float t = static_cast<float>(ir) / static_cast<float>(nr);
    for (int i = 0; i < nCap; ++i) {
      const Eigen::Vector3f &p0 = cap.nodes[static_cast<size_t>(i)];
      float ax = nodeAx[static_cast<size_t>(i)];
      float rSrc = nodeR[static_cast<size_t>(i)];
      // 每个节点按自身半径比例扫到 r1，保持相对厚度均匀
      float r =
          (rSrc > 1e-8f) ? (rSrc + t * (r1 - rSrc)) : (r0 + t * (r1 - r0));
      Eigen::Vector3f p;
      if (ir == 0) {
        // 源层：离群点拉回内壁半径，其余贴原表面
        if (clampedSrc > 0 && nodeR[static_cast<size_t>(i)] <= block.radialInner + 1e-4f &&
            RadialOf(p0, origin, axis) > rBandHi) {
          p = origin + ax * axis +
              nodeRadDir[static_cast<size_t>(i)] * nodeR[static_cast<size_t>(i)];
        } else {
          p = p0;
        }
      } else {
        p = origin + ax * axis + nodeRadDir[static_cast<size_t>(i)] * r;
      }
      // 外层：只沿原射线贴合目标半径，不改角度
      if (ir == nr) {
        Eigen::Vector3f snapped;
        if (ProjectRadiallyToSurface(p, origin, axis, r1, mesh, blockIndex, 7,
                                     snapped)) {
          p = snapped;
        }
      }
      out.nodes.push_back(p);
    }
  }

  auto base = [&](int ir) { return ir * nCap; };
  int flipped = 0;
  int rejected = 0;
  for (int ir = 0; ir < nr; ++ir) {
    int b0 = base(ir);
    int b1 = base(ir + 1);
    for (const auto &q : quads) {
      auto hex = MakeHexFromSweep(b0 + q[0], b0 + q[1], b0 + q[2], b0 + q[3],
                                  b1 + q[0], b1 + q[1], b1 + q[2], b1 + q[3],
                                  out.nodes);
      if (hex[1] == b0 + q[3]) {
        ++flipped;
      }
      std::array<Eigen::Vector3f, 8> pv;
      for (int k = 0; k < 8; ++k) {
        pv[static_cast<size_t>(k)] =
            out.nodes[static_cast<size_t>(hex[static_cast<size_t>(k)])];
      }
      if (HexSignedVolume(pv) <= 1e-12f) {
        ++rejected;
        continue;
      }
      // 侧棱应近似径向：角向漂移过大则丢弃（错误连接）
      bool badMap = false;
      for (int k = 0; k < 4; ++k) {
        const Eigen::Vector3f &pi =
            out.nodes[static_cast<size_t>(hex[static_cast<size_t>(k)])];
        const Eigen::Vector3f &po =
            out.nodes[static_cast<size_t>(hex[static_cast<size_t>(k + 4)])];
        Eigen::Vector3f ri = (pi - origin) - AxialOf(pi, origin, axis) * axis;
        Eigen::Vector3f ro = (po - origin) - AxialOf(po, origin, axis) * axis;
        if (ri.norm() > 1e-8f && ro.norm() > 1e-8f) {
          float cosang = ri.normalized().dot(ro.normalized());
          if (cosang < 0.85f) {
            badMap = true;
            break;
          }
        }
      }
      if (badMap) {
        ++rejected;
        continue;
      }
      out.hexes.push_back(hex);
    }
    for (const auto &tr : tris) {
      out.wedges.push_back(MakeWedgeFromSweep(
          b0 + tr[0], b0 + tr[1], b0 + tr[2], b1 + tr[0], b1 + tr[1],
          b1 + tr[2], out.nodes));
    }
  }

  std::cout << "[SweepSurfaceMesher] cylindrical wall sweep: r=[" << r0 << ","
            << r1 << "] hex=" << out.hexes.size()
            << " wedge=" << out.wedges.size() << " volFlip=" << flipped
            << " rejected=" << rejected << " splitBadQuad=" << splitBad
            << " clampedSrc=" << clampedSrc << "\n";
  return !out.hexes.empty() || !out.wedges.empty();
}
