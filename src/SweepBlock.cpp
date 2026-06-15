#include "SweepBlock.h"
#include "Mesh/iterators.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

#ifdef ENABLE_OMP
#include <omp.h>
#endif

namespace {

const int kDx[] = {-1, 1, 0, 0, 0, 0};
const int kDy[] = {0, 0, -1, 1, 0, 0};
const int kDz[] = {0, 0, 0, 0, -1, 1};

struct IsoVoxel {
  float field;
  VoxelIndex idx;
};

struct IsoVoxelGreater {
  bool operator()(const IsoVoxel &a, const IsoVoxel &b) const {
    return a.field > b.field;
  }
};

Eigen::Vector3f FootOnAxis(const Eigen::Vector3f &p,
                           const Eigen::Vector3f &axisRef,
                           const Eigen::Vector3f &axisDir) {
  return axisRef + (p - axisRef).dot(axisDir) * axisDir;
}

Eigen::Vector3f RadialComponent(const Eigen::Vector3f &p,
                                const Eigen::Vector3f &axisRef,
                                const Eigen::Vector3f &axisDir) {
  return p - FootOnAxis(p, axisRef, axisDir);
}

void BuildCylinderCrossFrame(const Eigen::Vector3f &axisDir,
                             Eigen::Vector3f &crossY, Eigen::Vector3f &crossZ) {
  Eigen::Vector3f arbitrary =
      (std::abs(axisDir.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axisDir.dot(arbitrary) * axisDir).normalized();
  crossZ = axisDir.cross(crossY).normalized();
}

} // namespace

SweepDecomposer::SweepDecomposer(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<std::vector<float>>> &Field,
    const std::vector<std::vector<std::vector<int>>> &FieldLabel,
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField,
    const std::vector<PrimeData> &primes, float threshold,
    bool cylinderPairsOnly, MeshLib::CTMesh *mesh)
    : Coord(Coord), Field(Field), FieldLabel(FieldLabel), GradField(GradField),
      primes(primes), mesh(mesh), threshold(threshold),
      cylinderPairsOnly(cylinderPairsOnly) {

  D1 = static_cast<int>(Field.size());
  D2 = (D1 > 0) ? static_cast<int>(Field[0].size()) : 0;
  D3 = (D2 > 0) ? static_cast<int>(Field[0][0].size()) : 0;
  stepSize = (D1 > 0 && D2 > 0 && D3 > 0)
                 ? (Coord[0][0][0] - Coord[0][0][1]).norm()
                 : 0.1f;

  std::cout << "[SweepDecomposer] grid " << D1 << "x" << D2 << "x" << D3
            << " step=" << stepSize << " threshold=" << threshold << " rad\n";

  if (this->cylinderPairsOnly) {
    CollectInteriorVoxels();
    this->cylinderPairs = FindMatchingCylinderPairs();
    auto groups = BuildCylinderSweepGroups();
    for (const auto &group : groups) {
      std::cout << "[SweepDecomposer] cylinder sweep group:";
      for (int pid : group) {
        std::cout << " " << pid;
      }
      std::cout << " [energy-grown region]\n";
      GrowCylinderSweepRegion(group);
    }
    selectedBlocks = candidateBlocks;
    std::cout << "[SweepDecomposer] " << selectedBlocks.size()
              << " cylinder sweep block(s)\n";
    return;
  }

  CollectInteriorVoxels();

  std::set<int> pairedPatches;
  this->cylinderPairs = FindMatchingCylinderPairs();
  for (const auto &pair : this->cylinderPairs) {
    std::cout << "[SweepDecomposer] cylinder pair (gradient-matched): "
              << pair.first << " <-> " << pair.second << "\n";
    GrowBlockBetweenCylinderPair(pair.first, pair.second);
    pairedPatches.insert(pair.first);
    pairedPatches.insert(pair.second);
  }

  if (!this->cylinderPairsOnly) {
    for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
      if (primes[static_cast<size_t>(i)].params.size() < 10) {
        continue;
      }
      if (pairedPatches.count(i)) {
        continue;
      }
      if (IsCylinderPatch(i)) {
        continue;
      }
      GrowBlockFromPatch(i);
    }
  }

  std::cout << "[SweepDecomposer] " << candidateBlocks.size()
            << " candidates from " << primes.size() << " patches\n";

  GreedySetCover();

  std::cout << "[SweepDecomposer] selected " << selectedBlocks.size()
            << " blocks\n";
}

void SweepDecomposer::CollectInteriorVoxels() {
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] < 0.0f) {
          allInteriorVoxels.insert({x, y, z});
        }
      }
    }
  }
  std::cout << "[SweepDecomposer] interior voxels: "
            << allInteriorVoxels.size() << "\n";
}

bool SweepDecomposer::IsCylinderPatch(int patchId) const {
  if (patchId < 0 || patchId >= static_cast<int>(primes.size())) {
    return false;
  }
  const PrimeData &p = primes[static_cast<size_t>(patchId)];
  if (p.isPlane || p.params.size() < 10) {
    return false;
  }

  const auto &q = p.params;
  Eigen::Matrix3f hessian;
  hessian << static_cast<float>(2 * q[7]), static_cast<float>(q[4]),
      static_cast<float>(q[5]), static_cast<float>(q[4]),
      static_cast<float>(2 * q[8]), static_cast<float>(q[6]),
      static_cast<float>(q[5]), static_cast<float>(q[6]),
      static_cast<float>(2 * q[9]);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(hessian);
  const auto &evals = solver.eigenvalues();
  float absEvals[3] = {std::abs(evals[0]), std::abs(evals[1]), std::abs(evals[2])};
  float maxEv = std::max({absEvals[0], absEvals[1], absEvals[2]});
  if (maxEv < 1e-6f) {
    return false;
  }

  int nearZeroCount = 0;
  int largeCount = 0;
  for (float ev : absEvals) {
    if (ev / maxEv < 0.2f) {
      nearZeroCount++;
    }
    if (ev / maxEv > 0.5f) {
      largeCount++;
    }
  }

  // 柱面：Hessian 一特征值近零、其余两个同量级
  if (nearZeroCount >= 1 && largeCount >= 2) {
    return true;
  }

  double qxx = q[7], qyy = q[8], qzz = q[9];
  double cross = std::abs(q[4]) + std::abs(q[5]) + std::abs(q[6]);
  double spread =
      std::max({std::abs(qxx - qyy), std::abs(qyy - qzz), std::abs(qxx - qzz)});

  if (qxx > 1e-8 && qyy > 1e-8 && qzz > 1e-8 &&
      spread < 0.25 * (qxx + qyy + qzz)) {
    return false;
  }

  return cross > 1e-6 ||
         spread < 0.35 * (std::abs(qxx) + std::abs(qyy) + std::abs(qzz));
}

Eigen::Vector3f SweepDecomposer::ExtractCylinderAxis(int patchId) const {
  const auto &q = primes[static_cast<size_t>(patchId)].params;
  Eigen::Matrix3f hessian;
  hessian << static_cast<float>(2 * q[7]), static_cast<float>(q[4]),
      static_cast<float>(q[5]), static_cast<float>(q[4]),
      static_cast<float>(2 * q[8]), static_cast<float>(q[6]),
      static_cast<float>(q[5]), static_cast<float>(q[6]),
      static_cast<float>(2 * q[9]);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(hessian);
  const auto &evals = solver.eigenvalues();
  float maxAbs = std::max({std::abs(evals[0]), std::abs(evals[1]), std::abs(evals[2])});
  int minIdx = 0;
  float bestScore = std::numeric_limits<float>::max();
  for (int i = 0; i < 3; ++i) {
    float score = maxAbs > 1e-8f ? std::abs(evals[i]) / maxAbs : std::abs(evals[i]);
    if (score < bestScore) {
      bestScore = score;
      minIdx = i;
    }
  }

  Eigen::Vector3f axis = solver.eigenvectors().col(minIdx);
  if (axis.norm() < 1e-8f) {
    return Eigen::Vector3f::UnitY();
  }
  return axis.normalized();
}

Eigen::Vector3f SweepDecomposer::ComputePatchCentroid(int patchId) const {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f centroid = Eigen::Vector3f::Zero();
  int count = 0;
  const float surfaceBand = 2.5f * stepSize;

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (FieldLabel[x][y][z] == primeId &&
            std::abs(Field[x][y][z]) < surfaceBand) {
          centroid += Coord[x][y][z];
          count++;
        }
      }
    }
  }

  if (count > 0) {
    centroid /= static_cast<float>(count);
  }
  return centroid;
}

Eigen::Vector3f
SweepDecomposer::ComputePatchAvgGradient(int patchId, int *sampleCount) const {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f avg = Eigen::Vector3f::Zero();
  int count = 0;
  const float surfaceBand = 2.5f * stepSize;

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (FieldLabel[x][y][z] != primeId) {
          continue;
        }
        if (std::abs(Field[x][y][z]) >= surfaceBand) {
          continue;
        }
        const Eigen::Vector3f &g = GradField[x][y][z];
        if (g.norm() < 1e-8f) {
          continue;
        }
        avg += g.normalized();
        count++;
      }
    }
  }

  if (sampleCount) {
    *sampleCount = count;
  }
  if (count > 0) {
    avg /= static_cast<float>(count);
  }
  return avg;
}

Eigen::Vector3f
SweepDecomposer::ComputePatchCentroidFromMesh(int patchId) const {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f centroid = Eigen::Vector3f::Zero();
  int count = 0;
  if (!mesh) {
    return centroid;
  }

  for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
    auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
    if (v->label() != primeId || v->FeaturePoint()) {
      continue;
    }
    centroid += Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
    count++;
  }
  if (count > 0) {
    centroid /= static_cast<float>(count);
  }
  return centroid;
}

Eigen::Vector3f SweepDecomposer::EvalPrimeGradient(
    int patchId, const Eigen::Vector3f &p) const {
  const auto &q = primes[static_cast<size_t>(patchId)].params;
  return Eigen::Vector3f(
      static_cast<float>(q[1] + q[4] * p.y() + q[5] * p.z() + 2 * q[7] * p.x()),
      static_cast<float>(q[2] + q[4] * p.x() + q[6] * p.z() + 2 * q[8] * p.y()),
      static_cast<float>(q[3] + q[5] * p.x() + q[6] * p.y() + 2 * q[9] * p.z()));
}

Eigen::Vector3f SweepDecomposer::ComputePatchAvgAnalyticGradientFromMesh(
    int patchId, int *sampleCount) const {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f avg = Eigen::Vector3f::Zero();
  int count = 0;
  if (!mesh) {
    if (sampleCount) {
      *sampleCount = 0;
    }
    return avg;
  }

  for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
    auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
    if (v->label() != primeId || v->FeaturePoint()) {
      continue;
    }
    Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
    Eigen::Vector3f g = EvalPrimeGradient(patchId, p);
    if (g.norm() < 1e-8f) {
      continue;
    }
    avg += g.normalized();
    count++;
  }

  if (sampleCount) {
    *sampleCount = count;
  }
  if (count > 0) {
    avg /= static_cast<float>(count);
  }
  return avg;
}

Eigen::Vector3f
SweepDecomposer::ComputePatchAvgGradientFromMesh(int patchId,
                                                 int *sampleCount) const {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f avg = Eigen::Vector3f::Zero();
  int count = 0;
  if (!mesh) {
    if (sampleCount) {
      *sampleCount = 0;
    }
    return avg;
  }

  for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
    auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
    if (v->label() != primeId || v->FeaturePoint()) {
      continue;
    }
    Eigen::Vector3f n(v->normal()[0], v->normal()[1], v->normal()[2]);
    if (n.norm() < 1e-8f) {
      continue;
    }
    avg += n.normalized();
    count++;
  }

  if (sampleCount) {
    *sampleCount = count;
  }
  if (count > 0) {
    avg /= static_cast<float>(count);
  }
  return avg;
}

std::map<int, Eigen::Vector3f> SweepDecomposer::BuildHexFromPoints(
    const std::vector<Eigen::Vector3f> &points, const Eigen::Vector3f &axis,
    const Eigen::Vector3f &origin) const {
  Eigen::Vector3f dirX = axis.normalized();
  Eigen::Vector3f arbitrary =
      (std::abs(dirX.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  Eigen::Vector3f dirY =
      (arbitrary - dirX.dot(arbitrary) * dirX).normalized();
  Eigen::Vector3f dirZ = dirX.cross(dirY).normalized();

  float minX = std::numeric_limits<float>::max();
  float maxX = std::numeric_limits<float>::lowest();
  float minY = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::lowest();
  float minZ = std::numeric_limits<float>::max();
  float maxZ = std::numeric_limits<float>::lowest();

  const float margin = std::max(stepSize, 0.05f);
  for (const auto &p : points) {
    float px = (p - origin).dot(dirX);
    float py = (p - origin).dot(dirY);
    float pz = (p - origin).dot(dirZ);
    minX = std::min(minX, px);
    maxX = std::max(maxX, px);
    minY = std::min(minY, py);
    maxY = std::max(maxY, py);
    minZ = std::min(minZ, pz);
    maxZ = std::max(maxZ, pz);
  }

  minX -= margin;
  maxX += margin;
  minY -= margin;
  maxY += margin;
  minZ -= margin;
  maxZ += margin;

  Eigen::Matrix3f A;
  A.row(0) = dirX.transpose();
  A.row(1) = dirY.transpose();
  A.row(2) = dirZ.transpose();

  std::map<int, Eigen::Vector3f> vertices;
  float bx[] = {minX, maxX};
  float by[] = {minY, maxY};
  float bz[] = {minZ, maxZ};
  int idx = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        Eigen::Vector3f rhs = origin + bx[i] * dirX + by[j] * dirY + bz[k] * dirZ;
        vertices[idx++] = rhs;
      }
    }
  }
  return vertices;
}

std::map<int, Eigen::Vector3f> SweepDecomposer::BuildCylinderRadialHex(
    const std::vector<Eigen::Vector3f> &points, const Eigen::Vector3f &axisDir,
    const Eigen::Vector3f &axisOrigin, float rInner, float rOuter, float axMin,
    float axMax, const Eigen::Vector3f &crossY,
    const Eigen::Vector3f &crossZ) const {
  float minY = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::lowest();
  float minZ = std::numeric_limits<float>::max();
  float maxZ = std::numeric_limits<float>::lowest();
  for (const auto &p : points) {
    Eigen::Vector3f rad = RadialComponent(p, axisOrigin, axisDir);
    minY = std::min(minY, rad.dot(crossY));
    maxY = std::max(maxY, rad.dot(crossY));
    minZ = std::min(minZ, rad.dot(crossZ));
    maxZ = std::max(maxZ, rad.dot(crossZ));
  }

  const float margin = std::max(stepSize, 0.05f);
  axMin -= margin;
  axMax += margin;
  minY -= margin;
  maxY += margin;
  minZ -= margin;
  maxZ += margin;
  rOuter = std::max(rOuter + margin, rInner + margin);

  auto corner = [&](float ax, float cy, float cz) {
    Eigen::Vector3f foot = axisOrigin + ax * axisDir;
    return foot + cy * crossY + cz * crossZ;
  };

  std::map<int, Eigen::Vector3f> vertices;
  float axVals[] = {axMin, axMax};
  float yVals[] = {minY, maxY};
  float zVals[] = {minZ, maxZ};
  int idx = 0;
  for (float ax : axVals) {
    for (float cy : yVals) {
      for (float cz : zVals) {
        vertices[idx++] = corner(ax, cy, cz);
      }
    }
  }
  return vertices;
}

CylinderPairViz SweepDecomposer::BuildCylinderPairViz(
    int patchA, int patchB, const Eigen::Vector3f &axis,
    const Eigen::Vector3f &origin, const std::vector<Eigen::Vector3f> &ptsA,
    const std::vector<Eigen::Vector3f> &ptsB) const {
  CylinderPairViz viz;
  viz.patchIdA = patchA;
  viz.patchIdB = patchB;
  viz.primeIdA = primes[static_cast<size_t>(patchA)].id;
  viz.primeIdB = primes[static_cast<size_t>(patchB)].id;
  viz.sweepAxis = axis.normalized();
  viz.sweepOrigin = origin;
  BuildCylinderCrossFrame(viz.sweepAxis, viz.crossDirY, viz.crossDirZ);

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float rMax = 0.0f;
  float oMinY = std::numeric_limits<float>::max();
  float oMaxY = std::numeric_limits<float>::lowest();
  float oMinZ = std::numeric_limits<float>::max();
  float oMaxZ = std::numeric_limits<float>::lowest();

  auto accumulate = [&](const Eigen::Vector3f &p) {
    float ax = (p - origin).dot(viz.sweepAxis);
    Eigen::Vector3f rad = RadialComponent(p, origin, viz.sweepAxis);
    float r = rad.norm();
    axMin = std::min(axMin, ax);
    axMax = std::max(axMax, ax);
    rMax = std::max(rMax, r);
    oMinY = std::min(oMinY, rad.dot(viz.crossDirY));
    oMaxY = std::max(oMaxY, rad.dot(viz.crossDirY));
    oMinZ = std::min(oMinZ, rad.dot(viz.crossDirZ));
    oMaxZ = std::max(oMaxZ, rad.dot(viz.crossDirZ));
  };

  for (const auto &p : ptsA) {
    accumulate(p);
  }
  for (const auto &p : ptsB) {
    accumulate(p);
  }

  const float margin = std::max(0.02f, 0.5f * stepSize);
  viz.axialMid = 0.5f * (axMin + axMax);
  viz.radialOuter = rMax + margin;
  viz.radialInner = std::max(margin, 0.08f * viz.radialOuter);

  viz.outerCrossMinY = oMinY;
  viz.outerCrossMaxY = oMaxY;
  viz.outerCrossMinZ = oMinZ;
  viz.outerCrossMaxZ = oMaxZ;

  const float scale =
      viz.radialOuter > 1e-6f ? viz.radialInner / viz.radialOuter : 0.2f;
  viz.innerCrossMinY = oMinY * scale;
  viz.innerCrossMaxY = oMaxY * scale;
  viz.innerCrossMinZ = oMinZ * scale;
  viz.innerCrossMaxZ = oMaxZ * scale;
  return viz;
}

void SweepDecomposer::BuildCylinderSinglePatchBlockFromMesh(
    int patchId, int pairedPatchId, const Eigen::Vector3f &axis,
    const Eigen::Vector3f &origin, const std::vector<Eigen::Vector3f> &patchPts,
    const CylinderPairViz &pairViz) {
  if (patchPts.size() < 2) {
    return;
  }

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float rMax = 0.0f;
  for (const auto &p : patchPts) {
    float ax = (p - origin).dot(axis);
    axMin = std::min(axMin, ax);
    axMax = std::max(axMax, ax);
    rMax = std::max(rMax, RadialComponent(p, origin, axis).norm());
  }

  const float margin = std::max(0.02f, 0.5f * stepSize);
  float rInner = std::max(margin, 0.08f * (rMax + margin));
  float rOuter = rMax + margin;

  auto hex = BuildCylinderRadialHex(patchPts, axis, origin, rInner, rOuter, axMin,
                                    axMax, pairViz.crossDirY, pairViz.crossDirZ);

  SweepBlockRegion block;
  block.patchId = patchId;
  block.pairedPatchId = pairedPatchId;
  block.primeId = primes[static_cast<size_t>(patchId)].id;
  block.pairedPrimeId =
      pairedPatchId >= 0 ? primes[static_cast<size_t>(pairedPatchId)].id : -1;
  block.kind = SweepKind::CylindricalBase;
  block.minIsoValue = -1.0f;
  block.maxIsoValue = 0.0f;
  block.axialLower = axMin;
  block.axialUpper = axMax;
  block.radialInner = rInner;
  block.radialOuter = rOuter;
  block.sweepAxis = axis;
  block.sweepOrigin = origin;
  block.crossDirY = pairViz.crossDirY;
  block.crossDirZ = pairViz.crossDirZ;
  block.isValid = true;

  candidateBlocks.push_back(std::move(block));
  blockHexVertices.push_back(std::move(hex));
}

void SweepDecomposer::BuildCylinderPatchBlockFromMesh(int patchId) {
  int primeId = primes[static_cast<size_t>(patchId)].id;
  Eigen::Vector3f axis = ExtractCylinderAxis(patchId).normalized();
  Eigen::Vector3f cen = ComputePatchCentroidFromMesh(patchId);
  if (cen.norm() < 1e-8f) {
    cen = ComputePatchCentroid(patchId);
  }
  Eigen::Vector3f origin = FootOnAxis(cen, cen, axis);

  std::vector<Eigen::Vector3f> pts;
  if (mesh) {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
      if (v->FeaturePoint() || v->label() != primeId) {
        continue;
      }
      pts.emplace_back(v->point()[0], v->point()[1], v->point()[2]);
    }
  }
  if (pts.size() < 2) {
    return;
  }

  CylinderPairViz viz;
  viz.patchIdA = patchId;
  viz.patchIdB = -1;
  viz.primeIdA = primeId;
  viz.primeIdB = -1;
  viz.sweepAxis = axis;
  viz.sweepOrigin = origin;
  BuildCylinderCrossFrame(axis, viz.crossDirY, viz.crossDirZ);

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float rMax = 0.0f;
  float oMinY = std::numeric_limits<float>::max();
  float oMaxY = std::numeric_limits<float>::lowest();
  float oMinZ = std::numeric_limits<float>::max();
  float oMaxZ = std::numeric_limits<float>::lowest();
  for (const auto &p : pts) {
    float ax = (p - origin).dot(axis);
    Eigen::Vector3f rad = RadialComponent(p, origin, axis);
    float r = rad.norm();
    axMin = std::min(axMin, ax);
    axMax = std::max(axMax, ax);
    rMax = std::max(rMax, r);
    oMinY = std::min(oMinY, rad.dot(viz.crossDirY));
    oMaxY = std::max(oMaxY, rad.dot(viz.crossDirY));
    oMinZ = std::min(oMinZ, rad.dot(viz.crossDirZ));
    oMaxZ = std::max(oMaxZ, rad.dot(viz.crossDirZ));
  }
  const float margin = std::max(0.02f, 0.5f * stepSize);
  viz.axialMid = 0.5f * (axMin + axMax);
  viz.radialOuter = rMax + margin;
  viz.radialInner = std::max(margin, 0.08f * viz.radialOuter);
  viz.outerCrossMinY = oMinY;
  viz.outerCrossMaxY = oMaxY;
  viz.outerCrossMinZ = oMinZ;
  viz.outerCrossMaxZ = oMaxZ;
  const float scale =
      viz.radialOuter > 1e-6f ? viz.radialInner / viz.radialOuter : 0.2f;
  viz.innerCrossMinY = oMinY * scale;
  viz.innerCrossMaxY = oMaxY * scale;
  viz.innerCrossMinZ = oMinZ * scale;
  viz.innerCrossMaxZ = oMaxZ * scale;

  BuildCylinderSinglePatchBlockFromMesh(patchId, -1, axis, origin, pts, viz);
}

void SweepDecomposer::BuildCylinderPairBlockFromMesh(int patchA, int patchB) {
  const PrimeData &primeA = primes[static_cast<size_t>(patchA)];
  const PrimeData &primeB = primes[static_cast<size_t>(patchB)];
  int idA = primeA.id;
  int idB = primeB.id;

  Eigen::Vector3f cenA = ComputePatchCentroidFromMesh(patchA);
  Eigen::Vector3f cenB = ComputePatchCentroidFromMesh(patchB);
  if (cenA.norm() < 1e-8f && cenB.norm() < 1e-8f) {
    cenA = ComputePatchCentroid(patchA);
    cenB = ComputePatchCentroid(patchB);
  }

  Eigen::Vector3f axis = ExtractCylinderAxis(patchA).normalized();
  Eigen::Vector3f mid = 0.5f * (cenA + cenB);
  Eigen::Vector3f origin = FootOnAxis(mid, cenA, axis);

  std::vector<Eigen::Vector3f> ptsA;
  std::vector<Eigen::Vector3f> ptsB;
  if (mesh) {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
      if (v->FeaturePoint()) {
        continue;
      }
      Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
      int lbl = v->label();
      if (lbl == idA) {
        ptsA.push_back(p);
      } else if (lbl == idB) {
        ptsB.push_back(p);
      }
    }
  }
  if (ptsA.size() < 2 || ptsB.size() < 2) {
    return;
  }

  CylinderPairViz viz =
      BuildCylinderPairViz(patchA, patchB, axis, origin, ptsA, ptsB);
  cylinderPairViz.push_back(viz);

  BuildCylinderSinglePatchBlockFromMesh(patchA, patchB, axis, origin, ptsA, viz);
  BuildCylinderSinglePatchBlockFromMesh(patchB, patchA, axis, origin, ptsB, viz);
}

bool SweepDecomposer::GradientsMatch(int patchA, int patchB) const {
  if (!IsCylinderPatch(patchA) || !IsCylinderPatch(patchB)) {
    return false;
  }

  Eigen::Vector3f axisA = ExtractCylinderAxis(patchA);
  Eigen::Vector3f axisB = ExtractCylinderAxis(patchB);
  if (std::abs(axisA.dot(axisB)) < 0.85f) {
    return false;
  }

  Eigen::Vector3f cenA = ComputePatchCentroidFromMesh(patchA);
  Eigen::Vector3f cenB = ComputePatchCentroidFromMesh(patchB);
  if (cenA.squaredNorm() + cenB.squaredNorm() < 1e-12f) {
    cenA = ComputePatchCentroid(patchA);
    cenB = ComputePatchCentroid(patchB);
  }

  Eigen::Vector3f axis = axisA.normalized();
  float axialSep = std::abs((cenB - cenA).dot(axis));
  if (axialSep <= std::max(0.15f, 2.0f * stepSize)) {
    return false;
  }

  int countA = 0;
  int countB = 0;
  Eigen::Vector3f gradA;
  Eigen::Vector3f gradB;
  if (mesh) {
    gradA = ComputePatchAvgAnalyticGradientFromMesh(patchA, &countA);
    gradB = ComputePatchAvgAnalyticGradientFromMesh(patchB, &countB);
  }
  if (!mesh || countA < 2 || countB < 2) {
    gradA = ComputePatchAvgGradient(patchA, &countA);
    gradB = ComputePatchAvgGradient(patchB, &countB);
  }
  if (countA < 2 || countB < 2) {
    return false;
  }

  float na = gradA.norm();
  float nb = gradB.norm();
  if (na < 1e-6f || nb < 1e-6f) {
    return false;
  }

  float gradAlign = std::abs(gradA.dot(gradB) / (na * nb));
  return gradAlign >= 0.55f;
}

std::vector<std::pair<int, int>>
SweepDecomposer::FindMatchingCylinderPairs() const {
  std::vector<int> cylinders;
  for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
    if (IsCylinderPatch(i)) {
      cylinders.push_back(i);
    }
  }

  std::vector<std::pair<int, int>> pairs;
  std::set<int> used;

  for (size_t i = 0; i < cylinders.size(); ++i) {
    if (used.count(cylinders[i])) {
      continue;
    }
    int bestJ = -1;
    float bestAlign = 0.0f;
    for (size_t j = i + 1; j < cylinders.size(); ++j) {
      if (used.count(cylinders[j])) {
        continue;
      }
      if (!GradientsMatch(cylinders[i], cylinders[j])) {
        continue;
      }
      Eigen::Vector3f axis = ExtractCylinderAxis(cylinders[i]).normalized();
      Eigen::Vector3f cenA = ComputePatchCentroidFromMesh(cylinders[i]);
      Eigen::Vector3f cenB = ComputePatchCentroidFromMesh(cylinders[j]);
      float axialSep = std::abs((cenB - cenA).dot(axis));
      int ca = 0, cb = 0;
      Eigen::Vector3f ga = ComputePatchAvgAnalyticGradientFromMesh(cylinders[i], &ca);
      Eigen::Vector3f gb = ComputePatchAvgAnalyticGradientFromMesh(cylinders[j], &cb);
      if (ca < 2 || cb < 2) {
        ga = ComputePatchAvgGradient(cylinders[i], &ca);
        gb = ComputePatchAvgGradient(cylinders[j], &cb);
      }
      float align = 0.0f;
      if (ga.norm() > 1e-8f && gb.norm() > 1e-8f) {
        align = std::abs(ga.normalized().dot(gb.normalized()));
      }
      float score = axialSep * (0.5f + align);
      if (score > bestAlign) {
        bestAlign = score;
        bestJ = static_cast<int>(j);
      }
    }
    if (bestJ >= 0) {
      pairs.emplace_back(cylinders[i], cylinders[static_cast<size_t>(bestJ)]);
      used.insert(cylinders[i]);
      used.insert(cylinders[static_cast<size_t>(bestJ)]);
    }
  }

  if (cylinders.size() >= 2 && pairs.empty()) {
    std::cout << "[SweepDecomposer] cylinder patches without gradient-matched pair:";
    for (int cid : cylinders) {
      int sc = 0;
      Eigen::Vector3f g = mesh ? ComputePatchAvgAnalyticGradientFromMesh(cid, &sc)
                               : ComputePatchAvgGradient(cid, &sc);
      std::cout << " " << cid << "(n=" << sc << ",g=" << g.transpose() << ")";
    }
    std::cout << "\n";
  }

  return pairs;
}

SweepKind SweepDecomposer::ClassifyPatch(int patchId) const {
  const PrimeData &p = primes[static_cast<size_t>(patchId)];
  if (p.isPlane) {
    return SweepKind::Translational;
  }

  const auto &q = p.params;
  double qxx = q[7], qyy = q[8], qzz = q[9];
  double cross = std::abs(q[4]) + std::abs(q[5]) + std::abs(q[6]);

  double spread =
      std::max({std::abs(qxx - qyy), std::abs(qyy - qzz), std::abs(qxx - qzz)});
  if (qxx > 1e-8 && qyy > 1e-8 && qzz > 1e-8 &&
      spread < 0.25 * (qxx + qyy + qzz)) {
    return SweepKind::Radial;
  }

  if (IsCylinderPatch(patchId)) {
    return SweepKind::CylindricalBase;
  }

  if (cross > 1e-6 ||
      spread < 0.35 * (std::abs(qxx) + std::abs(qyy) + std::abs(qzz))) {
    return SweepKind::Rotational;
  }

  return SweepKind::Translational;
}

Eigen::Vector3f SweepDecomposer::ComputePatchFrame(int patchId, SweepKind kind,
                                                   Eigen::Vector3f &outOrigin) const {
  const PrimeData &p = primes[patchId];
  const auto &params = p.params;
  outOrigin = Eigen::Vector3f::Zero();
  int count = 0;

  int primeId = p.id;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (FieldLabel[x][y][z] == primeId &&
            std::abs(Field[x][y][z]) < 2.5f * stepSize) {
          outOrigin += Coord[x][y][z];
          count++;
        }
      }
    }
  }
  if (count > 0) {
    outOrigin /= static_cast<float>(count);
  }

  if (kind == SweepKind::CylindricalBase) {
    return ExtractCylinderAxis(patchId);
  }

  if (kind == SweepKind::Translational || p.isPlane) {
    Eigen::Vector3f n(static_cast<float>(params[1]),
                      static_cast<float>(params[2]),
                      static_cast<float>(params[3]));
    if (n.norm() > 1e-8f) {
      return n.normalized();
    }
    return Eigen::Vector3f::UnitZ();
  }

  if (kind == SweepKind::Radial) {
    Eigen::Vector3f avgN = Eigen::Vector3f::Zero();
    int nCount = 0;
    for (int x = 0; x < D1; ++x) {
      for (int y = 0; y < D2; ++y) {
        for (int z = 0; z < D3; ++z) {
          if (FieldLabel[x][y][z] == primeId &&
              std::abs(Field[x][y][z]) < 2.5f * stepSize) {
            Eigen::Vector3f n = ComputeIsoSurfaceNormal(x, y, z);
            if (n.norm() > 1e-6f) {
              avgN += n;
              nCount++;
            }
          }
        }
      }
    }
    if (nCount > 0 && avgN.norm() > 1e-8f) {
      return avgN.normalized();
    }
    return Eigen::Vector3f::UnitZ();
  }

  // Rotational: axis ≈ 二次型最小特征方向；用梯度叉积估计
  Eigen::Vector3f axis = Eigen::Vector3f::Zero();
  int axisCount = 0;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (FieldLabel[x][y][z] == primeId &&
            std::abs(Field[x][y][z]) < 2.5f * stepSize) {
          Eigen::Vector3f g = ComputeIsoSurfaceNormal(x, y, z);
          Eigen::Vector3f r = Coord[x][y][z] - outOrigin;
          if (g.norm() > 1e-6f && r.norm() > stepSize) {
            axis += g.cross(r);
            axisCount++;
          }
        }
      }
    }
  }
  if (axisCount > 0 && axis.norm() > 1e-8f) {
    return axis.normalized();
  }

  Eigen::Vector3f fallback(static_cast<float>(params[1]),
                             static_cast<float>(params[2]),
                             static_cast<float>(params[3]));
  if (fallback.norm() > 1e-8f) {
    return fallback.normalized();
  }
  return Eigen::Vector3f::UnitZ();
}

Eigen::Vector3f SweepDecomposer::LocalSweepDirection(
    SweepKind kind, int x, int y, int z, const Eigen::Vector3f &axis,
    const Eigen::Vector3f &origin) const {
  const Eigen::Vector3f &p = Coord[x][y][z];

  if (kind == SweepKind::Translational) {
    return axis;
  }

  Eigen::Vector3f r = p - origin;
  if (kind == SweepKind::CylindricalBase || kind == SweepKind::Radial) {
    Eigen::Vector3f radial = r - r.dot(axis) * axis;
    if (radial.norm() > 1e-8f) {
      return radial.normalized();
    }
    if (kind == SweepKind::Radial && r.norm() > 1e-8f) {
      return r.normalized();
    }
    return axis;
  }

  // Rotational: 切向 = (r × axis)，再取与轴垂直的分量
  Eigen::Vector3f tangent = r.cross(axis);
  if (tangent.norm() > 1e-8f) {
    return tangent.normalized();
  }
  Eigen::Vector3f perp = r - r.dot(axis) * axis;
  if (perp.norm() > 1e-8f) {
    return perp.normalized();
  }
  return axis;
}

Eigen::Vector3f SweepDecomposer::ComputeIsoSurfaceNormal(int x, int y,
                                                         int z) const {
  float dfdx = 0.0f, dfdy = 0.0f, dfdz = 0.0f;

  if (x > 0 && x < D1 - 1) {
    dfdx = (Field[x + 1][y][z] - Field[x - 1][y][z]) / (2.0f * stepSize);
  } else if (x == 0 && D1 > 1) {
    dfdx = (Field[x + 1][y][z] - Field[x][y][z]) / stepSize;
  } else if (x == D1 - 1 && D1 > 1) {
    dfdx = (Field[x][y][z] - Field[x - 1][y][z]) / stepSize;
  }

  if (y > 0 && y < D2 - 1) {
    dfdy = (Field[x][y + 1][z] - Field[x][y - 1][z]) / (2.0f * stepSize);
  } else if (y == 0 && D2 > 1) {
    dfdy = (Field[x][y + 1][z] - Field[x][y][z]) / stepSize;
  } else if (y == D2 - 1 && D2 > 1) {
    dfdy = (Field[x][y][z] - Field[x][y - 1][z]) / stepSize;
  }

  if (z > 0 && z < D3 - 1) {
    dfdz = (Field[x][y][z + 1] - Field[x][y][z - 1]) / (2.0f * stepSize);
  } else if (z == 0 && D3 > 1) {
    dfdz = (Field[x][y][z + 1] - Field[x][y][z]) / stepSize;
  } else if (z == D3 - 1 && D3 > 1) {
    dfdz = (Field[x][y][z] - Field[x][y][z - 1]) / stepSize;
  }

  Eigen::Vector3f grad(dfdx, dfdy, dfdz);
  float norm = grad.norm();
  if (norm < 1e-8f) {
    return Eigen::Vector3f::Zero();
  }
  return grad / norm;
}

bool SweepDecomposer::CheckDirectionConstraint(
    int x, int y, int z, const Eigen::Vector3f &isoNormal) const {

  const Eigen::Vector3f &dir = GradField[x][y][z];
  float dirNorm = dir.norm();
  if (dirNorm < 1e-8f) {
    return false;
  }

  if (isoNormal.norm() < 1e-8f) {
    return true;
  }

  Eigen::Vector3f nd = dir / dirNorm;
  float cosAngle = std::abs(nd.dot(isoNormal));
  cosAngle = std::min(1.0f, std::max(0.0f, cosAngle));
  float angle = std::acos(cosAngle);

  float devPerp =
      std::abs(angle - static_cast<float>(M_PI) / 2.0f);
  float devTan = angle;

  return (devPerp < threshold) || (devTan < threshold);
}

float SweepDecomposer::SweepParameter(SweepKind kind, const Eigen::Vector3f &p,
                                      const Eigen::Vector3f &axis,
                                      const Eigen::Vector3f &origin) const {
  Eigen::Vector3f r = p - origin;

  if (kind == SweepKind::Translational) {
    return (p - origin).dot(axis);
  }

  if (kind == SweepKind::CylindricalBase || kind == SweepKind::Radial) {
    Eigen::Vector3f radial = r - r.dot(axis) * axis;
    return radial.norm();
  }

  // Rotational: 绕 axis 的极角
  Eigen::Vector3f radial = r - r.dot(axis) * axis;
  if (radial.norm() < 1e-8f) {
    return 0.0f;
  }
  Eigen::Vector3f ref =
      (std::abs(axis.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  Eigen::Vector3f u = (ref - ref.dot(axis) * axis).normalized();
  Eigen::Vector3f v = axis.cross(u);
  Eigen::Vector3f rn = radial.normalized();
  return std::atan2(rn.dot(v), rn.dot(u));
}

std::vector<std::vector<int>> SweepDecomposer::BuildCylinderSweepGroups() const {
  std::vector<int> cylinders;
  for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
    if (IsCylinderPatch(i)) {
      cylinders.push_back(i);
    }
  }
  if (cylinders.empty()) {
    return {};
  }

  std::vector<std::vector<int>> groups;
  std::set<int> used;
  for (int seed : cylinders) {
    if (used.count(seed)) {
      continue;
    }
    Eigen::Vector3f axisS = ExtractCylinderAxis(seed).normalized();
    Eigen::Vector3f cenS = ComputePatchCentroidFromMesh(seed);
    if (cenS.norm() < 1e-8f) {
      cenS = ComputePatchCentroid(seed);
    }

    std::vector<int> group;
    for (int other : cylinders) {
      if (used.count(other)) {
        continue;
      }
      Eigen::Vector3f axisO = ExtractCylinderAxis(other).normalized();
      if (std::abs(axisS.dot(axisO)) < 0.9f) {
        continue;
      }
      Eigen::Vector3f cenO = ComputePatchCentroidFromMesh(other);
      if (cenO.norm() < 1e-8f) {
        cenO = ComputePatchCentroid(other);
      }
      Eigen::Vector3f delta = cenO - cenS;
      Eigen::Vector3f perp = delta - delta.dot(axisS) * axisS;
      if (perp.norm() > 8.0f * stepSize) {
        continue;
      }
      group.push_back(other);
      used.insert(other);
    }
    if (!group.empty()) {
      groups.push_back(std::move(group));
    }
  }
  return groups;
}

float SweepDecomposer::CylinderSweepEnergy(
    int x, int y, int z, const Eigen::Vector3f &sweepDir) const {
  if (sweepDir.norm() < 1e-8f) {
    return 0.0f;
  }
  const Eigen::Vector3f &g = GradField[x][y][z];
  if (g.norm() < 1e-8f) {
    return 0.0f;
  }
  float c = std::abs(g.normalized().dot(sweepDir.normalized()));
  c = std::min(1.0f, std::max(0.0f, c));
  float ang = std::acos(c);
  float devPerp = std::abs(ang - static_cast<float>(M_PI) / 2.0f);
  if (devPerp < threshold || ang < threshold) {
    return 1.0f;
  }
  return 0.0f;
}

Eigen::Vector3f SweepDecomposer::LocalCylinderSweepDirection(
    int x, int y, int z, const Eigen::Vector3f &axis,
    const Eigen::Vector3f &origin, int bottomPrimeId, float axBottom,
    const Eigen::Vector3f &bottomSweepDir) const {
  const Eigen::Vector3f &p = Coord[x][y][z];
  int fl = FieldLabel[x][y][z];
  float ax = (p - origin).dot(axis);
  const float capBand = 2.5f * stepSize;

  if (bottomPrimeId >= 0 &&
      (fl == bottomPrimeId || std::abs(ax - axBottom) < capBand)) {
    if (bottomSweepDir.norm() > 1e-8f) {
      return bottomSweepDir.normalized();
    }
    return axis.normalized();
  }

  Eigen::Vector3f radial = RadialComponent(p, origin, axis);
  if (radial.norm() > 1e-8f) {
    return radial.normalized();
  }
  return axis.normalized();
}

void SweepDecomposer::GrowCylinderSweepRegion(
    const std::vector<int> &groupPatchIds) {
  if (groupPatchIds.empty()) {
    return;
  }

  int refPatch = groupPatchIds.front();
  Eigen::Vector3f axis = ExtractCylinderAxis(refPatch).normalized();
  Eigen::Vector3f refCen = ComputePatchCentroidFromMesh(refPatch);
  if (refCen.norm() < 1e-8f) {
    refCen = ComputePatchCentroid(refPatch);
  }
  Eigen::Vector3f origin = FootOnAxis(refCen, refCen, axis);

  std::vector<int> memberPrimeIds;
  memberPrimeIds.reserve(groupPatchIds.size());
  std::set<int> memberPrimeSet;
  for (int pid : groupPatchIds) {
    int id = primes[static_cast<size_t>(pid)].id;
    memberPrimeIds.push_back(id);
    memberPrimeSet.insert(id);
  }

  int bottomPatchId = groupPatchIds.front();
  float axBottom = std::numeric_limits<float>::max();
  for (int pid : groupPatchIds) {
    Eigen::Vector3f cen = ComputePatchCentroidFromMesh(pid);
    if (cen.norm() < 1e-8f) {
      cen = ComputePatchCentroid(pid);
    }
    float ax = (cen - origin).dot(axis);
    if (ax < axBottom) {
      axBottom = ax;
      bottomPatchId = pid;
    }
  }
  int bottomPrimeId = primes[static_cast<size_t>(bottomPatchId)].id;

  int bgCount = 0;
  Eigen::Vector3f bottomSweepDir =
      ComputePatchAvgAnalyticGradientFromMesh(bottomPatchId, &bgCount);
  if (bgCount < 2) {
    bottomSweepDir = ComputePatchAvgGradient(bottomPatchId, &bgCount);
  }
  if (bgCount < 2 || bottomSweepDir.norm() < 1e-8f) {
    bottomSweepDir = axis;
  } else {
    bottomSweepDir.normalize();
  }

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float rOuter = 0.0f;
  const float surfaceBand = 2.0f * stepSize;
  std::vector<Eigen::Vector3f> allPatchPts;

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (!memberPrimeSet.count(FieldLabel[x][y][z])) {
          continue;
        }
        if (std::abs(Field[x][y][z]) >= surfaceBand) {
          continue;
        }
        const Eigen::Vector3f &p = Coord[x][y][z];
        float ax = (p - origin).dot(axis);
        axMin = std::min(axMin, ax);
        axMax = std::max(axMax, ax);
        rOuter = std::max(rOuter, RadialComponent(p, origin, axis).norm());
        allPatchPts.push_back(p);
      }
    }
  }

  if (mesh) {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
      if (v->FeaturePoint() || !memberPrimeSet.count(v->label())) {
        continue;
      }
      Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
      allPatchPts.push_back(p);
      float ax = (p - origin).dot(axis);
      axMin = std::min(axMin, ax);
      axMax = std::max(axMax, ax);
      rOuter = std::max(rOuter, RadialComponent(p, origin, axis).norm());
    }
  }

  if (rOuter < stepSize || allPatchPts.size() < 4) {
    return;
  }

  axMin -= 2.0f * stepSize;
  axMax += 2.0f * stepSize;
  rOuter += 2.0f * stepSize;
  const float rInner = std::max(stepSize, 0.06f * rOuter);

  const float isoBand = 0.5f * stepSize;
  std::unordered_set<VoxelIndex, VoxelIndexHash> visited;
  std::unordered_set<VoxelIndex, VoxelIndexHash> blockVoxels;
  std::priority_queue<IsoVoxel, std::vector<IsoVoxel>, IsoVoxelGreater> pq;

  auto inAxialRange = [&](int x, int y, int z) {
    float ax = (Coord[x][y][z] - origin).dot(axis);
    return ax >= axMin && ax <= axMax;
  };
  auto radialDist = [&](int x, int y, int z) {
    return RadialComponent(Coord[x][y][z], origin, axis).norm();
  };

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] >= 0.0f || !inAxialRange(x, y, z)) {
          continue;
        }
        if (radialDist(x, y, z) > rInner) {
          continue;
        }
        VoxelIndex vi{x, y, z};
        visited.insert(vi);
        pq.push({Field[x][y][z], vi});
      }
    }
  }

  if (pq.empty()) {
    return;
  }

  float minIso = 0.0f;
  float maxIso = -isoBand;

  while (!pq.empty()) {
    IsoVoxel cur = pq.top();
    pq.pop();
    if (blockVoxels.count(cur.idx)) {
      continue;
    }

    int x = cur.idx.x, y = cur.idx.y, z = cur.idx.z;
    if (Field[x][y][z] >= 0.0f || !inAxialRange(x, y, z)) {
      continue;
    }
    if (radialDist(x, y, z) > rOuter + stepSize) {
      continue;
    }

    Eigen::Vector3f sweepDir = LocalCylinderSweepDirection(
        x, y, z, axis, origin, bottomPrimeId, axBottom, bottomSweepDir);
    if (CylinderSweepEnergy(x, y, z, sweepDir) < 0.5f) {
      continue;
    }

    Eigen::Vector3f isoN = ComputeIsoSurfaceNormal(x, y, z);
    if (!CheckDirectionConstraint(x, y, z, isoN)) {
      continue;
    }

    blockVoxels.insert(cur.idx);
    minIso = std::min(minIso, Field[x][y][z]);
    maxIso = std::max(maxIso, Field[x][y][z]);

    for (int d = 0; d < 6; ++d) {
      int nx = x + kDx[d], ny = y + kDy[d], nz = z + kDz[d];
      if (nx < 0 || nx >= D1 || ny < 0 || ny >= D2 || nz < 0 || nz >= D3) {
        continue;
      }
      if (Field[nx][ny][nz] >= 0.0f || !inAxialRange(nx, ny, nz)) {
        continue;
      }
      if (radialDist(nx, ny, nz) > rOuter + stepSize) {
        continue;
      }
      if (Field[nx][ny][nz] < Field[x][y][z] - isoBand) {
        continue;
      }
      VoxelIndex nvi{nx, ny, nz};
      if (visited.count(nvi)) {
        continue;
      }
      visited.insert(nvi);
      pq.push({Field[nx][ny][nz], nvi});
    }
  }

  if (blockVoxels.empty()) {
    return;
  }

  Eigen::Vector3f crossY, crossZ;
  BuildCylinderCrossFrame(axis, crossY, crossZ);

  int innerPatch = groupPatchIds.front();
  int outerPatch = groupPatchIds.front();
  float minR = std::numeric_limits<float>::max();
  float maxR = 0.0f;
  for (int pid : groupPatchIds) {
    if (pid == bottomPatchId) {
      continue;
    }
    float avgR = 0.0f;
    int cnt = 0;
    int primeId = primes[static_cast<size_t>(pid)].id;
    if (mesh) {
      for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
        auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
        if (v->label() != primeId) {
          continue;
        }
        avgR += RadialComponent(
                   Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]),
                   origin, axis)
                   .norm();
        cnt++;
      }
    }
    if (cnt > 0) {
      avgR /= static_cast<float>(cnt);
      if (avgR < minR) {
        minR = avgR;
        innerPatch = pid;
      }
      if (avgR > maxR) {
        maxR = avgR;
        outerPatch = pid;
      }
    }
  }

  std::map<int, Eigen::Vector3f> hex;

  std::vector<Eigen::Vector3f> ptsA;
  std::vector<Eigen::Vector3f> ptsB;
  int idInner = primes[static_cast<size_t>(innerPatch)].id;
  int idOuter = primes[static_cast<size_t>(outerPatch)].id;
  if (mesh) {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
      if (v->FeaturePoint()) {
        continue;
      }
      Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
      if (v->label() == idInner) {
        ptsA.push_back(p);
      } else if (v->label() == idOuter) {
        ptsB.push_back(p);
      }
    }
  }
  if (ptsA.size() >= 2 && ptsB.size() >= 2) {
    cylinderPairViz.push_back(
        BuildCylinderPairViz(innerPatch, outerPatch, axis, origin, ptsA, ptsB));
  }

  SweepBlockRegion block;
  block.patchId = innerPatch;
  block.pairedPatchId = outerPatch;
  block.bottomPatchId = bottomPatchId;
  block.primeId = idInner;
  block.pairedPrimeId = idOuter;
  block.memberPrimeIds = memberPrimeIds;
  block.kind = SweepKind::CylindricalBase;
  block.minIsoValue = minIso;
  block.maxIsoValue = maxIso;
  block.axialLower = axMin;
  block.axialUpper = axMax;
  block.radialInner = rInner;
  block.radialOuter = rOuter;
  block.sweepAxis = axis;
  block.sweepOrigin = origin;
  block.crossDirY = crossY;
  block.crossDirZ = crossZ;
  block.isValid = true;

  for (const auto &v : blockVoxels) {
    if (allInteriorVoxels.empty() || allInteriorVoxels.count(v)) {
      block.coveredVoxels.insert(v);
    }
  }

  if (!block.coveredVoxels.empty()) {
    float grownAxMin = std::numeric_limits<float>::max();
    float grownAxMax = std::numeric_limits<float>::lowest();
    float grownRMin = std::numeric_limits<float>::max();
    float grownRMax = 0.0f;
    std::vector<Eigen::Vector3f> grownPts;
    grownPts.reserve(block.coveredVoxels.size());
    for (const auto &v : block.coveredVoxels) {
      const Eigen::Vector3f &p = Coord[v.x][v.y][v.z];
      grownPts.push_back(p);
      float ax = (p - origin).dot(axis);
      float rad = RadialComponent(p, origin, axis).norm();
      grownAxMin = std::min(grownAxMin, ax);
      grownAxMax = std::max(grownAxMax, ax);
      grownRMin = std::min(grownRMin, rad);
      grownRMax = std::max(grownRMax, rad);
    }
    const float grownMargin = std::max(stepSize, 0.05f);
    block.axialLower = grownAxMin - grownMargin;
    block.axialUpper = grownAxMax + grownMargin;
    block.radialInner = std::max(stepSize, grownRMin - grownMargin);
    block.radialOuter = grownRMax + grownMargin;
    hex = BuildCylinderRadialHex(grownPts, axis, origin, block.radialInner,
                                 block.radialOuter, block.axialLower,
                                 block.axialUpper, crossY, crossZ);

    const size_t grownVoxels = block.coveredVoxels.size();
    candidateBlocks.push_back(std::move(block));
    blockHexVertices.push_back(std::move(hex));
    std::cout << "[SweepDecomposer]   grown voxels=" << grownVoxels
              << " inner=" << innerPatch << " outer=" << outerPatch
              << " bottom=" << bottomPatchId << std::endl;
  }
}

void SweepDecomposer::GrowBlockBetweenCylinderPair(int patchA, int patchB) {
  const PrimeData &primeA = primes[static_cast<size_t>(patchA)];
  const PrimeData &primeB = primes[static_cast<size_t>(patchB)];
  int idA = primeA.id;
  int idB = primeB.id;

  Eigen::Vector3f cenA = ComputePatchCentroid(patchA);
  Eigen::Vector3f cenB = ComputePatchCentroid(patchB);
  Eigen::Vector3f axis = ExtractCylinderAxis(patchA).normalized();
  Eigen::Vector3f origin = FootOnAxis(0.5f * (cenA + cenB), cenA, axis);

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float rOuter = 0.0f;
  const float surfaceBand = 2.0f * stepSize;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        int fl = FieldLabel[x][y][z];
        if (fl != idA && fl != idB) {
          continue;
        }
        if (std::abs(Field[x][y][z]) >= surfaceBand) {
          continue;
        }
        const Eigen::Vector3f &p = Coord[x][y][z];
        float ax = (p - origin).dot(axis);
        axMin = std::min(axMin, ax);
        axMax = std::max(axMax, ax);
        rOuter = std::max(rOuter, RadialComponent(p, origin, axis).norm());
      }
    }
  }
  if (rOuter < stepSize) {
    return;
  }

  axMin -= 2.0f * stepSize;
  axMax += 2.0f * stepSize;
  rOuter += 2.0f * stepSize;
  const float rInner = std::max(stepSize, 0.08f * rOuter);

  const float isoBand = 0.5f * stepSize;
  std::unordered_set<VoxelIndex, VoxelIndexHash> visited;
  std::unordered_set<VoxelIndex, VoxelIndexHash> blockVoxels;
  std::priority_queue<IsoVoxel, std::vector<IsoVoxel>, IsoVoxelGreater> pq;

  auto inAxialRange = [&](int x, int y, int z) {
    float ax = (Coord[x][y][z] - origin).dot(axis);
    return ax >= axMin && ax <= axMax;
  };

  auto radialDist = [&](int x, int y, int z) {
    return RadialComponent(Coord[x][y][z], origin, axis).norm();
  };

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] >= 0.0f || !inAxialRange(x, y, z)) {
          continue;
        }
        if (radialDist(x, y, z) > rInner) {
          continue;
        }
        VoxelIndex vi{x, y, z};
        visited.insert(vi);
        pq.push({Field[x][y][z], vi});
      }
    }
  }

  if (pq.empty()) {
    return;
  }

  float minIso = 0.0f;
  float maxIso = -isoBand;

  while (!pq.empty()) {
    IsoVoxel cur = pq.top();
    pq.pop();

    if (blockVoxels.count(cur.idx)) {
      continue;
    }

    int x = cur.idx.x, y = cur.idx.y, z = cur.idx.z;
    if (Field[x][y][z] >= 0.0f || !inAxialRange(x, y, z)) {
      continue;
    }
    if (radialDist(x, y, z) > rOuter + stepSize) {
      continue;
    }

    Eigen::Vector3f isoN = ComputeIsoSurfaceNormal(x, y, z);
    if (!CheckDirectionConstraint(x, y, z, isoN)) {
      continue;
    }

    Eigen::Vector3f localSweep =
        LocalSweepDirection(SweepKind::CylindricalBase, x, y, z, axis, origin);
    if (localSweep.norm() > 1e-8f && GradField[x][y][z].norm() > 1e-8f) {
      Eigen::Vector3f fd = GradField[x][y][z].normalized();
      float c = std::abs(fd.dot(localSweep.normalized()));
      c = std::min(1.0f, std::max(0.0f, c));
      float ang = std::acos(c);
      float devPerp = std::abs(ang - static_cast<float>(M_PI) / 2.0f);
      if (devPerp >= threshold && ang >= threshold) {
        continue;
      }
    }

    blockVoxels.insert(cur.idx);
    minIso = std::min(minIso, Field[x][y][z]);
    maxIso = std::max(maxIso, Field[x][y][z]);

    for (int d = 0; d < 6; ++d) {
      int nx = x + kDx[d], ny = y + kDy[d], nz = z + kDz[d];
      if (nx < 0 || nx >= D1 || ny < 0 || ny >= D2 || nz < 0 || nz >= D3) {
        continue;
      }
      if (Field[nx][ny][nz] >= 0.0f || !inAxialRange(nx, ny, nz)) {
        continue;
      }
      if (radialDist(nx, ny, nz) > rOuter + stepSize) {
        continue;
      }
      // 从柱轴附近向外生长：仅向更接近表面的体素扩展
      if (Field[nx][ny][nz] < Field[x][y][z] - isoBand) {
        continue;
      }

      VoxelIndex nvi{nx, ny, nz};
      if (visited.count(nvi)) {
        continue;
      }
      visited.insert(nvi);
      pq.push({Field[nx][ny][nz], nvi});
    }
  }

  if (blockVoxels.empty()) {
    return;
  }

  Eigen::Vector3f crossY, crossZ;
  BuildCylinderCrossFrame(axis, crossY, crossZ);

  SweepBlockRegion block;
  block.patchId = patchA;
  block.pairedPatchId = patchB;
  block.primeId = idA;
  block.pairedPrimeId = idB;
  block.kind = SweepKind::CylindricalBase;
  block.minIsoValue = minIso;
  block.maxIsoValue = maxIso;
  block.axialLower = axMin;
  block.axialUpper = axMax;
  block.radialInner = rInner;
  block.radialOuter = rOuter;
  block.sweepAxis = axis;
  block.sweepOrigin = origin;
  block.crossDirY = crossY;
  block.crossDirZ = crossZ;
  block.isValid = true;

  for (const auto &v : blockVoxels) {
    if (allInteriorVoxels.count(v)) {
      block.coveredVoxels.insert(v);
    }
  }

  if (!block.coveredVoxels.empty() && ValidateSweepTopology(block)) {
    candidateBlocks.push_back(std::move(block));
  }
}

void SweepDecomposer::GrowBlockFromPatch(int patchId) {
  const PrimeData &prime = primes[patchId];
  int primeId = prime.id;
  SweepKind kind = ClassifyPatch(patchId);
  Eigen::Vector3f origin;
  Eigen::Vector3f axis = ComputePatchFrame(patchId, kind, origin);

  const float surfaceBand = 2.0f * stepSize;
  const float isoBand = 0.5f * stepSize;

  std::unordered_set<VoxelIndex, VoxelIndexHash> visited;
  std::unordered_set<VoxelIndex, VoxelIndexHash> blockVoxels;

  // 优先队列：从表面（Field 接近 0 的内部体素）向内生长
  std::priority_queue<IsoVoxel, std::vector<IsoVoxel>, IsoVoxelGreater> pq;

  auto isOnPatchSurface = [&](int x, int y, int z) {
    return FieldLabel[x][y][z] == primeId &&
           std::abs(Field[x][y][z]) < surfaceBand;
  };

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] >= 0.0f) {
          continue;
        }
        bool adjacentSurface = false;
        for (int d = 0; d < 6; ++d) {
          int nx = x + kDx[d], ny = y + kDy[d], nz = z + kDz[d];
          if (nx < 0 || nx >= D1 || ny < 0 || ny >= D2 || nz < 0 || nz >= D3) {
            continue;
          }
          if (isOnPatchSurface(nx, ny, nz)) {
            adjacentSurface = true;
            break;
          }
        }
        if (adjacentSurface) {
          VoxelIndex vi{x, y, z};
          visited.insert(vi);
          pq.push({Field[x][y][z], vi});
        }
      }
    }
  }

  if (pq.empty()) {
    return;
  }

  float minIso = 0.0f;
  float maxIso = -isoBand;

  while (!pq.empty()) {
    IsoVoxel cur = pq.top();
    pq.pop();

    if (blockVoxels.count(cur.idx)) {
      continue;
    }

    int x = cur.idx.x, y = cur.idx.y, z = cur.idx.z;
    if (Field[x][y][z] >= 0.0f) {
      continue;
    }

    Eigen::Vector3f isoN = ComputeIsoSurfaceNormal(x, y, z);
    if (!CheckDirectionConstraint(x, y, z, isoN)) {
      continue;
    }

    // 广义扫掠：场方向还应与局部扫掠方向一致（允许垂直/相切于扫掠轨迹）
    Eigen::Vector3f localSweep =
        LocalSweepDirection(kind, x, y, z, axis, origin);
    if (localSweep.norm() > 1e-8f && GradField[x][y][z].norm() > 1e-8f) {
      Eigen::Vector3f fd = GradField[x][y][z].normalized();
      float c = std::abs(fd.dot(localSweep.normalized()));
      c = std::min(1.0f, std::max(0.0f, c));
      float ang = std::acos(c);
      float devPerp = std::abs(ang - static_cast<float>(M_PI) / 2.0f);
      if (devPerp >= threshold && ang >= threshold) {
        continue;
      }
    }

    blockVoxels.insert(cur.idx);
    minIso = std::min(minIso, Field[x][y][z]);
    maxIso = std::max(maxIso, Field[x][y][z]);

    for (int d = 0; d < 6; ++d) {
      int nx = x + kDx[d], ny = y + kDy[d], nz = z + kDz[d];
      if (nx < 0 || nx >= D1 || ny < 0 || ny >= D2 || nz < 0 || nz >= D3) {
        continue;
      }
      if (Field[nx][ny][nz] >= 0.0f) {
        continue;
      }
      // 仅向更深处（更负的等值面）扩展
      if (Field[nx][ny][nz] > Field[x][y][z] + isoBand) {
        continue;
      }

      VoxelIndex nvi{nx, ny, nz};
      if (visited.count(nvi)) {
        continue;
      }
      visited.insert(nvi);
      pq.push({Field[nx][ny][nz], nvi});
    }
  }

  if (blockVoxels.empty()) {
    return;
  }

  SweepBlockRegion block;
  block.patchId = patchId;
  block.primeId = primeId;
  block.kind = kind;
  block.minIsoValue = minIso;
  block.maxIsoValue = maxIso;
  block.sweepAxis = axis;
  block.sweepOrigin = origin;
  block.isValid = true;

  for (const auto &v : blockVoxels) {
    if (allInteriorVoxels.count(v)) {
      block.coveredVoxels.insert(v);
    }
  }

  if (!block.coveredVoxels.empty() && ValidateSweepTopology(block)) {
    candidateBlocks.push_back(std::move(block));
  }
}

bool SweepDecomposer::ValidateSweepTopology(
    const SweepBlockRegion &block) const {
  if (block.coveredVoxels.size() < 4) {
    return false;
  }

  float sMin = std::numeric_limits<float>::max();
  float sMax = std::numeric_limits<float>::lowest();

  for (const auto &v : block.coveredVoxels) {
    float s = SweepParameter(block.kind, Coord[v.x][v.y][v.z], block.sweepAxis,
                             block.sweepOrigin);
    sMin = std::min(sMin, s);
    sMax = std::max(sMax, s);
  }

  float extent = sMax - sMin;
  if (extent < 1.5f * stepSize) {
    return false;
  }

  const float sliceTol = stepSize * 1.2f;
  int outerBase = 0;
  int innerBase = 0;
  int sideMid = 0;
  float sMid = 0.5f * (sMin + sMax);

  for (const auto &v : block.coveredVoxels) {
    float s = SweepParameter(block.kind, Coord[v.x][v.y][v.z], block.sweepAxis,
                             block.sweepOrigin);
    if (std::abs(s - sMin) < sliceTol) {
      outerBase++;
    }
    if (std::abs(s - sMax) < sliceTol) {
      innerBase++;
    }
    if (std::abs(s - sMid) < sliceTol) {
      sideMid++;
    }
  }

  if (block.kind == SweepKind::CylindricalBase && block.pairedPatchId >= 0) {
    int idA = primes[static_cast<size_t>(block.patchId)].id;
    int idB = primes[static_cast<size_t>(block.pairedPatchId)].id;
    int innerHits = 0;
    int outerA = 0;
    int outerB = 0;

    for (const auto &v : block.coveredVoxels) {
      float s = SweepParameter(block.kind, Coord[v.x][v.y][v.z], block.sweepAxis,
                               block.sweepOrigin);
      int fl = FieldLabel[v.x][v.y][v.z];
      if (std::abs(s - sMin) < sliceTol) {
        innerHits++;
      }
      if (std::abs(s - sMax) < sliceTol) {
        if (fl == idA) {
          outerA++;
        }
        if (fl == idB) {
          outerB++;
        }
      }
    }

    bool hasInnerCore = innerHits >= 1;
    bool hasOuterPair = outerA >= 1 && outerB >= 1;
    return hasInnerCore && hasOuterPair && sideMid >= 1;
  }

  return outerBase >= 2 && innerBase >= 2 && sideMid >= 1;
}

void SweepDecomposer::GreedySetCover() {
  selectedBlocks.clear();

  if (allInteriorVoxels.empty() || candidateBlocks.empty()) {
    return;
  }

  std::set<VoxelIndex> uncovered = allInteriorVoxels;
  std::vector<bool> used(candidateBlocks.size(), false);

  while (!uncovered.empty()) {
    int bestIdx = -1;
    size_t bestGain = 0;

    for (int i = 0; i < static_cast<int>(candidateBlocks.size()); ++i) {
      if (used[i]) {
        continue;
      }
      size_t gain = 0;
      for (const auto &v : candidateBlocks[i].coveredVoxels) {
        if (uncovered.count(v)) {
          gain++;
        }
      }
      if (gain > bestGain) {
        bestGain = gain;
        bestIdx = i;
      }
    }

    if (bestIdx < 0 || bestGain == 0) {
      std::cout << "[SweepDecomposer] warning: " << uncovered.size()
                << " interior voxels uncovered\n";
      break;
    }

    used[bestIdx] = true;
    selectedBlocks.push_back(candidateBlocks[bestIdx]);
    for (const auto &v : candidateBlocks[bestIdx].coveredVoxels) {
      uncovered.erase(v);
    }
  }
}

std::vector<std::map<int, Eigen::Vector3f>>
SweepDecomposer::GetBlockHexVertices() const {
  if (!blockHexVertices.empty() &&
      blockHexVertices.size() == selectedBlocks.size()) {
    return blockHexVertices;
  }

  std::vector<std::map<int, Eigen::Vector3f>> result;

  for (const auto &block : selectedBlocks) {
    if (block.coveredVoxels.empty()) {
      continue;
    }

    Eigen::Vector3f dirX = block.sweepAxis;
    if (block.kind == SweepKind::Radial) {
      dirX = LocalSweepDirection(block.kind, block.coveredVoxels.begin()->x,
                                 block.coveredVoxels.begin()->y,
                                 block.coveredVoxels.begin()->z,
                                 block.sweepAxis, block.sweepOrigin);
    }

    Eigen::Vector3f arbitrary =
        (std::abs(dirX.dot(Eigen::Vector3f::UnitX())) < 0.9f)
            ? Eigen::Vector3f::UnitX()
            : Eigen::Vector3f::UnitY();
    Eigen::Vector3f dirY =
        (arbitrary - dirX.dot(arbitrary) * dirX).normalized();
    Eigen::Vector3f dirZ = dirX.cross(dirY).normalized();

    float minX = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float minY = std::numeric_limits<float>::max();
    float maxY = std::numeric_limits<float>::lowest();
    float minZ = std::numeric_limits<float>::max();
    float maxZ = std::numeric_limits<float>::lowest();

    for (const auto &v : block.coveredVoxels) {
      const Eigen::Vector3f &p = Coord[v.x][v.y][v.z];
      float px = p.dot(dirX), py = p.dot(dirY), pz = p.dot(dirZ);
      minX = std::min(minX, px);
      maxX = std::max(maxX, px);
      minY = std::min(minY, py);
      maxY = std::max(maxY, py);
      minZ = std::min(minZ, pz);
      maxZ = std::max(maxZ, pz);
    }

    Eigen::Matrix3f A;
    A.row(0) = dirX.transpose();
    A.row(1) = dirY.transpose();
    A.row(2) = dirZ.transpose();

    std::map<int, Eigen::Vector3f> vertices;
    float bx[] = {minX, maxX};
    float by[] = {minY, maxY};
    float bz[] = {minZ, maxZ};

    int idx = 0;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          Eigen::Vector3f rhs(bx[i], by[j], bz[k]);
          vertices[idx++] = A.colPivHouseholderQr().solve(rhs);
        }
      }
    }

    result.push_back(vertices);
  }

  return result;
}
