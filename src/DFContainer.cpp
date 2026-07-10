#include "DFContainer.h"
#include "CTMesh.h"
#include "ColorImplementer.h"
#include "SweepBlock.h"
#include "SweepDirDetector.h"
#include "SweepDirFilter.h"
#include "SweepDirSpliter.h"
#include "SweepFaceImprinter.h"
#include "SweepHexMesher.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <unordered_map>

#ifdef ENABLE_CUDA
#include "DistanceFieldCUDA.cuh"
#endif

#ifdef ENABLE_METAL
#include "MetalNearestPoint.h"
#endif

DistanceField::DistanceField() { this->primes.clear(); };

DistanceField::DistanceField(MeshLib::CTMesh *mesh) {
  this->primes.clear();
  SetMesh(mesh);
}

void DistanceField::exportPlanesToFile(const std::string &filename) {
  std::ofstream outFile(filename);
  if (!outFile.is_open())
    return;

  std::set<std::string> seenPlanes;

  outFile << "# Blender Cutting Planes Data\n";
  outFile << "planes_params = [\n";

  const int faceIndices[6][4] = {{0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
                                 {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5}};

  for (const auto &hex : CuttingHexLists) {
    for (int f_idx = 0; f_idx < 6; ++f_idx) {
      const Eigen::Vector3f &pA = hex.at(faceIndices[f_idx][0]);
      const Eigen::Vector3f &pB = hex.at(faceIndices[f_idx][1]);
      const Eigen::Vector3f &pC = hex.at(faceIndices[f_idx][2]);

      Eigen::Vector3f nQuad = (pB - pA).cross(pC - pA);
      if (nQuad.norm() < 1e-6f)
        continue;
      nQuad.normalize();

      float D = -nQuad.dot(pA);

      char buf[128];
      sprintf(buf, "%.4f,%.4f,%.4f,%.4f", nQuad.x(), nQuad.y(), nQuad.z(), D);
      std::string planeId(buf);

      sprintf(buf, "%.4f,%.4f,%.4f,%.4f", -nQuad.x(), -nQuad.y(), -nQuad.z(),
              -D);
      std::string invPlaneId(buf);

      if (seenPlanes.find(planeId) == seenPlanes.end() &&
          seenPlanes.find(invPlaneId) == seenPlanes.end()) {
        outFile << "    (" << std::fixed << std::setprecision(6) << nQuad.x()
                << ", " << nQuad.y() << ", " << nQuad.z() << ", " << D
                << "),\n";
        seenPlanes.insert(planeId);
      }
    }
  }

  outFile << "]\n";
  outFile.close();
}

void DistanceField::SetMesh(MeshLib::CTMesh *mesh) {
  this->mesh = mesh;
  this->PointList.clear();
  this->PointIDList.clear();
  this->VertexPtrList.clear();
  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); mfiter++) {
    MeshLib::CToolFace *face =
        static_cast<MeshLib::CToolFace *>(mfiter.value());
    auto normal = (face->halfedge()->target()->point() -
                   face->halfedge()->source()->point()) ^
                  (face->halfedge()->he_next()->target()->point() -
                   face->halfedge()->he_next()->source()->point());
    face->area() = 0.5 * abs(normal.norm());
    normal /= normal.norm();
    face->normal() = normal / face->area();
  }
  for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
    MeshLib::CToolVertex *tv =
        static_cast<MeshLib::CToolVertex *>(viter.value());
    tv->FeaturePoint() = false;
    std::vector<float> pointcoord;
    pointcoord.push_back(viter.value()->point()[0]);
    pointcoord.push_back(viter.value()->point()[1]);
    pointcoord.push_back(viter.value()->point()[2]);
    this->PointList.push_back(pointcoord);
    this->PointIDList.push_back(viter.value()->id());
    this->VertexPtrList.push_back(tv);
    Eigen::Vector3f normal = Eigen::Vector3f(0.0, 0.0, 0.0);
    MeshLib::CTMesh::CVertex *v = *viter;
    int label = -2;
    for (MeshLib::CTMesh::VertexFaceIterator vfiter(v); !vfiter.end();
         vfiter++) {
      MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(vfiter.value());
      int facelabel = f->label();
      if (f->area() < 1e-10)
        continue;
      normal += Eigen::Vector3f(f->normal()[0] * f->area(),
                                f->normal()[1] * f->area(),
                                f->normal()[2] * f->area());
      if (label == -2)
        label = facelabel;
      else if (label != facelabel)
        tv->FeaturePoint() = true;
    }
    normal.normalize();
    v->normal()[0] = normal[0];
    v->normal()[1] = normal[1];
    v->normal()[2] = normal[2];
  }
}

void DistanceField::GridScalar(int MinScatter) {
  float MinX = FLT_MAX;
  float MaxX = FLT_MIN;
  float MaxY = FLT_MIN;
  float MinY = FLT_MAX;
  float MaxZ = FLT_MIN;
  float MinZ = FLT_MAX;

  for (int i = 0; i < PointList.size(); i++) {
    if (MinX > PointList[i][0]) {
      MinX = PointList[i][0];
    }
    if (MaxX < PointList[i][0]) {
      MaxX = PointList[i][0];
    }
    if (MinY > PointList[i][1]) {
      MinY = PointList[i][1];
    }
    if (MaxY < PointList[i][1]) {
      MaxY = PointList[i][1];
    }
    if (MinZ > PointList[i][2]) {
      MinZ = PointList[i][2];
    }
    if (MaxZ < PointList[i][2]) {
      MaxZ = PointList[i][2];
    }
  }

  float GapX = MaxX - MinX;
  float GapY = MaxY - MinY;
  float GapZ = MaxZ - MinZ;
  MinX -= 0.3 * GapX;
  MaxX += 0.3 * GapX;
  MinY -= 0.3 * GapY;
  MaxY += 0.3 * GapY;
  MinZ -= 0.3 * GapZ;
  MaxZ += 0.3 * GapZ;

  float MinPatch;
  if (MaxZ - MinZ > MaxY - MinY)
    MinPatch = MaxY - MinY;
  else
    MinPatch = MaxZ - MinZ;
  if (MaxX - MinX < MinPatch)
    MinPatch = MaxX - MinX;
  this->PatchSize = MinPatch / MinScatter;
  for (int i = 0; i < (MaxX - MinX) / PatchSize; i++) {
    std::vector<std::vector<float>> FieldX;
    std::vector<std::vector<Eigen::Vector3f>> CoordX;
    for (int j = 0; j < (MaxY - MinY) / PatchSize; j++) {
      std::vector<float> FieldXY;
      std::vector<Eigen::Vector3f> CoordXY;
      for (int k = 0; k < (MaxZ - MinZ) / PatchSize; k++) {
        float FieldXYZ = MinX + i * PatchSize;
        Eigen::Vector3f CoordXYZ;
        CoordXYZ[0] = MinX + i * PatchSize;
        CoordXYZ[1] = MinY + j * PatchSize;
        CoordXYZ[2] = MinZ + k * PatchSize;
        FieldXY.push_back(FieldXYZ);
        CoordXY.push_back(CoordXYZ);
      }
      FieldX.push_back(FieldXY);
      CoordX.push_back(CoordXY);
    }

    this->Field.push_back(FieldX);
    this->Coord.push_back(CoordX);
  }
}

void DistanceField::ExtractSweepDir() {
  if (this->mesh == NULL)
    return;
  this->SweepDir.clear();

  SweepDirDetector detector(this->mesh, this->primes);
  this->SweepDir = detector.GetSweepDir();

  std::cout << "Logging Potential Sweep Direction..." << std::endl;

  for (int i = 0; i < this->SweepDir.size(); i++) {
    std::cout << "Potential SweepDir: " << this->SweepDir[i] << std::endl;
  }

  return;
}

void DistanceField::BuildOctree() {
  if (PointList.empty())
    return;

  Eigen::Vector3f minPoint =
      Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
  Eigen::Vector3f maxPoint =
      Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());

  for (const auto &point : PointList) {
    Eigen::Vector3f p(point[0], point[1], point[2]);
    minPoint = minPoint.cwiseMin(p);
    maxPoint = maxPoint.cwiseMax(p);
  }

  Eigen::Vector3f center = (minPoint + maxPoint) * 0.5f;
  float halfSize = (maxPoint - minPoint).norm() * 0.5f + 0.1f;

  octreeRoot = std::make_shared<OctreeNode>(center, halfSize);

  std::vector<int> allIndices(PointList.size());
  for (int i = 0; i < PointList.size(); ++i) {
    allIndices[i] = i;
  }

  BuildOctreeRecursive(octreeRoot, allIndices, 0);
}

void DistanceField::BuildOctreeRecursive(std::shared_ptr<OctreeNode> node,
                                         const std::vector<int> &pointIndices,
                                         int depth) {
  if (pointIndices.empty())
    return;

  if (pointIndices.size() <= maxPointsPerNode || depth >= maxDepth) {
    node->pointIndices = pointIndices;
    node->isLeaf = true;
    return;
  }

  SubdivideNode(node);
  node->isLeaf = false;

  std::vector<std::vector<int>> childIndices(8);
  for (int idx : pointIndices) {
    const auto &point = PointList[idx];
    Eigen::Vector3f p(point[0], point[1], point[2]);

    int childIndex = 0;
    if (p.x() > node->center.x())
      childIndex |= 1;
    if (p.y() > node->center.y())
      childIndex |= 2;
    if (p.z() > node->center.z())
      childIndex |= 4;

    childIndices[childIndex].push_back(idx);
  }

  for (int i = 0; i < 8; ++i) {
    if (!childIndices[i].empty()) {
      BuildOctreeRecursive(node->children[i], childIndices[i], depth + 1);
    }
  }
}

void DistanceField::SubdivideNode(std::shared_ptr<OctreeNode> node) {
  float childHalfSize = node->halfSize * 0.5f;

  for (int i = 0; i < 8; ++i) {
    Eigen::Vector3f childCenter = node->center;
    childCenter.x() += (i & 1) ? childHalfSize : -childHalfSize;
    childCenter.y() += (i & 2) ? childHalfSize : -childHalfSize;
    childCenter.z() += (i & 4) ? childHalfSize : -childHalfSize;

    node->children.push_back(
        std::make_shared<OctreeNode>(childCenter, childHalfSize));
  }
}

int DistanceField::FindNearestPointInOctree(
    const Eigen::Vector3f &point, std::shared_ptr<OctreeNode> node,
    float &bestDist) {
  if (!node)
    return -1;

  Eigen::Vector3f diff = (point - node->center).cwiseAbs();
  float distToNode =
      (diff - Eigen::Vector3f::Constant(node->halfSize)).cwiseMax(0.0f).norm();

  if (distToNode > bestDist)
    return -1;

  int bestIdx = -1;

  if (node->isLeaf) {
    for (int idx : node->pointIndices) {
      const auto &p = PointList[idx];
      float d = (point - Eigen::Vector3f(p[0], p[1], p[2])).norm();
      if (d < bestDist) {
        bestDist = d;
        bestIdx = idx;
      }
    }
    return bestIdx;
  }

  for (const auto &child : node->children) {
    int found = FindNearestPointInOctree(point, child, bestDist);
    if (found >= 0)
      bestIdx = found;
  }
  return bestIdx;
}

Eigen::Vector3f DistanceField::ClosestPointOnTriangle(
    const Eigen::Vector3f &point, const Eigen::Vector3f &v0,
    const Eigen::Vector3f &v1, const Eigen::Vector3f &v2) {
  Eigen::Vector3f edge0 = v1 - v0;
  Eigen::Vector3f edge1 = v2 - v0;
  Eigen::Vector3f v0ToPoint = point - v0;

  float a = edge0.dot(edge0);
  float b = edge0.dot(edge1);
  float c = edge1.dot(edge1);
  float d = edge0.dot(v0ToPoint);
  float e = edge1.dot(v0ToPoint);

  float det = a * c - b * b;
  float s = b * e - c * d;
  float t = b * d - a * e;

  if (s + t < det) {
    if (s < 0.0f) {
      if (t < 0.0f) {
        if (d < 0.0f) {
          s = std::clamp(-d / a, 0.0f, 1.0f);
          t = 0.0f;
        } else {
          s = 0.0f;
          t = std::clamp(-e / c, 0.0f, 1.0f);
        }
      } else {
        s = 0.0f;
        t = std::clamp(-e / c, 0.0f, 1.0f);
      }
    } else if (t < 0.0f) {
      s = std::clamp(-d / a, 0.0f, 1.0f);
      t = 0.0f;
    } else {
      float invDet = 1.0f / det;
      s *= invDet;
      t *= invDet;
    }
  } else {
    if (s < 0.0f) {
      float tmp0 = b + d;
      float tmp1 = c + e;
      if (tmp1 > tmp0) {
        float numer = tmp1 - tmp0;
        float denom = a - 2 * b + c;
        s = std::clamp(numer / denom, 0.0f, 1.0f);
        t = 1 - s;
      } else {
        t = std::clamp(-e / c, 0.0f, 1.0f);
        s = 0.0f;
      }
    } else if (t < 0.0f) {
      if (a + d > b + e) {
        float numer = c + e - b - d;
        float denom = a - 2 * b + c;
        s = std::clamp(numer / denom, 0.0f, 1.0f);
        t = 1 - s;
      } else {
        s = std::clamp(-e / c, 0.0f, 1.0f);
        t = 0.0f;
      }
    } else {
      float numer = c + e - b - d;
      float denom = a - 2 * b + c;
      s = std::clamp(numer / denom, 0.0f, 1.0f);
      t = 1.0f - s;
    }
  }

  return v0 + s * edge0 + t * edge1;
}

float DistanceField::PointToTriangleDistance(const Eigen::Vector3f &point,
                                             const Eigen::Vector3f &v0,
                                             const Eigen::Vector3f &v1,
                                             const Eigen::Vector3f &v2) {
  Eigen::Vector3f closestPoint = ClosestPointOnTriangle(point, v0, v1, v2);
  return (point - closestPoint).norm();
}

Eigen::Vector4f DistanceField::ComputeVertexDistance(
    const Eigen::Vector3f &point, MeshLib::CToolVertex *nearestVertex,
    int x, int y, int z) {
  Eigen::Vector4f DistanceVector;
  DistanceVector.setZero();
  if (!nearestVertex)
    return DistanceVector;

  float minDistance =
      (point - Eigen::Vector3f(nearestVertex->point()[0],
                               nearestVertex->point()[1],
                               nearestVertex->point()[2]))
          .norm();

  if (!this->primes.empty()) {
    bool isFeature = nearestVertex->FeaturePoint();
    if (isFeature) {
      Eigen::Vector3f pointOnMesh(nearestVertex->point()[0],
                                  nearestVertex->point()[1],
                                  nearestVertex->point()[2]);
      Eigen::Vector3f NormalOnMesh(nearestVertex->normal()[0],
                                   nearestVertex->normal()[1],
                                   nearestVertex->normal()[2]);
      NormalOnMesh.normalize();
      if ((pointOnMesh - point).dot(NormalOnMesh) < 0)
        minDistance = -abs(minDistance);
      else
        minDistance = abs(minDistance);
      DistanceVector[0] = 0;
      DistanceVector[1] = 0;
      DistanceVector[2] = 0;
      DistanceVector[3] = minDistance;
      this->FieldLabel[x][y][z] = -1;
    } else {
      minDistance = this->DisCompute(point, nearestVertex->label());
      this->FieldLabel[x][y][z] = nearestVertex->label();
      Eigen::Vector3f VertexNormal(nearestVertex->normal()[0],
                                   nearestVertex->normal()[1],
                                   nearestVertex->normal()[2]);
      Eigen::Vector3f vertexPoint(nearestVertex->point()[0],
                                  nearestVertex->point()[1],
                                  nearestVertex->point()[2]);
      if ((vertexPoint - point).dot(VertexNormal) < 0) {
        DistanceVector[0] = -nearestVertex->normal()[0];
        DistanceVector[1] = -nearestVertex->normal()[1];
        DistanceVector[2] = -nearestVertex->normal()[2];
        DistanceVector[3] = -minDistance;
      } else {
        DistanceVector[0] = nearestVertex->normal()[0];
        DistanceVector[1] = nearestVertex->normal()[1];
        DistanceVector[2] = nearestVertex->normal()[2];
        DistanceVector[3] = minDistance;
      }
    }
  } else {
    Eigen::Vector3f pointOnMesh(nearestVertex->point()[0],
                                nearestVertex->point()[1],
                                nearestVertex->point()[2]);
    Eigen::Vector3f NormalOnMesh(nearestVertex->normal()[0],
                                 nearestVertex->normal()[1],
                                 nearestVertex->normal()[2]);
    if ((pointOnMesh - point).dot(NormalOnMesh) < 0)
      minDistance = -minDistance;
    DistanceVector[0] = nearestVertex->normal()[0];
    DistanceVector[1] = nearestVertex->normal()[1];
    DistanceVector[2] = nearestVertex->normal()[2];
    DistanceVector[3] = minDistance;
  }
  return DistanceVector;
}

bool DistanceField::HasNonPlanarPrimes() const {
  for (const auto &p : this->primes) {
    if (p.params.size() >= 10 && !p.isPlane) {
      return true;
    }
  }
  return false;
}

bool DistanceField::PrimeLabelValid(int label) const {
  return label >= 0 && label < static_cast<int>(this->primes.size());
}

const PrimeData *DistanceField::GetPrimeByLabel(int label) const {
  if (!PrimeLabelValid(label)) {
    return nullptr;
  }
  return &this->primes[label];
}

void DistanceField::ReindexPrimesById() {
  if (this->primes.empty()) {
    return;
  }

  int maxId = 0;
  for (const auto &p : this->primes) {
    maxId = std::max(maxId, p.id);
  }

  std::vector<PrimeData> byId(static_cast<size_t>(maxId) + 1);
  for (const auto &p : this->primes) {
    if (p.id >= 0 && p.id <= maxId) {
      byId[static_cast<size_t>(p.id)] = p;
    }
  }
  this->primes = std::move(byId);
}

double DistanceField::DisCompute(Eigen::Vector3f point, int label) {
  const PrimeData *prime = GetPrimeByLabel(label);
  if (!prime) {
    return std::numeric_limits<double>::infinity();
  }
  auto &m_params = prime->params;

  if (prime->isPlane) {
    double a = m_params[1], b = m_params[2], c = m_params[3], d = m_params[0];
    double norm = sqrt(a * a + b * b + c * c);

    if (norm < 1e-10) {
      return std::numeric_limits<double>::infinity();
    }

    return std::abs(a * point[0] + b * point[1] + c * point[2] + d) / norm;
  }

  int max_iter = 20;
  double lambda = 0.0;
  Eigen::Vector3f q = point;

  for (int i = 0; i < max_iter; ++i) {
    double x = q(0), y = q(1), z = q(2);

    double F = m_params[0] + m_params[1] * x + m_params[2] * y +
               m_params[3] * z + m_params[4] * x * y + m_params[5] * x * z +
               m_params[6] * y * z + m_params[7] * x * x + m_params[8] * y * y +
               m_params[9] * z * z;

    Eigen::Vector3f gradF;
    gradF << m_params[1] + m_params[4] * y + m_params[5] * z +
                 2 * m_params[7] * x,
        m_params[2] + m_params[4] * x + m_params[6] * z + 2 * m_params[8] * y,
        m_params[3] + m_params[5] * x + m_params[6] * y + 2 * m_params[9] * z;

    Eigen::Vector4f residual;
    residual.head<3>() = q - point + lambda * gradF;
    residual(3) = F;

    Eigen::Matrix4f J;
    Eigen::Matrix3f hessian;
    hessian << 2 * m_params[7], m_params[4], m_params[5], m_params[4],
        2 * m_params[8], m_params[6], m_params[5], m_params[6], 2 * m_params[9];

    J.block<3, 3>(0, 0) = Eigen::Matrix3f::Identity() + lambda * hessian;
    J.block<3, 1>(0, 3) = gradF;
    J.block<1, 3>(3, 0) = gradF.transpose();
    J(3, 3) = 0.0;

    Eigen::Vector4f delta = J.colPivHouseholderQr().solve(-residual);

    if (!delta.allFinite()) {
      break;
    }

    q += delta.head<3>();
    lambda += delta(3);

    if (delta.norm() < 1e-6 && std::abs(F) < 1e-6) {
      break;
    }
  }

  return (q - point).norm();
};

void ComputeNearestPointsCPU(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<float>> &PointList,
    std::vector<std::vector<std::vector<int>>> &NearestIndex) {
  int xSize = (int)Coord.size();
  if (xSize == 0)
    return;
  int ySize = (int)Coord[0].size();
  int zSize = (int)Coord[0][0].size();
  int numPoints = (int)PointList.size();

  NearestIndex.resize(xSize);
  for (int i = 0; i < xSize; ++i) {
    NearestIndex[i].resize(ySize);
    for (int j = 0; j < ySize; ++j)
      NearestIndex[i][j].resize(zSize, -1);
  }

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        const Eigen::Vector3f &query = Coord[i][j][k];
        float bestDist = std::numeric_limits<float>::max();
        int bestIdx = -1;
        for (int p = 0; p < numPoints; ++p) {
          float dx = query.x() - PointList[p][0];
          float dy = query.y() - PointList[p][1];
          float dz = query.z() - PointList[p][2];
          float d2 = dx * dx + dy * dy + dz * dz;
          if (d2 < bestDist) {
            bestDist = d2;
            bestIdx = p;
          }
        }
        NearestIndex[i][j][k] = bestIdx;
      }
    }
  }
}

void DistanceField::ComputeDistanceField() {
  if (PointList.empty() || Field.empty() || Coord.empty()) {
    return;
  }

  int xSize = Field.size();
  int ySize = Field[0].size();
  int zSize = Field[0][0].size();
  GradianceCount.resize(xSize);
  GradianceField.resize(xSize);
  FieldLabel.resize(xSize);
  for (int i = 0; i < xSize; ++i) {
    GradianceCount[i].resize(ySize);
    GradianceField[i].resize(ySize);
    FieldLabel[i].resize(ySize);
    for (int j = 0; j < ySize; ++j) {
      GradianceCount[i][j].resize(zSize, 0);
      GradianceField[i][j].resize(zSize);
      FieldLabel[i][j].resize(zSize);
    }
  }

  std::vector<std::vector<std::vector<int>>> NearestPoint;

#ifdef ENABLE_CUDA
  std::cout << "Using CUDA for nearest point computation..." << std::endl;
  ComputeNearestPointsCUDA(this->Coord, this->PointList, NearestPoint);
#elif defined(ENABLE_METAL)
  std::cout << "Using Metal GPU for nearest point computation..." << std::endl;
  ComputeNearestPointsMetal(this->Coord, this->PointList, NearestPoint);
#else
  std::cout << "Using CPU for nearest point computation..." << std::endl;
  ComputeNearestPointsCPU(this->Coord, this->PointList, NearestPoint);
#endif

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        Eigen::Vector3f point = this->Coord[i][j][k];
        MeshLib::CToolVertex *nearestVertex =
            this->VertexPtrList[NearestPoint[i][j][k]];
        Eigen::Vector4f distance =
            ComputeVertexDistance(point, nearestVertex, i, j, k);
        GradianceField[i][j][k] = distance.head(3);
        Field[i][j][k] = distance[3];
      }
    }
  }

  GradianceCount.resize(xSize);
  GradianceField.resize(xSize);
  GradianceDiff.resize(xSize);
  for (int i = 0; i < xSize; ++i) {
    GradianceCount[i].resize(ySize);
    GradianceField[i].resize(ySize);
    GradianceDiff[i].resize(ySize);
    for (int j = 0; j < ySize; ++j) {
      GradianceCount[i][j].resize(zSize);
      GradianceField[i][j].resize(zSize);
      GradianceDiff[i][j].resize(zSize);
    }
  }
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        if (this->Field[i][j][k] < 0) {
          GradianceCount[i][j][k] = 0;
          GradianceDiff[i][j][k] = 0;
          continue;
        }
        float MaxGradianceDiff = 0;

        if (i > 0) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i - 1][j][k]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i - 1][j][k]));
        }

        if (i < xSize - 1) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i + 1][j][k]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i + 1][j][k]));
        }
        if (j > 0) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i][j - 1][k]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i][j - 1][k]));
        }
        if (j < ySize - 1) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i][j + 1][k]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i][j + 1][k]));
        }
        if (k > 0) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i][j][k - 1]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i][j][k - 1]));
        }
        if (k < zSize - 1) {
          MaxGradianceDiff =
              MaxGradianceDiff > acos(this->GradianceField[i][j][k].dot(
                                     this->GradianceField[i][j][k + 1]))
                  ? MaxGradianceDiff
                  : acos(this->GradianceField[i][j][k].dot(
                        this->GradianceField[i][j][k + 1]));
        }
        if (this->GradianceField[i][j][k].norm() < 1e-16)
          MaxGradianceDiff = 0;
        this->GradianceDiff[i][j][k] = MaxGradianceDiff;
        int currentLabel = this->FieldLabel[i][j][k];
        if (currentLabel < 0) {
          GradianceCount[i][j][k] = 0;
        }
        bool hasDifferent = false;

        if (i > 0 && this->FieldLabel[i - 1][j][k] != currentLabel &&
            this->FieldLabel[i - 1][j][k] != -1) {
          hasDifferent = true;
        } else if (i < xSize - 1 &&
                   this->FieldLabel[i + 1][j][k] != currentLabel &&
                   this->FieldLabel[i + 1][j][k] != -1) {
          hasDifferent = true;
        }

        if (!hasDifferent) {
          if (j > 0 && this->FieldLabel[i][j - 1][k] != currentLabel &&
              this->FieldLabel[i][j - 1][k] != -1) {
            hasDifferent = true;
          } else if (j < ySize - 1 &&
                     this->FieldLabel[i][j + 1][k] != currentLabel &&
                     this->FieldLabel[i][j + 1][k] != -1) {
            hasDifferent = true;
          }
        }

        if (!hasDifferent) {
          if (k > 0 && this->FieldLabel[i][j][k - 1] != -1 &&
              this->FieldLabel[i][j][k - 1] != currentLabel) {
            hasDifferent = true;
          } else if (k < zSize - 1 && this->FieldLabel[i][j][k + 1] != -1 &&
                     this->FieldLabel[i][j][k + 1] != currentLabel) {
            hasDifferent = true;
          }
        }

        GradianceCount[i][j][k] = hasDifferent ? 1 : 0;
      }
    }
  }

  InitForbiddenBoundaryPoints();
}

void DistanceField::InitForbiddenBoundaryPoints() {
  if (Field.empty() || FieldLabel.empty()) {
    return;
  }

  int D1 = static_cast<int>(Field.size());
  int D2 = static_cast<int>(Field[0].size());
  int D3 = static_cast<int>(Field[0][0].size());
  float patch = (Coord[0][0][0] - Coord[0][0][1]).norm();
  const float fieldThreshold = patch;

  ForbiddenBoundaryPoints.assign(
      D1, std::vector<std::vector<bool>>(D2, std::vector<bool>(D3, false)));

  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        int fl = FieldLabel[x][y][z];
        if (fl < 0 || fl >= static_cast<int>(primes.size())) {
          continue;
        }
        if (!primes[static_cast<size_t>(fl)].isPlane &&
            std::abs(Field[x][y][z]) < fieldThreshold) {
          ForbiddenBoundaryPoints[x][y][z] = true;
        }
      }
    }
  }
}

Eigen::Vector3f DistanceField::RandomSweepColor(int seed) {
  auto hue01 = [](int s, int salt) -> float {
    uint32_t x = static_cast<uint32_t>(s * 374761393 + salt * 668265263);
    x = (x ^ (x >> 13)) * 1274126177u;
    x ^= x >> 16;
    return static_cast<float>(x % 1000) / 1000.0f;
  };

  float h = hue01(seed, 17);
  float s = 0.55f + 0.35f * hue01(seed, 29);
  float v = 0.75f + 0.2f * hue01(seed, 41);
  int hi = static_cast<int>(h * 6.0f) % 6;
  float f = h * 6.0f - static_cast<float>(hi);
  float p = v * (1.0f - s);
  float q = v * (1.0f - f * s);
  float t = v * (1.0f - (1.0f - f) * s);
  switch (hi) {
  case 0:
    return Eigen::Vector3f(v, t, p);
  case 1:
    return Eigen::Vector3f(q, v, p);
  case 2:
    return Eigen::Vector3f(p, v, t);
  case 3:
    return Eigen::Vector3f(p, q, v);
  case 4:
    return Eigen::Vector3f(t, p, v);
  default:
    return Eigen::Vector3f(v, p, q);
  }
}

void DistanceField::EnsureSweepBlockColors() {
  while (sweepBlockColors.size() < CuttingHexLists.size()) {
    sweepBlockColors.push_back(
        RandomSweepColor(static_cast<int>(sweepBlockColors.size())));
  }
  while (sweepBlockNonPlanar.size() < CuttingHexLists.size()) {
    sweepBlockNonPlanar.push_back(false);
  }
}

namespace {

Eigen::Vector3f RadialFromAxis(const Eigen::Vector3f &p,
                               const Eigen::Vector3f &origin,
                               const Eigen::Vector3f &axis) {
  Eigen::Vector3f rel = p - origin;
  return rel - rel.dot(axis) * axis;
}

} // namespace

int DistanceField::HexIndexToSweepBlockIndex(int hexIdx) const {
  if (hexIdx < 0) {
    return -1;
  }
  int planarHexCount = static_cast<int>(CuttingHexLists.size()) -
                       static_cast<int>(sweepBlocks.size());
  if (hexIdx < planarHexCount) {
    return -1;
  }
  int sweepIdx = hexIdx - planarHexCount;
  if (sweepIdx >= static_cast<int>(sweepBlocks.size())) {
    return -1;
  }
  return sweepIdx;
}

bool DistanceField::IsPointInCylinderSweepBlock(
    int hexIdx, const Eigen::Vector3f &position) const {
  int sweepIdx = HexIndexToSweepBlockIndex(hexIdx);
  if (sweepIdx < 0) {
    return false;
  }
  const SweepBlockRegion &block = sweepBlocks[static_cast<size_t>(sweepIdx)];
  if (block.kind != SweepKind::CylindricalBase) {
    return false;
  }
  Eigen::Vector3f axis = block.sweepAxis.normalized();
  float ax = (position - block.sweepOrigin).dot(axis);
  if (ax < block.axialLower || ax > block.axialUpper) {
    return false;
  }
  float radial = RadialFromAxis(position, block.sweepOrigin, axis).norm();
  const float margin = std::max(STEP_SIZE * 0.5f, 1e-3f);
  return radial >= block.radialInner - margin &&
         radial <= block.radialOuter + margin;
}

int DistanceField::SweepBlockToHexIndex(int sweepBlockIdx) const {
  if (sweepBlockIdx < 0) {
    return -1;
  }
  int planarHexCount = static_cast<int>(CuttingHexLists.size()) -
                       static_cast<int>(sweepBlocks.size());
  if (planarHexCount < 0) {
    return -1;
  }
  int hexIdx = planarHexCount + sweepBlockIdx;
  if (hexIdx >= static_cast<int>(CuttingHexLists.size())) {
    return -1;
  }
  return hexIdx;
}

float DistanceField::CuttingHexVolume(
    const std::map<int, Eigen::Vector3f> &verticesMap) {
  if (verticesMap.size() != 8) {
    return std::numeric_limits<float>::max();
  }
  const Eigen::Vector3f &v0 = verticesMap.at(0);
  const Eigen::Vector3f &v1 = verticesMap.at(1);
  const Eigen::Vector3f &v2 = verticesMap.at(2);
  const Eigen::Vector3f &v4 = verticesMap.at(4);
  float vol = std::abs((v4 - v0).dot((v2 - v0).cross(v1 - v0)));
  return vol > 1e-12f ? vol : std::numeric_limits<float>::max();
}

int DistanceField::FindBasePlanarHexIndex() const {
  if (CuttingHexLists.empty()) {
    return -1;
  }

  Eigen::Vector3f baseCen = Eigen::Vector3f::Zero();
  bool foundBase = false;
  if (mesh) {
    for (const auto &prime : primes) {
      if (!prime.isPlane) {
        continue;
      }
      Eigen::Vector3f cen = Eigen::Vector3f::Zero();
      int count = 0;
      for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
        auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
        if (v->label() != prime.id) {
          continue;
        }
        cen += Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
        count++;
      }
      if (count == 0) {
        continue;
      }
      cen /= static_cast<float>(count);
      if (!foundBase || cen.z() < baseCen.z()) {
        baseCen = cen;
        foundBase = true;
      }
    }
  }

  if (foundBase) {
    for (int i = 0; i < static_cast<int>(CuttingHexLists.size()); ++i) {
      if (i < static_cast<int>(sweepBlockNonPlanar.size()) &&
          sweepBlockNonPlanar[static_cast<size_t>(i)]) {
        continue;
      }
      if (CuttingHexLists[static_cast<size_t>(i)].size() != 8) {
        continue;
      }
      if (insideCuttingBox(baseCen, CuttingHexLists[static_cast<size_t>(i)])) {
        return i;
      }
    }
  }

  for (int i = 0; i < static_cast<int>(CuttingHexLists.size()); ++i) {
    if (i < static_cast<int>(sweepBlockNonPlanar.size()) &&
        sweepBlockNonPlanar[static_cast<size_t>(i)]) {
      continue;
    }
    if (CuttingHexLists[static_cast<size_t>(i)].size() == 8) {
      return i;
    }
  }
  return -1;
}

std::vector<int> DistanceField::getDisplayHexIndices() const {
  std::vector<int> indices;
  int baseHex = FindBasePlanarHexIndex();
  if (baseHex >= 0) {
    indices.push_back(baseHex);
  }
  for (int i = 0; i < static_cast<int>(CuttingHexLists.size()); ++i) {
    if (i < static_cast<int>(sweepBlockNonPlanar.size()) &&
        sweepBlockNonPlanar[static_cast<size_t>(i)]) {
      indices.push_back(i);
    }
  }
  return indices;
}

int DistanceField::FindSweepBlockForPrimeLabel(int primeLabel) const {
  if (primeLabel < 0) {
    return -1;
  }
  for (int i = 0; i < static_cast<int>(sweepBlocks.size()); ++i) {
    const SweepBlockRegion &block = sweepBlocks[static_cast<size_t>(i)];
    for (int memberId : block.memberPrimeIds) {
      if (memberId == primeLabel) {
        return SweepBlockToHexIndex(i);
      }
    }
  }
  for (int i = 0; i < static_cast<int>(sweepBlocks.size()); ++i) {
    if (sweepBlocks[static_cast<size_t>(i)].primeId == primeLabel) {
      return SweepBlockToHexIndex(i);
    }
  }
  return -1;
}

int DistanceField::FindSweepBlockForFace(
    const Eigen::Vector3f &position,
    const std::unordered_map<int, int> &labelVotes) const {
  // 距离场体素归属优先（两扫掠体的分割完全由距离场决定）
  int voxelBlock = FindSweepBlockForPoint(position);
  if (voxelBlock >= 0) {
    return voxelBlock;
  }

  int blockIdx = -1;
  int bestVotes = 0;
  for (const auto &[lbl, votes] : labelVotes) {
    int candidate = FindSweepBlockForPrimeLabel(lbl);
    if (candidate >= 0 && votes > bestVotes) {
      blockIdx = candidate;
      bestVotes = votes;
    }
  }
  return blockIdx;
}

int DistanceField::FindSweepBlockForPoint(const Eigen::Vector3f &position) const {
  // 优先使用距离场体素归属（两扫掠体管线）
  if (!voxelToBlock.empty()) {
    int D1 = static_cast<int>(Coord.size());
    int D2 = D1 > 0 ? static_cast<int>(Coord[0].size()) : 0;
    int D3 = D2 > 0 ? static_cast<int>(Coord[0][0].size()) : 0;
    VoxelIndex c = WorldToVoxel(position);
    int best = -1;
    float bestDist = std::numeric_limits<float>::max();
    for (int r = 0; r <= 4; ++r) {
      for (int dx = -r; dx <= r; ++dx) {
        for (int dy = -r; dy <= r; ++dy) {
          for (int dz = -r; dz <= r; ++dz) {
            if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != r) {
              continue;
            }
            int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
            if (nx < 0 || nx >= D1 || ny < 0 || ny >= D2 || nz < 0 ||
                nz >= D3) {
              continue;
            }
            auto it = voxelToBlock.find({nx, ny, nz});
            if (it == voxelToBlock.end()) {
              continue;
            }
            float dist = (Coord[nx][ny][nz] - position).squaredNorm();
            if (dist < bestDist) {
              bestDist = dist;
              best = it->second;
            }
          }
        }
      }
      if (best >= 0) {
        return best;
      }
    }
    return -1;
  }

  // 回退：旧的平面盒包含测试（纯平面管线）
  for (int i = static_cast<int>(CuttingHexLists.size()) - 1; i >= 0; --i) {
    if (i < static_cast<int>(sweepBlockNonPlanar.size()) &&
        sweepBlockNonPlanar[static_cast<size_t>(i)] &&
        IsPointInCylinderSweepBlock(i, position)) {
      return i;
    }
  }

  int bestPlanar = -1;
  float bestVolume = std::numeric_limits<float>::max();
  for (int i = 0; i < static_cast<int>(CuttingHexLists.size()); ++i) {
    if (i < static_cast<int>(sweepBlockNonPlanar.size()) &&
        sweepBlockNonPlanar[static_cast<size_t>(i)]) {
      continue;
    }
    const auto &hex = CuttingHexLists[static_cast<size_t>(i)];
    if (hex.size() != 8) {
      continue;
    }
    if (!insideCuttingBox(position, hex)) {
      continue;
    }
    float vol = CuttingHexVolume(hex);
    if (vol < bestVolume) {
      bestVolume = vol;
      bestPlanar = i;
    }
  }
  return bestPlanar;
}

Eigen::Vector3f
DistanceField::SweepDirectionAt(int blockIdx,
                              const Eigen::Vector3f &position) const {
  if (blockIdx < 0) {
    return Eigen::Vector3f::Zero();
  }
  if (blockIdx < static_cast<int>(sweepBlocks.size())) {
    const SweepBlockRegion &block = sweepBlocks[static_cast<size_t>(blockIdx)];
    if (block.kind == SweepKind::CylindricalBase) {
      Eigen::Vector3f r = position - block.sweepOrigin;
      Eigen::Vector3f radial = r - r.dot(block.sweepAxis) * block.sweepAxis;
      if (radial.norm() > 1e-6f) {
        return radial.normalized();
      }
    }
  }
  if (blockIdx < static_cast<int>(SweepDir.size())) {
    const Eigen::Vector3f &d = SweepDir[static_cast<size_t>(blockIdx)];
    if (d.norm() > 1e-6f) {
      return d.normalized();
    }
  }
  return Eigen::Vector3f::UnitZ();
}

Eigen::Vector3f
DistanceField::EncodeSweepDirColor(const Eigen::Vector3f &dir) {
  if (dir.norm() < 1e-8f) {
    return Eigen::Vector3f(0.5f, 0.5f, 0.5f);
  }
  Eigen::Vector3f n = dir.normalized();
  return 0.5f * (n + Eigen::Vector3f::Ones());
}

void DistanceField::AppendSweepBlocks(
    const std::vector<SweepBlockRegion> &blocks,
    const std::vector<std::map<int, Eigen::Vector3f>> &hexes, bool nonPlanar) {
  size_t n = std::min(blocks.size(), hexes.size());
  for (size_t i = 0; i < n; ++i) {
    const auto &block = blocks[i];
    this->sweepBlocks.push_back(block);
    this->SweepDir.push_back(block.sweepAxis);
    bool blockNonPlanar =
        nonPlanar || block.kind == SweepKind::CylindricalBase;
    this->sweepBlockNonPlanar.push_back(blockNonPlanar);
    this->sweepBlockColors.push_back(
        RandomSweepColor(static_cast<int>(this->sweepBlockColors.size())));
    this->CuttingHexLists.push_back(hexes[i]);
  }
}

void DistanceField::ApplySweepVisualization() {
  if (!this->mesh || this->CuttingHexLists.empty()) {
    return;
  }

  EnsureSweepBlockColors();

  auto blockColor = [&](int idx) -> Eigen::Vector3f {
    if (idx >= 0 && idx < static_cast<int>(sweepBlockColors.size())) {
      return sweepBlockColors[static_cast<size_t>(idx)];
    }
    return RandomSweepColor(idx);
  };

  auto applyBlockColor = [&](MeshLib::CToolVertex *v, int blockIdx) {
    v->cflabel() = blockIdx;
    if (blockIdx >= 0) {
      Eigen::Vector3f c = blockColor(blockIdx);
      v->rgb()[0] = c.x();
      v->rgb()[1] = c.y();
      v->rgb()[2] = c.z();
    }
  };

  auto applyFaceBlockColor = [&](MeshLib::CToolFace *f, int blockIdx) {
    f->sweeplabel() = blockIdx;
    if (blockIdx >= 0) {
      Eigen::Vector3f c = blockColor(blockIdx);
      f->rgb()[0] = c.x();
      f->rgb()[1] = c.y();
      f->rgb()[2] = c.z();
      if (blockIdx < static_cast<int>(sweepBlockNonPlanar.size()) &&
          sweepBlockNonPlanar[static_cast<size_t>(blockIdx)]) {
        f->sweepFaceType() = 4;
      }
    }
  };

  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolVertex *v =
        static_cast<MeshLib::CToolVertex *>(mviter.value());
    Eigen::Vector3f position =
        Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
    // 体素归属优先（按距离场空间位置，而非 prime 标签），避免穿过两体的
    // 同一柱面 prime 被整片强行归到管体。
    int blockIdx = FindSweepBlockForPoint(position);
    if (blockIdx < 0) {
      blockIdx = FindSweepBlockForPrimeLabel(v->label());
    }
    applyBlockColor(v, blockIdx);
  }

  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); ++mfiter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mfiter.value());
    int count = 0;
    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    std::unordered_map<int, int> labelVotes;
    for (MeshLib::CTMesh::FaceVertexIterator fviter(f); !fviter.end();
         ++fviter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(fviter.value());
      count++;
      position += Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
      labelVotes[v->label()]++;
    }
    position /= std::max(count, 1);

    applyFaceBlockColor(f, FindSweepBlockForFace(position, labelVotes));
  }

  std::cout << "Starting Implementer..." << std::endl;
  Implementer implementer(this->mesh);
  std::cout << "Implementer done, assigning sweep block colors..." << std::endl;

  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); ++mfiter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mfiter.value());
    std::unordered_map<int, int> labelVotes;
    for (MeshLib::CTMesh::FaceVertexIterator fviter(f); !fviter.end();
         ++fviter) {
      labelVotes[fviter.value()->label()]++;
    }
    int count = 0;
    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    for (MeshLib::CTMesh::FaceVertexIterator fviter(f); !fviter.end();
         ++fviter) {
      count++;
      position += Eigen::Vector3f(fviter.value()->point()[0],
                                  fviter.value()->point()[1],
                                  fviter.value()->point()[2]);
    }
    position /= std::max(count, 1);
    int blockIdx = FindSweepBlockForFace(position, labelVotes);
    if (blockIdx >= 0) {
      applyFaceBlockColor(f, blockIdx);
    } else if (f->sweeplabel() >= 0) {
      applyFaceBlockColor(f, f->sweeplabel());
    }
  }

  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolVertex *v =
        static_cast<MeshLib::CToolVertex *>(mviter.value());
    Eigen::Vector3f position(v->point()[0], v->point()[1], v->point()[2]);
    // 空间体素归属优先，保证穿过两体的同一柱面 prime 按位置分别着色
    int spatialBlock = FindSweepBlockForPoint(position);
    if (spatialBlock >= 0) {
      applyBlockColor(v, spatialBlock);
      continue;
    }
    std::unordered_map<int, int> votes;
    for (MeshLib::CTMesh::VertexFaceIterator vfiter(v); !vfiter.end();
         ++vfiter) {
      auto *f = static_cast<MeshLib::CToolFace *>(vfiter.value());
      if (f->sweeplabel() >= 0) {
        votes[f->sweeplabel()]++;
      }
    }
    if (!votes.empty()) {
      int bestBlock = votes.begin()->first;
      int bestCount = votes.begin()->second;
      for (const auto &[bid, cnt] : votes) {
        if (cnt > bestCount) {
          bestBlock = bid;
          bestCount = cnt;
        }
      }
      applyBlockColor(v, bestBlock);
    }
  }

  std::map<int, int> faceBlockCounts;
  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); ++mfiter) {
    auto *f = static_cast<MeshLib::CToolFace *>(mfiter.value());
    if (f->sweeplabel() >= 0) {
      faceBlockCounts[f->sweeplabel()]++;
    }
  }
  std::cout << "[ApplySweepVisualization] face block distribution:";
  for (const auto &[blockId, count] : faceBlockCounts) {
    std::cout << " " << blockId << "(" << count << ")";
  }
  std::cout << std::endl;

  for (MeshLib::MeshFaceIterator mviter(this->mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolFace *face =
        static_cast<MeshLib::CToolFace *>(mviter.value());
    int sl = face->sweeplabel();
    if (sl >= 0 && sl < static_cast<int>(sweepBlockNonPlanar.size()) &&
        sweepBlockNonPlanar[static_cast<size_t>(sl)]) {
      face->sweepFaceType() = 4;
      continue;
    }

    CPoint p1 = (face->halfedge()->target()->point() -
                 face->halfedge()->source()->point());
    CPoint p2 = (face->halfedge()->he_next()->target()->point() -
                 face->halfedge()->he_next()->source()->point());
    CPoint normal = p1 ^ p2;
    if (normal.norm() < 1e-12) {
      face->sweepFaceType() = 0;
      continue;
    }
    normal /= normal.norm();
    face->normal() = normal;

    face->sweepFaceType() = 0;
    if (sl < 0) {
      continue;
    }

    int count = 0;
    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    for (MeshLib::CTMesh::FaceVertexIterator fviter(face); !fviter.end();
         ++fviter) {
      count++;
      position += Eigen::Vector3f(fviter.value()->point()[0],
                                  fviter.value()->point()[1],
                                  fviter.value()->point()[2]);
    }
    position /= std::max(count, 1);

    Eigen::Vector3f sweepDir = SweepDirectionAt(sl, position);
    if (sweepDir.norm() < 1e-8f) {
      continue;
    }
    sweepDir.normalize();
    double angle = normal[0] * sweepDir.x() + normal[1] * sweepDir.y() +
                   normal[2] * sweepDir.z();
    if (angle > 0.8) {
      face->sweepFaceType() = 1;
    } else if (std::abs(angle) < 0.1) {
      face->sweepFaceType() = 2;
    } else if (angle < -0.8) {
      face->sweepFaceType() = 3;
    }
  }
}

void DistanceField::DFS(MeshLib::CToolVertex *vert, int label) {
  if (vert->label() == label)
    vert->marked() = true;
  else
    return;
  for (MeshLib::CTMesh::VertexVertexIterator vviter(vert); !vviter.end();
       vviter++) {
    MeshLib::CToolVertex *vertex =
        static_cast<MeshLib::CToolVertex *>(vviter.value());
    if (vertex->label() == label && vertex->marked() == false) {
      DFS(vertex, label);
    }
  }
}

/**
 * @brief Checks if a point is inside the hexahedron defined by its 8
 * vertices.
 * * Uses Point-to-Plane test with corrected external normal vectors.
 */
bool DistanceField::insideCuttingBox(
    Eigen::Vector3f point,
    const std::map<int, Eigen::Vector3f> &verticesMap) const {

  if (verticesMap.size() != 8) {
    std::cerr << "Error in insideCuttingBox: Vertex map size is not 8."
              << std::endl;
    return false;
  }

  const Eigen::Vector3f &v0 = verticesMap.at(0);
  const Eigen::Vector3f &v1 = verticesMap.at(1);
  const Eigen::Vector3f &v2 = verticesMap.at(2);
  const Eigen::Vector3f &v4 = verticesMap.at(4);

  Eigen::Vector3f dirX_vec = v4 - v0;
  Eigen::Vector3f dirY_vec = v2 - v0;
  Eigen::Vector3f dirZ_vec = v1 - v0;

  const float EPSILON = 1e-6f;

  Eigen::Vector3f nX = dirY_vec.cross(dirZ_vec).normalized();
  Eigen::Vector3f nY = dirX_vec.cross(dirZ_vec).normalized();
  Eigen::Vector3f nZ = dirX_vec.cross(dirY_vec).normalized();

  bool nX_points_to_maxX = (nX.dot(dirX_vec) > 0);

  Eigen::Vector3f n_minX = nX_points_to_maxX ? -nX : nX;
  float d_minX = -n_minX.dot(v0);

  Eigen::Vector3f n_maxX = nX_points_to_maxX ? nX : -nX;
  float d_maxX = -n_maxX.dot(v4);

  if (n_minX.dot(point) + d_minX > EPSILON)
    return false; // P 在 MinX 外侧
  if (n_maxX.dot(point) + d_maxX > EPSILON)
    return false; // P 在 MaxX 外侧

  bool nY_points_to_maxY = (nY.dot(dirY_vec) > 0);

  Eigen::Vector3f n_minY = nY_points_to_maxY ? -nY : nY;
  float d_minY = -n_minY.dot(v0);

  Eigen::Vector3f n_maxY = nY_points_to_maxY ? nY : -nY;
  float d_maxY = -n_maxY.dot(v2);

  if (n_minY.dot(point) + d_minY > EPSILON)
    return false; // P 在 MinY 外侧
  if (n_maxY.dot(point) + d_maxY > EPSILON)
    return false; // P 在 MaxY 外侧
  bool nZ_points_to_maxZ = (nZ.dot(dirZ_vec) > 0);

  Eigen::Vector3f n_minZ = nZ_points_to_maxZ ? -nZ : nZ;
  float d_minZ = -n_minZ.dot(v0);

  Eigen::Vector3f n_maxZ = nZ_points_to_maxZ ? nZ : -nZ;
  float d_maxZ = -n_maxZ.dot(v1);

  if (n_minZ.dot(point) + d_minZ > EPSILON)
    return false; // P 在 MinZ 外侧
  if (n_maxZ.dot(point) + d_maxZ > EPSILON)
    return false; // P 在 MaxZ 外侧

  return true; // 如果点在所有六个平面的内侧，则在六面体内部
}

void DistanceField::readPrime(string primefile) {
  std::ifstream file(primefile);
  if (!file.is_open()) {
    std::cerr << "Error: Unable to open the file " << primefile << std::endl;
    return;
  }

  std::string line;
  PrimeData current_prime;

  while (std::getline(file, line)) {
    if (line.find("m_primes[") != std::string::npos &&
        line.find("->GetParams():") != std::string::npos) {
      int id;
      if (std::sscanf(line.c_str(), "m_primes[%d]->GetParams():", &id) != 1) {
        std::cerr << "Warn: id resolve failed,Content: " << line << std::endl;
        continue;
      }
      current_prime.id = id;
      current_prime.params.clear();

      for (int i = 0; i < 10; ++i) {
        if (!std::getline(file, line)) {
          std::cerr << "Warn: param count not match 10，id=" << id << std::endl;
          break;
        }
        double param;
        if (std::sscanf(line.c_str(), "%lf", &param) != 1) {
          std::cerr << "Warn: param resolve failed ,Content: " << line
                    << std::endl;
          param = 0.0;
        }
        current_prime.params.push_back(param);
      }

      while (current_prime.params.size() < 10) {
        current_prime.params.push_back(0.0);
      }

    } else if (line.find("of Rank:") != std::string::npos &&
               line.find("Residual:") != std::string::npos) {
      int rank;
      double residual;
      if (std::sscanf(line.c_str(), "of Rank: %d,Residual: %lf", &rank,
                      &residual) != 2) {
        std::cerr << "Warn: rank/residual resolve failed, Content: " << line
                  << std::endl;
        rank = -1;
        residual = 0.0;
      }
      current_prime.rank = rank;
      current_prime.residual = residual;
      bool isPlane =
          (current_prime.params[4] == 0 && current_prime.params[5] == 0 &&
           current_prime.params[6] == 0 && current_prime.params[7] == 0 &&
           current_prime.params[8] == 0 && current_prime.params[9] == 0);
      current_prime.isPlane = isPlane;

      primes.push_back(current_prime);
    }
  }

  file.close();
  ReindexPrimesById();

  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolVertex *v =
        static_cast<MeshLib::CToolVertex *>(mviter.value());

    int lbl = v->label();
    if (!PrimeLabelValid(lbl)) {
      continue;
    }

    bool FeaturePoint = false;
    for (MeshLib::CTMesh::VertexVertexIterator vviter(v); !vviter.end();
         ++vviter) {
      if (static_cast<MeshLib::CToolVertex *>(vviter.value())->label() !=
          v->label()) {
        FeaturePoint = true;
        break;
      }
    }
    if (FeaturePoint)
      continue;
    auto params = primes[lbl].params;
    auto vertPoint = v->point();

    Eigen::Vector3f Percise_normal = Eigen::Vector3f(0, 0, 0);
    Percise_normal[0] = params[1] + vertPoint[1] * params[4] +
                        vertPoint[2] * params[5] + 2 * params[7] * vertPoint[0];
    Percise_normal[1] = params[2] + vertPoint[0] * params[4] +
                        vertPoint[2] * params[6] + 2 * params[8] * vertPoint[1];
    Percise_normal[2] = params[3] + vertPoint[0] * params[5] +
                        vertPoint[1] * params[6] + 2 * params[9] * vertPoint[2];
    Eigen::Vector3f FormerNormal =
        Eigen::Vector3f(v->normal()[0], v->normal()[1], v->normal()[2]);
    if (FormerNormal.dot(Percise_normal) < 0)
      Percise_normal = -Percise_normal;
    v->normal()[0] = Percise_normal[0];
    v->normal()[1] = Percise_normal[1];
    v->normal()[2] = Percise_normal[2];
  }
  int maxPrimeid = static_cast<int>(this->primes.size()) - 1;

  for (int i = 0; i <= maxPrimeid; ++i) {
    if (this->primes[static_cast<size_t>(i)].id != i) {
      continue;
    }

    MeshLib::CToolVertex *startvertex = NULL;
    for (MeshLib::MeshVertexIterator mviter(this->mesh); !mviter.end();
         ++mviter) {
      MeshLib::CToolVertex *v =
          static_cast<MeshLib::CToolVertex *>(mviter.value());
      if (startvertex == NULL && v->label() == i)
        startvertex = v;
      v->marked() = false;
    }
    if (startvertex == NULL)
      continue;

    DFS(startvertex, i);
    bool hasPoped = false;
    for (MeshLib::MeshVertexIterator mviter(this->mesh); !mviter.end();
         ++mviter) {
      MeshLib::CToolVertex *v =
          static_cast<MeshLib::CToolVertex *>(mviter.value());
      if (v->label() == i && v->marked() == false) {
        v->label() = maxPrimeid + 1;
        hasPoped = true;
      }
    }
    if (hasPoped) {
      PrimeData newprime = this->primes[static_cast<size_t>(i)];
      newprime.id = maxPrimeid + 1;
      maxPrimeid++;
      if (static_cast<int>(this->primes.size()) <= maxPrimeid) {
        this->primes.resize(static_cast<size_t>(maxPrimeid) + 1);
      }
      this->primes[static_cast<size_t>(maxPrimeid)] = newprime;
    }
  }

  for (int i = 0; i < static_cast<int>(this->primes.size()); i++) {
    auto &m_params = this->primes[static_cast<size_t>(i)].params;

    std::cout << "ID: " << this->primes[static_cast<size_t>(i)].id
              << " Params: " << m_params[0]
              << " " << m_params[1] << " " << m_params[2] << " " << m_params[3]
              << " " << m_params[4] << " " << m_params[5] << " " << m_params[6]
              << " " << m_params[7] << " " << m_params[8] << " " << m_params[9]
              << std::endl;
  }
}
void DistanceField::SaveFieldToBinary(const std::string &filename) {
  std::ofstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  int xSize = Field.size();
  int ySize = (xSize > 0) ? Field[0].size() : 0;
  int zSize = (ySize > 0) ? Field[0][0].size() : 0;

  file.write(reinterpret_cast<const char *>(&xSize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&ySize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&zSize), sizeof(int));

  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      file.write(reinterpret_cast<const char *>(Field[i][j].data()),
                 zSize * sizeof(float));
    }
  }

  file.close();
  std::cout << "Field data saved to " << filename << std::endl;
}
void DistanceField::SaveGradianceToBinary(const std::string &filename) {

  std::ofstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  int xSize = Field.size();
  int ySize = (xSize > 0) ? Field[0].size() : 0;
  int zSize = (ySize > 0) ? Field[0][0].size() : 0;

  file.write(reinterpret_cast<const char *>(&xSize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&ySize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&zSize), sizeof(int));

  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      file.write(reinterpret_cast<const char *>(
                     this->getGradianceCount()[i][j].data()),
                 zSize * sizeof(float));
    }
  }

  file.close();
  std::cout << "Field data saved to " << filename << std::endl;
}

void DistanceField::RunCuttingBoxPipeline(bool cutMesh) {
  SweepProjection_Regist(cutMesh);
}

int DistanceField::ComputeSweepDirectionEnergies() {
  this->ExtractSweepDir();

  this->SweepProjScalar.clear();
  this->SweepProjEnergy.clear();
  this->sweepEnergyNames.clear();

  if (this->getGradianceCount().size() == 0) {
    return 0;
  }
  if (this->GradianceField.size() == 0) {
    return 0;
  }
  if (this->SweepDir.empty()) {
    std::cerr << "[ComputeSweepDirectionEnergies] no sweep directions\n";
    return 0;
  }

  int xSize = static_cast<int>(Field.size());
  int ySize = (xSize > 0) ? static_cast<int>(Field[0].size()) : 0;
  int zSize = (ySize > 0) ? static_cast<int>(Field[0][0].size()) : 0;

  for (int DirCount = 0; DirCount < static_cast<int>(SweepDir.size());
       DirCount++) {
    std::vector<std::vector<std::vector<float>>> ProjScalar;
    ProjScalar.reserve(static_cast<size_t>(xSize));
    for (int i = 0; i < xSize; i++) {
      std::vector<std::vector<float>> ProjScalarX;
      ProjScalarX.reserve(static_cast<size_t>(ySize));
      for (int j = 0; j < ySize; j++) {
        std::vector<float> ProjScalarXY;
        ProjScalarXY.reserve(static_cast<size_t>(zSize));
        for (int k = 0; k < zSize; k++) {
          float angle = std::acos(
              abs(this->GradianceField[i][j][k].dot(SweepDir[DirCount]) /
                  (this->GradianceField[i][j][k].norm() *
                   SweepDir[DirCount].norm())));
          float ProjScalarXYZ =
              abs(angle) > abs(PI / 2 - angle) ? abs(PI / 2 - angle)
                                               : abs(angle);
          if (this->Field[i][j][k] < 0.0f)
            ProjScalarXYZ = EXTERIOR_SWEEP_ENERGY;
          ProjScalarXY.push_back(ProjScalarXYZ);
        }
        ProjScalarX.push_back(ProjScalarXY);
      }
      ProjScalar.push_back(ProjScalarX);
    }
    this->SweepProjScalar.push_back(ProjScalar);
  }

  SweepDirFilter sf(&this->SweepDir, &this->SweepProjScalar, this->FieldLabel);
  this->SweepProjEnergy = this->SweepProjScalar;
  std::cout << "SweepDirFilter done. SweepDir size: " << this->SweepDir.size()
            << ", SweepProjEnergy size: " << this->SweepProjEnergy.size()
            << std::endl;
  SweepDirSpliter sp(this->mesh, &this->SweepDir, &this->SweepProjEnergy,
                     this->FieldLabel);
  std::cout << "SweepDirSpliter done. SweepDir size: " << this->SweepDir.size()
            << std::endl;

  // Spliter 可能增删方向，同步标量场尺寸
  if (this->SweepProjScalar.size() != this->SweepProjEnergy.size()) {
    this->SweepProjScalar = this->SweepProjEnergy;
  }

  int DirSize = static_cast<int>(this->SweepDir.size());
  auto SweepEnergy = this->SweepProjEnergy;
  float patch = (this->Coord[0][0][0] - this->Coord[0][0][1]).norm();
  STEP_SIZE = patch;
  std::cout << "Starting energy computation..." << std::endl;
  for (int dirs = 0; dirs < static_cast<int>(this->SweepProjEnergy.size());
       dirs++) {
    for (int x = 0; x < static_cast<int>(this->SweepProjEnergy[dirs].size());
         x++) {
      for (int y = 0;
           y < static_cast<int>(this->SweepProjEnergy[dirs][x].size()); y++) {
        for (int z = 0;
             z < static_cast<int>(this->SweepProjEnergy[dirs][x][y].size());
             z++) {
          if (this->Field[x][y][z] < 0.0f) {
            SweepEnergy[dirs][x][y][z] = EXTERIOR_SWEEP_ENERGY;
            continue;
          }
          int fl = this->FieldLabel[x][y][z];
          if (fl < 0) {
            SweepEnergy[dirs][x][y][z] = -2e-4;
            continue;
          }
          const PrimeData *prime = GetPrimeByLabel(fl);
          if (prime && prime->isPlane &&
              abs(this->Field[x][y][z]) < 2 * patch) {
            SweepEnergy[dirs][x][y][z] = -2e-3;
            continue;
          }
          float OtherEnergy = -2e-3;
          for (int RestDirs = 0;
               RestDirs < static_cast<int>(this->SweepProjEnergy.size());
               RestDirs++) {
            if (RestDirs != dirs)
              OtherEnergy += -this->SweepProjEnergy[RestDirs][x][y][z];
          }
          float SelfEnergy = this->SweepProjEnergy[dirs][x][y][z];
          SelfEnergy *= (DirSize - 2);
          SweepEnergy[dirs][x][y][z] =
              Alpha * SelfEnergy + (1 - Alpha) * OtherEnergy;
        }
      }
    }
  }
  this->SweepProjEnergy = SweepEnergy;

  // 能量与方向必须一一对应
  if (this->SweepDir.size() != this->SweepProjEnergy.size()) {
    std::cerr << "[ComputeSweepDirectionEnergies] size mismatch SweepDir="
              << this->SweepDir.size() << " Energy="
              << this->SweepProjEnergy.size() << "\n";
    const size_t n =
        std::min(this->SweepDir.size(), this->SweepProjEnergy.size());
    this->SweepDir.resize(n);
    this->SweepProjEnergy.resize(n);
    this->SweepProjScalar.resize(n);
  }

  this->sweepEnergyNames.clear();
  this->sweepEnergyNames.reserve(this->SweepProjEnergy.size());
  for (size_t i = 0; i < this->SweepProjEnergy.size(); ++i) {
    this->sweepEnergyNames.push_back("Sweep Energy " + std::to_string(i));
  }
  std::cout << "[ComputeSweepDirectionEnergies] registered "
            << this->SweepProjEnergy.size()
            << " directional energy fields (= #boxes to build)\n";
  return static_cast<int>(this->SweepProjEnergy.size());
}

void DistanceField::SweepProjection_Regist(bool cutMesh) {
  // 1) 算能量  2) 有几个能量就建几个框
  if (ComputeSweepDirectionEnergies() <= 0) {
    return;
  }
  BuildCuttingBoxesFromEnergies(cutMesh);
}

void DistanceField::BuildCuttingBoxesFromEnergies(bool cutMesh) {
  const int nEnergy = static_cast<int>(this->SweepProjEnergy.size());
  if (nEnergy <= 0) {
    std::cerr << "[BuildCuttingBoxesFromEnergies] no energy fields; call "
                 "ComputeSweepDirectionEnergies first\n";
    return;
  }
  if (static_cast<int>(this->SweepDir.size()) != nEnergy) {
    std::cerr << "[BuildCuttingBoxesFromEnergies] SweepDir("
              << this->SweepDir.size() << ") != SweepProjEnergy(" << nEnergy
              << "); syncing to energy count\n";
    if (static_cast<int>(this->SweepDir.size()) > nEnergy) {
      this->SweepDir.resize(static_cast<size_t>(nEnergy));
    } else {
      while (static_cast<int>(this->SweepDir.size()) < nEnergy) {
        this->SweepDir.push_back(Eigen::Vector3f::UnitY());
      }
    }
  }
  if (static_cast<int>(this->SweepProjScalar.size()) != nEnergy) {
    this->SweepProjScalar = this->SweepProjEnergy;
  }
  if (static_cast<int>(this->sweepEnergyNames.size()) != nEnergy) {
    this->sweepEnergyNames.resize(static_cast<size_t>(nEnergy));
    for (int i = 0; i < nEnergy; ++i) {
      if (this->sweepEnergyNames[static_cast<size_t>(i)].empty()) {
        this->sweepEnergyNames[static_cast<size_t>(i)] =
            "Sweep Energy " + std::to_string(i);
      }
    }
  }

  this->CuttingHexLists.clear();
  this->sweepBlockNonPlanar.clear();
  this->sweepBlockColors.clear();

  std::cout << "[BuildCuttingBoxesFromEnergies] building " << nEnergy
            << " CuttingBoxes (1 per energy / sweep dir)\n";
  for (int i = 0; i < nEnergy; ++i) {
    CuttingBox cb(this->SweepDir, &this->SweepProjEnergy, this->Coord,
                  this->FieldLabel, this->Field, this->primes, i);
    this->ForbiddenBoundaryPoints = cb.GetForbiddenBoundaryPoints();
    this->CuttingHexLists.push_back(cb.GetBoxVertices());
    this->sweepBlockNonPlanar.push_back(false);
    this->sweepBlockColors.push_back(RandomSweepColor(i));
    std::cout << "  box " << i << " dir=(" << this->SweepDir[static_cast<size_t>(i)].transpose()
              << ") energy=\"" << this->sweepEnergyNames[static_cast<size_t>(i)]
              << "\"\n";
  }
  std::cout << "CuttingBoxes done (" << this->CuttingHexLists.size() << ")"
            << std::endl;

  if (cutMesh && this->mesh && !this->CuttingHexLists.empty()) {
    std::cout << "Running MeshCutter..." << std::endl;
    MeshCutter mc(this->mesh, this->CuttingHexLists);
    std::cout << "MeshCutter done." << std::endl;
  }
}

static std::vector<std::map<int, Eigen::Vector3f>>
FilterValidHexes(const std::vector<std::map<int, Eigen::Vector3f>> &hexes,
                 float stepSize) {
  std::vector<std::map<int, Eigen::Vector3f>> validHexes;
  for (const auto &hex : hexes) {
    if (hex.size() != 8) {
      continue;
    }
    float minEdge = std::numeric_limits<float>::max();
    for (int i = 0; i < 8; ++i) {
      for (int j = i + 1; j < 8; ++j) {
        minEdge = std::min(minEdge, (hex.at(i) - hex.at(j)).norm());
      }
    }
    if (minEdge > 1e-4f * stepSize) {
      validHexes.push_back(hex);
    }
  }
  return validHexes;
}

void DistanceField::GeneralizedSweepDecomposition(float angularThreshold,
                                                bool cylinderPairsOnly) {
  if (this->Field.empty() || this->GradianceField.empty() ||
      this->primes.empty()) {
    std::cerr << "[GeneralizedSweepDecomposition] Prerequisite data missing. "
              << "Ensure ComputeDistanceField() and readPrime() are called first."
              << std::endl;
    return;
  }

  std::cout << "[GeneralizedSweepDecomposition] threshold=" << angularThreshold
            << " rad"
            << (cylinderPairsOnly ? " (cylinder pairs only)" : "") << std::endl;

  if (!cylinderPairsOnly) {
    this->CuttingHexLists.clear();
    this->sweepBlocks.clear();
    this->SweepDir.clear();
    this->sweepBlockNonPlanar.clear();
    this->sweepBlockColors.clear();
  }

  float stepSize = (Coord[0][0][0] - Coord[0][0][1]).norm();
  SweepDecomposer decomposer(this->Coord, this->Field, this->FieldLabel,
                             this->GradianceField, this->primes,
                             angularThreshold, cylinderPairsOnly, this->mesh);

  auto blocks = decomposer.GetBlocks();
  if (cylinderPairsOnly) {
    this->cylinderPairViz = decomposer.GetCylinderPairViz();
  }
  auto hexes = FilterValidHexes(decomposer.GetBlockHexVertices(), stepSize);
  if (blocks.size() != hexes.size()) {
    size_t n = std::min(blocks.size(), hexes.size());
    blocks.resize(n);
    hexes.resize(n);
  }

  bool nonPlanar = cylinderPairsOnly;
  if (!cylinderPairsOnly) {
    for (const auto &block : blocks) {
      nonPlanar = block.kind == SweepKind::CylindricalBase;
      if (nonPlanar) {
        break;
      }
    }
  }

  AppendSweepBlocks(blocks, hexes, nonPlanar);

  int nonPlanarCount = 0;
  for (bool flag : this->sweepBlockNonPlanar) {
    if (flag) {
      nonPlanarCount++;
    }
  }
  std::cout << "[GeneralizedSweepDecomposition] " << blocks.size()
            << " blocks appended, total hexes "
            << this->CuttingHexLists.size() << " (non-planar "
            << nonPlanarCount << ", base planar "
            << FindBasePlanarHexIndex() << ")\n";

  if (!cylinderPairsOnly && !this->CuttingHexLists.empty() && this->mesh) {
    MeshCutter mc(this->mesh, this->CuttingHexLists);
    std::cout << "[GeneralizedSweepDecomposition] MeshCutter done." << std::endl;
  }
}

void DistanceField::AppendCylinderSweepDecomposition(float angularThreshold) {
  std::cout << "[AppendCylinderSweep] non-planar sweep for gradient-matched "
               "cylinder pairs only\n";
  // 柱面扫掠也必须先有方向能量场（与平面 CuttingBox 同一套）
  if (this->SweepProjEnergy.empty()) {
    std::cout << "[AppendCylinderSweep] computing directional sweep energies "
                 "first\n";
    ComputeSweepDirectionEnergies();
  }
  GeneralizedSweepDecomposition(angularThreshold, true);
}

VoxelIndex DistanceField::WorldToVoxel(const Eigen::Vector3f &p) const {
  const Eigen::Vector3f &o = Coord[0][0][0];
  int D1 = static_cast<int>(Coord.size());
  int D2 = D1 > 0 ? static_cast<int>(Coord[0].size()) : 0;
  int D3 = D2 > 0 ? static_cast<int>(Coord[0][0].size()) : 0;
  int i = static_cast<int>(std::lround((p.x() - o.x()) / PatchSize));
  int j = static_cast<int>(std::lround((p.y() - o.y()) / PatchSize));
  int k = static_cast<int>(std::lround((p.z() - o.z()) / PatchSize));
  i = std::max(0, std::min(D1 - 1, i));
  j = std::max(0, std::min(D2 - 1, j));
  k = std::max(0, std::min(D3 - 1, k));
  return {i, j, k};
}

std::map<int, Eigen::Vector3f>
DistanceField::BuildOrientedHex(const std::vector<Eigen::Vector3f> &pts,
                                const Eigen::Vector3f &axis,
                                const Eigen::Vector3f &origin) const {
  Eigen::Vector3f ax = axis.normalized();
  Eigen::Vector3f arbitrary =
      (std::abs(ax.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  Eigen::Vector3f crossY = (arbitrary - ax.dot(arbitrary) * ax).normalized();
  Eigen::Vector3f crossZ = ax.cross(crossY).normalized();

  float axMin = std::numeric_limits<float>::max();
  float axMax = std::numeric_limits<float>::lowest();
  float minY = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::lowest();
  float minZ = std::numeric_limits<float>::max();
  float maxZ = std::numeric_limits<float>::lowest();
  for (const auto &p : pts) {
    Eigen::Vector3f rel = p - origin;
    float a = rel.dot(ax);
    float cy = rel.dot(crossY);
    float cz = rel.dot(crossZ);
    axMin = std::min(axMin, a);
    axMax = std::max(axMax, a);
    minY = std::min(minY, cy);
    maxY = std::max(maxY, cy);
    minZ = std::min(minZ, cz);
    maxZ = std::max(maxZ, cz);
  }

  const float margin = std::max(PatchSize, 0.05f);
  axMin -= margin;
  axMax += margin;
  minY -= margin;
  maxY += margin;
  minZ -= margin;
  maxZ += margin;

  auto corner = [&](float a, float cy, float cz) {
    return origin + a * ax + cy * crossY + cz * crossZ;
  };

  std::map<int, Eigen::Vector3f> vertices;
  float axVals[] = {axMin, axMax};
  float yVals[] = {minY, maxY};
  float zVals[] = {minZ, maxZ};
  int idx = 0;
  for (float a : axVals) {
    for (float cy : yVals) {
      for (float cz : zVals) {
        vertices[idx++] = corner(a, cy, cz);
      }
    }
  }
  return vertices;
}

namespace {
std::vector<std::vector<std::vector<float>>> ComputeDirectionalSweepEnergy(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField,
    const std::vector<std::vector<std::vector<float>>> &Field,
    const Eigen::Vector3f &sweepDir) {
  const int D1 = static_cast<int>(Field.size());
  const int D2 = D1 > 0 ? static_cast<int>(Field[0].size()) : 0;
  const int D3 = D2 > 0 ? static_cast<int>(Field[0][0].size()) : 0;
  std::vector<std::vector<std::vector<float>>> grid(
      D1, std::vector<std::vector<float>>(
              D2, std::vector<float>(D3, EXTERIOR_SWEEP_ENERGY)));
  const Eigen::Vector3f dir = sweepDir.normalized();
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] < 0.0f) {
          continue;
        }
        const Eigen::Vector3f &g = GradField[x][y][z];
        if (g.norm() < 1e-8f) {
          grid[x][y][z] = static_cast<float>(PI) / 4.0f;
          continue;
        }
        float c = std::min(
            1.0f, std::abs(g.normalized().dot(dir)));
        float angle = std::acos(c);
        grid[x][y][z] =
            std::min(angle, std::abs(static_cast<float>(PI) / 2.0f - angle));
      }
    }
  }
  return grid;
}
} // namespace

void DistanceField::DecomposeIntoTwoSweepBodies(float angularThreshold) {
  std::cout << "[TwoSweepBodies] decomposing into cylinder-radial (tube) + "
               "vertical (cube) bodies\n";

  // --- 0. 柱面扫掠也首先计算各方向扫掠能量（与平面 CuttingBox 同一套） ---
  const int nDirEnergy = ComputeSweepDirectionEnergies();
  auto planarScalar = this->SweepProjScalar;
  auto planarEnergy = this->SweepProjEnergy;
  auto planarEnergyNames = this->sweepEnergyNames;
  auto planarSweepDirs = this->SweepDir;
  std::cout << "[TwoSweepBodies] kept " << nDirEnergy
            << " planar directional energies for VolumeGrid\n";

  // --- 1. 再用柱面分解拿到柱轴/原点/成员 prime（仅取几何参数） ---
  GeneralizedSweepDecomposition(angularThreshold, /*cylinderPairsOnly=*/true);

  Eigen::Vector3f axis = Eigen::Vector3f::UnitY();
  Eigen::Vector3f origin = Eigen::Vector3f::Zero();
  std::vector<int> memberPrimeIds;
  if (!this->sweepBlocks.empty()) {
    axis = this->sweepBlocks.front().sweepAxis.normalized();
    origin = this->sweepBlocks.front().sweepOrigin;
    memberPrimeIds = this->sweepBlocks.front().memberPrimeIds;
  }
  if (axis.norm() < 1e-6f) {
    axis = Eigen::Vector3f::UnitY();
  }

  // 清空分块状态，重建两个扫掠体；方向能量场稍后恢复，不丢弃
  this->CuttingHexLists.clear();
  this->sweepBlocks.clear();
  this->SweepDir.clear();
  this->sweepBlockNonPlanar.clear();
  this->sweepBlockColors.clear();
  this->cylinderPairViz.clear();
  this->voxelToBlock.clear();

  auto radialDist = [&](const Eigen::Vector3f &p) {
    Eigen::Vector3f rel = p - origin;
    return (rel - rel.dot(axis) * axis).norm();
  };
  auto axialPos = [&](const Eigen::Vector3f &p) {
    return (p - origin).dot(axis);
  };

  // --- 2. 从网格统计每片柱面 prime 的半径与轴向范围，定位“管壁”特征 ---
  std::set<int> memberSet(memberPrimeIds.begin(), memberPrimeIds.end());
  std::map<int, float> primeRadiusSum;
  std::map<int, int> primeCount;
  std::map<int, float> primeAxMin;
  std::map<int, float> primeAxMax;
  if (this->mesh) {
    for (MeshLib::MeshVertexIterator viter(this->mesh); !viter.end(); ++viter) {
      auto *v = static_cast<MeshLib::CToolVertex *>(viter.value());
      int lbl = v->label();
      if (!memberSet.count(lbl)) {
        continue;
      }
      Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
      primeRadiusSum[lbl] += radialDist(p);
      float a = axialPos(p);
      if (!primeCount.count(lbl)) {
        primeAxMin[lbl] = a;
        primeAxMax[lbl] = a;
      } else {
        primeAxMin[lbl] = std::min(primeAxMin[lbl], a);
        primeAxMax[lbl] = std::max(primeAxMax[lbl], a);
      }
      primeCount[lbl]++;
    }
  }

  // 管壁 = 轴向跨度大的柱面 prime（区别于柱孔底盖等短小成员）。
  float maxExtent = 0.0f;
  for (const auto &[lbl, axMin] : primeAxMin) {
    maxExtent = std::max(maxExtent, primeAxMax[lbl] - axMin);
  }
  int outerPrime = -1;
  float rOuter = 0.0f;
  float rInner = std::numeric_limits<float>::max();
  std::set<int> wallPrimes;
  for (const auto &[lbl, sum] : primeRadiusSum) {
    int cnt = primeCount[lbl];
    if (cnt <= 0) {
      continue;
    }
    float r = sum / static_cast<float>(cnt);
    float extent = primeAxMax[lbl] - primeAxMin[lbl];
    bool isWall = extent >= 0.4f * maxExtent;
    std::cout << "[TwoSweepBodies] cylinder prime " << lbl << " avgR=" << r
              << " ax=[" << primeAxMin[lbl] << ", " << primeAxMax[lbl]
              << "] extent=" << extent << (isWall ? " [wall]" : "") << "\n";
    if (!isWall) {
      continue;
    }
    wallPrimes.insert(lbl);
    if (r > rOuter) {
      rOuter = r;
      outerPrime = lbl;
    }
    rInner = std::min(rInner, r);
  }
  if (rInner > rOuter) {
    rInner = 0.0f;
  }

  // 在已算好的方向能量中，找与柱轴最对齐的那一路
  int axisEnergyIdx = -1;
  float bestAlign = -1.0f;
  for (int i = 0; i < static_cast<int>(planarSweepDirs.size()); ++i) {
    float align = std::abs(planarSweepDirs[static_cast<size_t>(i)]
                               .normalized()
                               .dot(axis));
    if (align > bestAlign) {
      bestAlign = align;
      axisEnergyIdx = i;
    }
  }
  if (axisEnergyIdx >= 0) {
    std::cout << "[TwoSweepBodies] using planar Sweep Energy " << axisEnergyIdx
              << " (align=" << bestAlign << ") with cylinder axis\n";
  }

  float tubeAxLow = std::numeric_limits<float>::lowest();
  float tubeAxHigh = std::numeric_limits<float>::max();
  std::map<int, Eigen::Vector3f> tubeHexFromCut;
  std::vector<std::vector<std::vector<float>>> radialEnergyField;
  bool haveTubeCut = false;
  bool haveRadialEnergy = false;
  if (outerPrime >= 0) {
    CylinderCuttingBox cyl(axis, origin, rInner, rOuter, this->Coord,
                           this->Field, this->FieldLabel, this->GradianceField,
                           wallPrimes, angularThreshold);
    rInner = cyl.GetMinR();
    rOuter = cyl.GetMaxR();
    tubeAxLow = cyl.GetMinAx();
    tubeAxHigh = cyl.GetMaxAx();
    tubeHexFromCut = cyl.GetBoxVertices();
    radialEnergyField = cyl.ComputeRadialEnergyField();
    haveTubeCut = true;
    haveRadialEnergy = true;
  }

  int D1 = static_cast<int>(Field.size());
  int D2 = D1 > 0 ? static_cast<int>(Field[0].size()) : 0;
  int D3 = D2 > 0 ? static_cast<int>(Field[0][0].size()) : 0;

  SweepBlockRegion tubeBody;
  tubeBody.kind = SweepKind::CylindricalBase;
  tubeBody.sweepAxis = axis;
  tubeBody.sweepOrigin = origin;
  tubeBody.memberPrimeIds = memberPrimeIds;
  tubeBody.primeId = memberPrimeIds.empty() ? -1 : memberPrimeIds.front();
  tubeBody.isValid = true;

  SweepBlockRegion cubeBody;
  cubeBody.kind = SweepKind::Translational;
  cubeBody.sweepAxis = axis;
  cubeBody.primeId = -1;
  cubeBody.isValid = true;

  std::vector<Eigen::Vector3f> tubePts;
  std::vector<Eigen::Vector3f> cubePts;
  Eigen::Vector3f tubeCen = Eigen::Vector3f::Zero();
  Eigen::Vector3f cubeCen = Eigen::Vector3f::Zero();

  const float rBand = 2.0f * PatchSize;
  const bool useAxisEnergy =
      axisEnergyIdx >= 0 &&
      axisEnergyIdx < static_cast<int>(planarEnergy.size());
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] <= 0.0f) {
          continue;
        }
        const Eigen::Vector3f &p = Coord[x][y][z];
        float a = axialPos(p);
        float r = radialDist(p);
        VoxelIndex vi{x, y, z};
        bool inTubeBand = (a >= tubeAxLow && a <= tubeAxHigh);
        bool inTubeWall =
            haveTubeCut && inTubeBand && r <= rOuter + rBand;
        if (inTubeWall && useAxisEnergy) {
          float e = planarEnergy[static_cast<size_t>(axisEnergyIdx)][x][y][z];
          if (e > angularThreshold * 2.0f && r < rInner - rBand) {
            inTubeWall = false;
          }
        }
        if (inTubeWall) {
          tubeBody.coveredVoxels.insert(vi);
          tubePts.push_back(p);
          tubeCen += p;
        } else {
          cubeBody.coveredVoxels.insert(vi);
          cubePts.push_back(p);
          cubeCen += p;
        }
      }
    }
  }

  if (!tubePts.empty()) {
    tubeCen /= static_cast<float>(tubePts.size());
    tubeBody.radialInner = std::max(0.0f, rInner);
    tubeBody.radialOuter = rOuter;
    tubeBody.axialLower = tubeAxLow;
    tubeBody.axialUpper = tubeAxHigh;
    std::map<int, Eigen::Vector3f> tubeHex =
        haveTubeCut ? tubeHexFromCut : BuildOrientedHex(tubePts, axis, origin);
    this->sweepBlocks.push_back(tubeBody);
    this->SweepDir.push_back(axis);
    this->sweepBlockNonPlanar.push_back(true);
    this->sweepBlockColors.push_back(
        RandomSweepColor(static_cast<int>(this->sweepBlockColors.size())));
    this->CuttingHexLists.push_back(tubeHex);
  }

  if (!cubePts.empty()) {
    cubeCen /= static_cast<float>(cubePts.size());
    cubeBody.sweepOrigin = cubeCen;
    std::map<int, Eigen::Vector3f> cubeHex =
        BuildOrientedHex(cubePts, axis, cubeCen);
    this->sweepBlocks.push_back(cubeBody);
    this->SweepDir.push_back(axis);
    this->sweepBlockNonPlanar.push_back(false);
    this->sweepBlockColors.push_back(
        RandomSweepColor(static_cast<int>(this->sweepBlockColors.size())));
    this->CuttingHexLists.push_back(cubeHex);
  }

  for (int i = static_cast<int>(this->sweepBlocks.size()) - 1; i >= 0; --i) {
    for (const auto &v :
         this->sweepBlocks[static_cast<size_t>(i)].coveredVoxels) {
      voxelToBlock[v] = i;
    }
  }

  std::cout << "[TwoSweepBodies] body count=" << this->sweepBlocks.size()
            << " (tube voxels=" << tubePts.size()
            << ", cube voxels=" << cubePts.size() << ") outerPrime="
            << outerPrime << " rInner=" << rInner << " rOuter=" << rOuter
            << " tubeAx=[" << tubeAxLow << ", " << tubeAxHigh << "]\n";

  // --- 6. 恢复平面方向能量，再追加柱面专用能量（不再清空） ---
  this->SweepProjScalar = planarScalar;
  this->SweepProjEnergy = planarEnergy;
  this->sweepEnergyNames = planarEnergyNames;
  if (haveRadialEnergy) {
    this->SweepProjScalar.push_back(radialEnergyField);
    this->SweepProjEnergy.push_back(radialEnergyField);
    this->sweepEnergyNames.push_back("Cylinder Radial Energy");
  }
  auto verticalEnergyField =
      ComputeDirectionalSweepEnergy(this->GradianceField, this->Field, axis);
  this->SweepProjScalar.push_back(verticalEnergyField);
  this->SweepProjEnergy.push_back(verticalEnergyField);
  this->sweepEnergyNames.push_back("Vertical Sweep Energy");
  std::cout << "[TwoSweepBodies] energy fields registered: "
            << this->sweepEnergyNames.size() << " (planar " << nDirEnergy
            << " + cylinder extras)\n";
}

void DistanceField::ReducePrimesToDevelopable(double curvatureThreshold,
                                              const std::string &exportPath) {
  if (this->primes.empty() || !this->mesh) {
    std::cerr << "[ReducePrimesToDevelopable] Need mesh and primes (readPrime "
                 "first).\n";
    return;
  }

  DevelopableReducer::Config cfg;
  cfg.curvatureThreshold = curvatureThreshold;
  DevelopableReducer::ReduceAll(this->primes, this->mesh, cfg);

  if (!exportPath.empty()) {
    std::ofstream out(exportPath);
    if (!out.is_open()) {
      std::cerr << "[ReducePrimesToDevelopable] Cannot write " << exportPath
                << "\n";
      return;
    }
    for (const auto &prime : this->primes) {
      out << "m_primes[" << prime.id << "]->GetParams():\n";
      for (int i = 0; i < 10; ++i) {
        double v = (i < static_cast<int>(prime.params.size())) ? prime.params[i]
                                                               : 0.0;
        out << v << "\n";
      }
      out << "of Rank: " << prime.rank << ",Residual: " << prime.residual
          << "\n";
    }
    out.close();
    std::cout << "[ReducePrimesToDevelopable] Wrote " << exportPath << "\n";
  }
}

void DistanceField::GenerateSweepHexMeshes(int divisionsU, int divisionsV,
                                           int divisionsW,
                                           float targetCellSize) {
  SweepHexMesherConfig cfg;
  cfg.divisionsU = divisionsU;
  cfg.divisionsV = divisionsV;
  cfg.divisionsW = divisionsW;
  cfg.targetCellSize = targetCellSize; // 0 = 使用显式 divisions，不自动加密

  std::vector<SweepBlockRegion> blocks = this->sweepBlocks;
  const auto &hexes = this->CuttingHexLists;

  // 平面 CuttingBox 路径没有 sweepBlocks 时，从切割盒合成平移扫掠块（仅 hex 路径使用）。
  if (blocks.empty() && !hexes.empty()) {
    std::cout << "[GenerateSweepHexMeshes] synthesizing " << hexes.size()
              << " translational blocks from planar CuttingBoxes\n";
    blocks.reserve(hexes.size());
    for (size_t i = 0; i < hexes.size(); ++i) {
      SweepBlockRegion b{};
      b.kind = SweepKind::Translational;
      b.isValid = hexes[i].size() == 8;
      b.patchId = static_cast<int>(i);
      b.primeId = -1;
      Eigen::Vector3f axis = Eigen::Vector3f::UnitY();
      if (i < this->SweepDir.size() &&
          this->SweepDir[i].norm() > 1e-8f) {
        axis = this->SweepDir[i].normalized();
      }
      b.sweepAxis = axis;
      Eigen::Vector3f origin = Eigen::Vector3f::Zero();
      for (const auto &kv : hexes[i]) {
        origin += kv.second;
      }
      if (!hexes[i].empty()) {
        origin /= static_cast<float>(hexes[i].size());
      }
      b.sweepOrigin = origin;
      SweepCapFrame frame;
      if (SweepFaceImprinter::BuildFrameFromHex(hexes[i], axis, origin,
                                                frame)) {
        b.axialLower = frame.axMin;
        b.axialUpper = frame.axMax;
        b.crossDirY = frame.crossY;
        b.crossDirZ = frame.crossZ;
      }
      blocks.push_back(b);
    }
  }

  this->sweepHexMeshes =
      SweepHexMesher::Generate(blocks, hexes, cfg, this->mesh);
}

bool DistanceField::WriteSweepHexMeshesVTK(const std::string &path) const {
  return SweepHexMesher::WriteVTK(this->sweepHexMeshes, path);
}
