#include "DFContainer.h"
#include "CTMesh.h"
#include "SweepDirDetector.h"
#include "SweepDirFilter.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

DistanceField::DistanceField() { this->primes.clear(); };

DistanceField::DistanceField(MeshLib::CTMesh *mesh) {
  this->primes.clear();
  SetMesh(mesh);
}

void DistanceField::SetMesh(MeshLib::CTMesh *mesh) {
  this->mesh = mesh;
  this->PointList.clear();
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
    std::vector<float> pointcoord;
    pointcoord.clear();
    pointcoord.push_back(viter.value()->point()[0]);
    pointcoord.push_back(viter.value()->point()[1]);
    pointcoord.push_back(viter.value()->point()[2]);
    this->PointList.push_back(pointcoord);
    Eigen::Vector3f normal = Eigen::Vector3f(0.0, 0.0, 0.0);
    MeshLib::CTMesh::CVertex *v = *viter;
    for (MeshLib::CTMesh::VertexFaceIterator vfiter(v); !vfiter.end();
         vfiter++) {
      MeshLib::CTMesh::CFace *f = *vfiter;
      if (f->area() < 1e-10)
        continue;
      normal += Eigen::Vector3f(f->normal()[0] * f->area(),
                                f->normal()[1] * f->area(),
                                f->normal()[2] * f->area());
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

void DistanceField::FindNearestPointsInOctree(
    const Eigen::Vector3f &point, std::shared_ptr<OctreeNode> node,
    std::vector<int> &candidateIndices) {
  if (!node)
    return;

  Eigen::Vector3f diff = (point - node->center).cwiseAbs();
  float distToNode =
      (diff - Eigen::Vector3f::Constant(node->halfSize)).cwiseMax(0.0f).norm();

  if (node->isLeaf) {
    candidateIndices.insert(candidateIndices.end(), node->pointIndices.begin(),
                            node->pointIndices.end());
    return;
  }

  for (const auto &child : node->children) {
    if (child) {
      FindNearestPointsInOctree(point, child, candidateIndices);
    }
  }
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

Eigen::Vector4f DistanceField::DistanceToMesh(int x, int y, int z) {
  Eigen::Vector3f point = this->Coord[x][y][z];
  Eigen::Vector4f DistanceVector;
  if (!mesh || PointList.empty()) {
    DistanceVector.setZero();
    return DistanceVector;
  }

  std::vector<int> candidateIndices;
  FindNearestPointsInOctree(point, octreeRoot, candidateIndices);

  if (candidateIndices.empty()) {
    DistanceVector.setZero();
    return DistanceVector;
  }
  std::vector<float> nearestPoint;
  nearestPoint.resize(3);
  float minDistance = std::numeric_limits<float>::max();
  for (int idx : candidateIndices) {
    const auto &meshPoint = PointList[idx];
    Eigen::Vector3f p(meshPoint[0], meshPoint[1], meshPoint[2]);
    float dist = (point - p).norm();
    if (dist < minDistance) {
      nearestPoint[0] = p[0];
      nearestPoint[1] = p[1];
      nearestPoint[2] = p[2];
      minDistance = dist;
    }
  }

  if (!this->primes.empty()) {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      MeshLib::CToolVertex *v =
          static_cast<MeshLib::CToolVertex *>(viter.value());
      if (abs(v->point()[0] - nearestPoint[0]) < 1e-16 &&
          abs(v->point()[1] - nearestPoint[1]) < 1e-16 &&
          abs(v->point()[2] - nearestPoint[2]) < 1e-16) {
        bool FeaturePoint = false;
        for (MeshLib::CTMesh::VertexVertexIterator vviter(v); !vviter.end();
             ++vviter) {
          if (static_cast<MeshLib::CToolVertex *>(vviter.value())->label() !=
              v->label()) {
            FeaturePoint = true;
            break;
          }
        }
        if (FeaturePoint) {
          Eigen::Vector3f pointOnMesh =
              Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
          Eigen::Vector3f pointOnGrid = this->Coord[x][y][z];
          Eigen::Vector3f NormalOnMesh =
              Eigen::Vector3f(v->normal()[0], v->normal()[1], v->normal()[2]);
          NormalOnMesh.normalize();

          if ((pointOnMesh - pointOnGrid).dot(NormalOnMesh) < 0)
            minDistance = -abs(minDistance);
          else
            minDistance = abs(minDistance);

          DistanceVector[0] = 0;
          DistanceVector[1] = 0;
          DistanceVector[2] = 0;
          DistanceVector[3] = minDistance;
          this->FieldLabel[x][y][z] = -1;

          break;
        }
        minDistance = this->DisCompute(point, v->label());

        this->FieldLabel[x][y][z] = v->label();
        Eigen::Vector3f VertexNormal =
            Eigen::Vector3f(v->normal()[0], v->normal()[1], v->normal()[2]);
        Eigen::Vector3f vertexPoint =
            Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
        if ((vertexPoint - point).dot(VertexNormal) < 0) {
          DistanceVector[0] = -v->normal()[0];
          DistanceVector[1] = -v->normal()[1];
          DistanceVector[2] = -v->normal()[2];
          DistanceVector[3] = -minDistance;
        } else {
          DistanceVector[0] = v->normal()[0];
          DistanceVector[1] = v->normal()[1];
          DistanceVector[2] = v->normal()[2];
          DistanceVector[3] = minDistance;
        }
        break;
      }
    }
  } else {
    for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
      MeshLib::CToolVertex *v =
          static_cast<MeshLib::CToolVertex *>(viter.value());
      if (abs(v->point()[0] - nearestPoint[0]) < 1e-16 &&
          abs(v->point()[1] - nearestPoint[1]) < 1e-16 &&
          abs(v->point()[2] - nearestPoint[2]) < 1e-16) {
        Eigen::Vector3f pointOnMesh =
            Eigen::Vector3f(v->point()[0], v->point()[1], v->point()[2]);
        Eigen::Vector3f pointOnGrid = this->Coord[x][y][z];
        Eigen::Vector3f NormalOnMesh =
            Eigen::Vector3f(v->normal()[0], v->normal()[1], v->normal()[2]);

        if ((pointOnMesh - pointOnGrid).dot(NormalOnMesh) < 0)
          minDistance = -minDistance;
        DistanceVector[0] = v->normal()[0];
        DistanceVector[1] = v->normal()[1];
        DistanceVector[2] = v->normal()[2];
        DistanceVector[3] = minDistance;
        break;
      }
    }
  }

  return DistanceVector;
}

double DistanceField::DisCompute(Eigen::Vector3f point, int label) {
  auto &m_params = this->primes[label].params;

  if (this->primes[label].isPlane) {
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

  BuildOctree();

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif // ENABLE_OMP
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        Eigen::Vector4f distance = DistanceToMesh(i, j, k);
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
#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif // ENABLE_OMP
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
        if (currentLabel == -1) {
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
  std::cout << "CrossField Build" << std::endl;
  if (this->primes.size() != 0)
    this->SweepProjection_Regist();
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

  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolVertex *v =
        static_cast<MeshLib::CToolVertex *>(mviter.value());

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
    auto params = primes[v->label()].params;
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

  for (int i = 0; i < this->primes.size(); i++) {
    auto &m_params = this->primes[i].params;

    std::cout << "Params: " << m_params[0] << " " << m_params[1] << " "
              << m_params[2] << " " << m_params[3] << " " << m_params[4] << " "
              << m_params[5] << " " << m_params[6] << " " << m_params[7] << " "
              << m_params[8] << " " << m_params[9] << std::endl;
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

void DistanceField::SweepProjection_Regist() {
  this->ExtractSweepDir();

  if (this->getGradianceCount().size() == 0) {
    return;
  }
  if (this->GradianceField.size() == 0) {
    return;
  }

  int xSize = Field.size();
  int ySize = (xSize > 0) ? Field[0].size() : 0;
  int zSize = (ySize > 0) ? Field[0][0].size() : 0;

  for (int DirCount = 0; DirCount < SweepDir.size(); DirCount++) {
    std::vector<std::vector<std::vector<float>>> ProjScalar;
    for (int i = 0; i < xSize; i++) {
      std::vector<std::vector<float>> ProjScalarX;
      for (int j = 0; j < ySize; j++) {
        std::vector<float> ProjScalarXY;
        for (int k = 0; k < zSize; k++) {
          float ProjScalarXYZ;
          float angle = std::acos(
              abs(this->GradianceField[i][j][k].dot(SweepDir[DirCount]) /
                  (this->GradianceField[i][j][k].norm() *
                   SweepDir[DirCount].norm())));
          ProjScalarXYZ = abs(angle) > abs(PI / 2 - angle) ? abs(PI / 2 - angle)
                                                           : abs(angle);
          if (this->Field[i][j][k] < 0)
            ProjScalarXYZ = 0;
          ProjScalarXY.push_back(ProjScalarXYZ);
        }
        ProjScalarX.push_back(ProjScalarXY);
      }
      ProjScalar.push_back(ProjScalarX);
    }
    this->SweepProjScalar.push_back(ProjScalar);
  }

  SweepDirFilter sf(&this->SweepDir, &this->SweepProjScalar);
}
