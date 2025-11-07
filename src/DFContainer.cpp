#include "DFContainer.h"
#include "CTMesh.h"
#include <Eigen/src/Core/Matrix.h>

DistanceField::DistanceField() { this->primes.clear(); };

DistanceField::DistanceField(MeshLib::CTMesh *mesh) {
  this->primes.clear();
  SetMesh(mesh);
}

void DistanceField::SetMesh(MeshLib::CTMesh *mesh) {
  this->mesh = mesh;
  this->PointList.clear();
  for (MeshLib::MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
    std::vector<float> pointcoord;
    pointcoord.clear();
    pointcoord.push_back(viter.value()->point()[0]);
    pointcoord.push_back(viter.value()->point()[1]);
    pointcoord.push_back(viter.value()->point()[2]);
    this->PointList.push_back(pointcoord);
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

  // Todo: Complete Sweep Dir Detection;
  Eigen::Vector3f Dir1;
  Dir1[0] = 1;
  Dir1[1] = 0;
  Dir1[2] = 0;

  Eigen::Vector3f Dir2;
  Dir2[0] = 0;
  Dir2[1] = 1;
  Dir2[2] = 0;

  Eigen::Vector3f Dir3;
  Dir3[0] = 0;
  Dir3[1] = 0;
  Dir3[2] = 1;

  this->SweepDir.push_back(Dir1);
  this->SweepDir.push_back(Dir2);
  this->SweepDir.push_back(Dir3);

  return;
}

void DistanceField::BuildOctree() {
  if (PointList.empty())
    return;

  // 计算包围盒
  Eigen::Vector3f minPoint =
      Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
  Eigen::Vector3f maxPoint =
      Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());

  for (const auto &point : PointList) {
    Eigen::Vector3f p(point[0], point[1], point[2]);
    minPoint = minPoint.cwiseMin(p);
    maxPoint = maxPoint.cwiseMax(p);
  }

  // 扩展一点避免边界问题
  Eigen::Vector3f center = (minPoint + maxPoint) * 0.5f;
  float halfSize = (maxPoint - minPoint).norm() * 0.5f + 0.1f;

  octreeRoot = std::make_shared<OctreeNode>(center, halfSize);

  // 收集所有点的索引
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

  // 如果点数少于阈值或达到最大深度，设为叶子节点
  if (pointIndices.size() <= maxPointsPerNode || depth >= maxDepth) {
    node->pointIndices = pointIndices;
    node->isLeaf = true;
    return;
  }

  // 否则细分节点
  SubdivideNode(node);
  node->isLeaf = false;

  // 将点分配到子节点中
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

  // 递归构建子节点
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

  // 计算点到节点包围盒的距离
  Eigen::Vector3f diff = (point - node->center).cwiseAbs();
  float distToNode =
      (diff - Eigen::Vector3f::Constant(node->halfSize)).cwiseMax(0.0f).norm();

  // 如果是叶子节点，添加所有点
  if (node->isLeaf) {
    candidateIndices.insert(candidateIndices.end(), node->pointIndices.begin(),
                            node->pointIndices.end());
    return;
  }

  // 递归搜索子节点
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

float DistanceField::DistanceToMesh(const Eigen::Vector3f &point) {
  if (!mesh || PointList.empty())
    return 0.0f;

  // 使用八叉树找到候选点
  std::vector<int> candidateIndices;
  FindNearestPointsInOctree(point, octreeRoot, candidateIndices);

  if (candidateIndices.empty())
    return 0.0f;
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
        minDistance = this->DisCompute(point, v->label());
        break;
      }
    }
  }

  return minDistance;
}

double DistanceField::DisCompute(Eigen::Vector3f point, int label) {
  auto &m_params = this->primes[label].params;

  // std::cout << "Params: " << m_params[0] << " " << m_params[1] << " "
  //           << m_params[2] << " " << m_params[3] << " " << m_params[4] << " "
  //           << m_params[5] << " " << m_params[6] << " " << m_params[7] << " "
  //           << m_params[8] << " " << m_params[9] << std::endl;

  // 更安全的平面检测
  bool isPlane = (m_params[4] == 0 && m_params[5] == 0 && m_params[6] == 0 &&
                  m_params[7] == 0 && m_params[8] == 0 && m_params[9] == 0);

  if (isPlane) {
    double a = m_params[1], b = m_params[2], c = m_params[3], d = m_params[0];
    double norm = sqrt(a * a + b * b + c * c);

    // 检查法向量是否为零
    if (norm < 1e-10) {
      return std::numeric_limits<double>::infinity();
    }

    return std::abs(a * point[0] + b * point[1] + c * point[2] + d) / norm;
  }

  // 二次曲面情况
  int max_iter = 20;
  double lambda = 0.0;
  Eigen::Vector3f q = point;

  for (int i = 0; i < max_iter; ++i) {
    double x = q(0), y = q(1), z = q(2);

    // 计算函数值和梯度
    double F = m_params[0] + m_params[1] * x + m_params[2] * y +
               m_params[3] * z + m_params[4] * x * y + m_params[5] * x * z +
               m_params[6] * y * z + m_params[7] * x * x + m_params[8] * y * y +
               m_params[9] * z * z;

    Eigen::Vector3f gradF;
    gradF << m_params[1] + m_params[4] * y + m_params[5] * z +
                 2 * m_params[7] * x,
        m_params[2] + m_params[4] * x + m_params[6] * z + 2 * m_params[8] * y,
        m_params[3] + m_params[5] * x + m_params[6] * y + 2 * m_params[9] * z;

    // 构造残差
    Eigen::Vector4f residual;
    residual.head<3>() = q - point + lambda * gradF;
    residual(3) = F;

    // 构造雅可比矩阵
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

  BuildOctree();

  int xSize = Field.size();
  int ySize = Field[0].size();
  int zSize = Field[0][0].size();

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif // ENABLE_OMP
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        const Eigen::Vector3f &samplePoint = Coord[i][j][k];
        float distance = DistanceToMesh(samplePoint);
        Field[i][j][k] = distance;
      }
    }
  }

  GradianceCount.resize(xSize);
  for (int i = 0; i < xSize; ++i) {
    GradianceCount[i].resize(ySize);
    for (int j = 0; j < ySize; ++j) {
      GradianceCount[i][j].resize(zSize, 0);
    }
  }
#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif // ENABLE_OMP
  for (int i = 0; i < xSize; ++i) {
    for (int j = 0; j < ySize; ++j) {
      for (int k = 0; k < zSize; ++k) {
        if (i <= 0 || i >= xSize - 1 || j <= 0 || j >= ySize - 1 || k <= 0 ||
            k >= zSize - 1) {
          GradianceCount[i][j][k] = 1;
          continue;
        }

        float dx = 1.0f;
        float dy = 1.0f;
        float dz = 1.0f;

        float dxx =
            (Field[i + 1][j][k] - 2 * Field[i][j][k] + Field[i - 1][j][k]) /
            (dx * dx);
        float dyy =
            (Field[i][j + 1][k] - 2 * Field[i][j][k] + Field[i][j - 1][k]) /
            (dy * dy);
        float dzz =
            (Field[i][j][k + 1] - 2 * Field[i][j][k] + Field[i][j][k - 1]) /
            (dz * dz);

        float dxy = (Field[i + 1][j + 1][k] - Field[i + 1][j - 1][k] -
                     Field[i - 1][j + 1][k] + Field[i - 1][j - 1][k]) /
                    (4 * dx * dy);
        float dxz = (Field[i + 1][j][k + 1] - Field[i + 1][j][k - 1] -
                     Field[i - 1][j][k + 1] + Field[i - 1][j][k - 1]) /
                    (4 * dx * dz);
        float dyz = (Field[i][j + 1][k + 1] - Field[i][j + 1][k - 1] -
                     Field[i][j - 1][k + 1] + Field[i][j - 1][k - 1]) /
                    (4 * dy * dz);

        Eigen::Matrix3f hessian;
        hessian << dxx, dxy, dxz, dxy, dyy, dyz, dxz, dyz, dzz;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(hessian);
        if (eigensolver.info() != Eigen::Success) {
          GradianceCount[i][j][k] = 1;
          continue;
        }

        Eigen::Vector3f eigenvalues = eigensolver.eigenvalues();

        eigenvalues.normalize();

        int nonZeroCount = 0;
        for (int idx = 0; idx < 3; ++idx) {
          if (std::abs(eigenvalues[idx]) > epsilon) {
            nonZeroCount++;
          }
        }

        if (nonZeroCount == 0) nonZeroCount = 1;

          GradianceCount[i][j][k] = nonZeroCount;
      }
    }
  }
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

      primes.push_back(current_prime);
    }
  }

  file.close();

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

  // 写入维度信息
  file.write(reinterpret_cast<const char *>(&xSize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&ySize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&zSize), sizeof(int));

  // 写入数据
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

  // 写入维度信息
  file.write(reinterpret_cast<const char *>(&xSize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&ySize), sizeof(int));
  file.write(reinterpret_cast<const char *>(&zSize), sizeof(int));

  // 写入数据
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
