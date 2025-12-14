#include "CuttingBox.h"

CuttingBox::CuttingBox(
    std::vector<Eigen::Vector3f> SweepDir,
    std::vector<std::vector<std::vector<std::vector<float>>>> *Energy,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord,
    std::vector<std::vector<std::vector<int>>> FieldLabel,
    std::vector<std::vector<std::vector<float>>> Field,
    std::vector<PrimeData> primes, int id) {
  this->SweepDir = SweepDir;
  this->Coord = Coord;
  this->Energy = Energy;
  this->FieldLabel = FieldLabel;
  this->Field = Field;
  this->primes = primes;
  this->id = id;
  STEP_SIZE = (this->Coord[0][0][0] - this->Coord[0][0][1]).norm();
  this->PrecomputeForbiddenPoints();
  this->DirCompute();
  this->PositionInit();
  this->TunePosition();
  this->EnergyUpdate();
}

void CuttingBox::PrecomputeForbiddenPoints() {
  int D1 = this->Field.size();
  int D2 = this->Field[0].size();
  int D3 = this->Field[0][0].size();

  ForbiddenBoundaryPoints.resize(
      D1, std::vector<std::vector<bool>>(D2, std::vector<bool>(D3, false)));

  const float FIELD_THRESHOLD = 2 * STEP_SIZE;

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (!this->primes[this->FieldLabel[x][y][z]].isPlane &&
            std::abs(this->Field[x][y][z]) < FIELD_THRESHOLD) {
          this->ForbiddenBoundaryPoints[x][y][z] = true;
        }
      }
    }
  }
}
// 新增函数实现 (在 .cpp 文件中)
bool CuttingBox::CheckForbiddenPointsInNewRegion(int bound_index,
                                                 float original_value,
                                                 float new_value) {
  // 确定检查的轴和边界值
  const Eigen::Vector3f *dir;
  if (bound_index <= 1)
    dir = &this->dirx;
  else if (bound_index <= 3)
    dir = &this->diry;
  else
    dir = &this->dirz;

  // 确定检查范围 (较小值到较大值)
  float min_check = std::min(original_value, new_value);
  float max_check = std::max(original_value, new_value);

  // 容忍度，用于判断点是否“在”边界上或“在”新区域内。
  // 由于移动是以 STEP_SIZE 为单位的，这里的容忍度需要比 STEP_SIZE 小很多。
  const float PROJECTION_TOLERANCE = STEP_SIZE * 1e-2f;

  int D1 = this->ForbiddenBoundaryPoints.size();
  int D2 = this->ForbiddenBoundaryPoints[0].size();
  int D3 = this->ForbiddenBoundaryPoints[0][0].size();

#ifdef ENABLE_OMP
// 这里的 OpenMP 并行需要注意，因为我们是 early exit (找到一个 forbidden point
// 就返回 false) 使用 #pragma omp critical
// 或原子操作会损失性能，但这里是非并行区域调用，故不加。
#endif
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (this->ForbiddenBoundaryPoints[x][y][z]) {
          const Eigen::Vector3f &P = this->Coord[x][y][z];
          float projection = P.dot(*dir);

          // 检查点是否位于移动后的新边界上，或者在旧边界和新边界之间的区域内
          // P·dir 应该在 [min_check, max_check] 之间 (考虑容忍度)
          if (projection >= min_check - PROJECTION_TOLERANCE &&
              projection <= max_check + PROJECTION_TOLERANCE) {

            // 如果是收缩操作，且点位于新旧边界之间，则禁止
            // 如果是扩张操作，且点位于新边界上，则禁止

            // 简化检查：只要禁止点位于新边界点 new_value 附近，或者位于
            // new_value 划定的新区域内，就禁止。 由于 OBB
            // 内部的点不影响边界有效性，我们只关心点是否在新边界 new_value
            // 附近。 最简单的策略是只检查点是否“落在了”新边界 new_value 附近
            if (std::abs(projection - new_value) < PROJECTION_TOLERANCE) {
              return false; // 找到一个禁止点
            }
          }
        }
      }
    }
  }

  return true; // 区域内没有禁止点
}
/**
 * @brief Computes MatchingEnergy for all sweep directions relative to the
 * region defined by this->id, and selects three non-coplanar directions:
 * 1. dirx is fixed to SweepDir[this->id].
 * 2. diry and dirz are selected from the remaining directions based on minimum
 * MatchingEnergy and non-coplanarity relative to the chosen axes.
 */
void CuttingBox::DirCompute() {
  this->dirx = Eigen::Vector3f::Zero();
  this->diry = Eigen::Vector3f::Zero();
  this->dirz = Eigen::Vector3f::Zero();

  if (this->id < 0 || this->id >= this->SweepDir.size()) {
    std::cerr << "Error: Invalid ID for SweepDir array." << std::endl;
    std::cout << "Single Sweep Volume;" << std::endl;
    return;
  }

  this->dirx = this->SweepDir[this->id];

  std::vector<std::pair<float, int>> energy_index_pairs;
  energy_index_pairs.reserve(this->SweepDir.size());

  for (int index = 0; index < this->SweepDir.size(); index++) {
    float ME = 0;
    if (index == this->id) {
      ME = FLT_MAX;
    } else {
      for (size_t x = 0; x < this->Energy[0][index].size(); x++) {
        for (size_t y = 0; y < this->Energy[0][index][0].size(); y++) {
          for (size_t z = 0; z < this->Energy[0][index][0][0].size(); z++) {
            if (this->Energy[0][id][x][y][z] < -2e-2) {
              ME += this->Energy[0][index][x][y][z];
            }
          }
        }
      }
    }

    if (index != this->id) {
      energy_index_pairs.push_back({ME, index});
    }
  }

  std::sort(energy_index_pairs.begin(), energy_index_pairs.end(),
            [](const std::pair<float, int> &a, const std::pair<float, int> &b) {
              return a.first < b.first;
            });

  std::vector<int> final_indices;
  final_indices.reserve(2);
  const float PLANE_TOLERANCE = 1e-4f;

  Eigen::Vector3f temp_diry = Eigen::Vector3f::Zero();

  for (const auto &pair : energy_index_pairs) {
    int current_index = pair.second;
    const Eigen::Vector3f &v_current = this->SweepDir[current_index];

    if (final_indices.empty()) {
      if ((this->dirx.cross(v_current)).norm() > PLANE_TOLERANCE) {
        final_indices.push_back(current_index);
        temp_diry = v_current;
      }
    } else if (final_indices.size() == 1) {
      float mixed_product =
          std::abs(this->dirx.dot(temp_diry.cross(v_current)));

      if (mixed_product > PLANE_TOLERANCE) {
        final_indices.push_back(current_index);
        this->diry = temp_diry;
        this->dirz = v_current;

        break;
      }
    }
  }

  if (final_indices.size() == 2) {

  } else if (final_indices.size() == 1) {
    this->diry = temp_diry; // 确认 diry

    Eigen::Vector3f cross_dir = this->dirx.cross(this->diry);

    if (cross_dir.norm() < PLANE_TOLERANCE) {
      final_indices.clear();
    } else {
      this->dirz = cross_dir;
    }
  }

  if (final_indices.size() < 1) {
    std::cout << "Single Sweep Volume;" << std::endl;

    this->dirx = this->SweepDir[this->id];
    this->diry = Eigen::Vector3f::Zero();
    this->dirz = Eigen::Vector3f::Zero();
    return;
  }

  this->dirx.normalize();
  this->diry.normalize();
  this->dirz.normalize();
}

void CuttingBox::TunePosition() {
  // Update box tune stratagy to enlarge the cutting box
  // New Strategy: In each iteration, find the single step change
  // (out of 12 possibilities) that yields the maximum energy reduction,
  // and apply only that change.

  float *bounds[] = {&this->MinX, &this->MaxX, &this->MinY,
                     &this->MaxY, &this->MinZ, &this->MaxZ};

  const char *names[] = {"MinX", "MaxX", "MinY", "MaxY", "MinZ", "MaxZ"};

  float current_energy = ComputeTotalEnergy();
  float previous_energy = std::numeric_limits<float>::max();
  int iteration = 0;

  std::cout << "Starting Energy: " << current_energy << std::endl;

  while (iteration < MAX_ITERATIONS &&
         std::abs(previous_energy - current_energy) > ENERGY_TOLERANCE) {

    previous_energy = current_energy;
    float best_energy_in_iteration = current_energy;
    float best_value = 0.0f;
    int best_param_index = -1;
    float *best_param_ptr = nullptr;
    bool improvement_found = false;

    for (int i = 0; i < 6; ++i) {
      float *current_param = bounds[i];
      float original_value = *current_param;

      float moves[] = {-STEP_SIZE, STEP_SIZE, -3 * STEP_SIZE, 3 * STEP_SIZE};

      for (float move : moves) {
        float new_value = original_value + move;
        *current_param = new_value;
        bool valid_move = true;
        if (i == 0 || i == 1) {
          if (this->MinX >= this->MaxX)
            valid_move = false;
        } else if (i == 2 || i == 3) {
          if (this->MinY >= this->MaxY)
            valid_move = false;
        } else if (i == 4 || i == 5) {
          if (this->MinZ >= this->MaxZ)
            valid_move = false;
        }
        if (valid_move && !ValidBox()) {
          valid_move = false;
        }
        if (valid_move) {
          if (!CheckForbiddenPointsInNewRegion(i, original_value, new_value)) {
            valid_move = false;
          }
        }
        if (valid_move) {
          float new_energy = ComputeTotalEnergy();

          if (new_energy < best_energy_in_iteration) {
            best_energy_in_iteration = new_energy;
            best_value = *current_param;
            best_param_index = i;
            best_param_ptr = current_param;
            improvement_found = true;
          }
        }

        *current_param = original_value;
      }
    }

    if (improvement_found) {
      *best_param_ptr = best_value;
      current_energy = best_energy_in_iteration;
      // std::cout << "Iteration " << iteration << ": Applied change to "
      //           << names[best_param_index] << ". New Energy: " <<
      //           current_energy
      //           << std::endl;
    } else {
      // std::cout << "Iteration " << iteration << ": No improving move found."
      //           << std::endl;
      break;
    }

    iteration++;
  }

  std::cout << "Tuned Energy: " << current_energy
            << " After iteration: " << iteration << std::endl;
}

float CuttingBox::ComputeTotalEnergy() {
  float total_energy = 0.0f;

  int D1 = this->Energy[0][this->id].size();
  int D2 = this->Energy[0][this->id][0].size();
  int D3 = this->Energy[0][this->id][0][0].size();

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        const Eigen::Vector3f &P = this->Coord[x][y][z];

        float projX = P.dot(this->dirx);
        float projY = P.dot(this->diry);
        float projZ = P.dot(this->dirz);

        bool is_inside = (projX >= this->MinX && projX <= this->MaxX) &&
                         (projY >= this->MinY && projY <= this->MaxY) &&
                         (projZ >= this->MinZ && projZ <= this->MaxZ);

        if (is_inside) {
          total_energy += this->Energy[0][this->id][x][y][z];
        }
      }
    }
  }

  return total_energy;
}

void CuttingBox::EnergyUpdate() {
  int D1 = this->Energy[0][this->id].size();
  int D2 = this->Energy[0][this->id][0].size();
  int D3 = this->Energy[0][this->id][0][0].size();

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        // 当前采样点的三维坐标
        const Eigen::Vector3f &P = this->Coord[x][y][z];

        float projX = P.dot(this->dirx);
        float projY = P.dot(this->diry);
        float projZ = P.dot(this->dirz);

        bool is_inside = (projX >= this->MinX && projX <= this->MaxX) &&
                         (projY >= this->MinY && projY <= this->MaxY) &&
                         (projZ >= this->MinZ && projZ <= this->MaxZ);

        if (is_inside) {
          for (int i = 0; i < this->SweepDir.size(); i++) {
            if ((this->SweepDir[i].cross(this->SweepDir[this->id])).norm() >
                1e-6)
              this->Energy[0][i][x][y][z] = 100;
          }
        }
      }
    }
  }
}
void CuttingBox::PositionInit() {
  this->MinX = std::numeric_limits<float>::max();
  this->MinY = std::numeric_limits<float>::max();
  this->MinZ = std::numeric_limits<float>::max();
  this->MaxX = std::numeric_limits<float>::lowest();
  this->MaxY = std::numeric_limits<float>::lowest();
  this->MaxZ = std::numeric_limits<float>::lowest();

  this->UpperBoundX = std::numeric_limits<float>::lowest();
  this->UpperBoundY = std::numeric_limits<float>::lowest();
  this->UpperBoundZ = std::numeric_limits<float>::lowest();
  this->DownBoundX = std::numeric_limits<float>::max();
  this->DownBoundY = std::numeric_limits<float>::max();
  this->DownBoundZ = std::numeric_limits<float>::max();
  int MinEnergyLabel = -1;
  float MinEnergy = std::numeric_limits<float>::max();
#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < this->Energy[0][id].size(); x++) {
    for (int y = 0; y < this->Energy[0][id][0].size(); y++) {
      for (int z = 0; z < this->Energy[0][id][0][0].size(); z++) {
        if (this->Energy[0][id][x][y][z] < MinEnergy) {
          MinEnergy = this->Energy[0][id][x][y][z];
          MinEnergyLabel = this->FieldLabel[x][y][z];
        }
      }
    }
  }

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < this->Energy[0][id].size(); x++) {
    for (int y = 0; y < this->Energy[0][id][0].size(); y++) {
      for (int z = 0; z < this->Energy[0][id][0][0].size(); z++) {
        if (x == 0 || y == 0 || z == 0 || x == this->Energy[0][id].size() - 1 ||
            y == this->Energy[0][id][0].size() - 1 ||
            z == this->Energy[0][id][0][0].size() - 1) {
          float dx = this->Coord[x][y][z].dot(dirx);
          float dy = this->Coord[x][y][z].dot(diry);
          float dz = this->Coord[x][y][z].dot(dirz);
          if (dx < this->DownBoundX)
            this->DownBoundX = dx;
          if (dx > this->UpperBoundX) {
            this->UpperBoundX = dx;
          }
          if (dy < this->DownBoundY)
            this->DownBoundY = dy;
          if (dy > this->UpperBoundY)
            this->UpperBoundY = dy;
          if (dz < this->DownBoundZ)
            this->DownBoundZ = dz;
          if (dz > this->UpperBoundZ)
            this->UpperBoundZ = dz;
        }
        if (this->Energy[0][id][x][y][z] > 99)
          continue;
        if (this->FieldLabel[x][y][z] != MinEnergyLabel ||
            abs(this->Field[x][y][z]) > 2 * STEP_SIZE)
          continue;
        float dx = this->Coord[x][y][z].dot(dirx);
        float dy = this->Coord[x][y][z].dot(diry);
        float dz = this->Coord[x][y][z].dot(dirz);
        if (dx < this->MinX)
          this->MinX = dx;
        if (dx > this->MaxX) {
          this->MaxX = dx;
        }
        if (dy < this->MinY)
          this->MinY = dy;
        if (dy > this->MaxY)
          this->MaxY = dy;
        if (dz < this->MinZ)
          this->MinZ = dz;
        if (dz > this->MaxZ)
          this->MaxZ = dz;
      }
    }
  }
}

/**
 * @brief Calculates the world coordinates of the eight vertices of the
 * bounding box.
 * * @return std::map<int, Eigen::Vector3f> A map where the key is the vertex
 * index (0-7) and the value is the vertex world coordinate.
 */
std::map<int, Eigen::Vector3f> CuttingBox::GetBoxVertices() {
  std::map<int, Eigen::Vector3f> vertices;

  // 存储所有的边界值
  float boundsX[] = {this->MinX, this->MaxX};
  float boundsY[] = {this->MinY, this->MaxY};
  float boundsZ[] = {this->MinZ, this->MaxZ};

  int vertex_index = 0;

  // 遍历所有 2*2*2 = 8 种边界组合
  for (int i = 0; i < 2; ++i) {     // i=0 corresponds to MinX, i=1 to MaxX
    for (int j = 0; j < 2; ++j) {   // j=0 corresponds to MinY, j=1 to MaxY
      for (int k = 0; k < 2; ++k) { // k=0 corresponds to MinZ, k=1 to MaxZ

        // 获取当前的边界值
        float boundX = boundsX[i];
        float boundY = boundsY[j];
        float boundZ = boundsZ[k];

        // 求解得到角点的世界坐标
        Eigen::Vector3f V = SolveVertex(boundX, boundY, boundZ);

        // 将结果存入 Map
        vertices[vertex_index] = V;
        vertex_index++;
      }
    }
  }

  return vertices;
}
/**
 * @brief Solves the linear system to find the world coordinate of a vertex
 * defined by three boundary projections.
 * * @param boundX The projection value on dirx (MinX or MaxX).
 * @param boundY The projection value on diry (MinY or MaxY).
 * @param boundZ The projection value on dirz (MinZ or MaxZ).
 * @return Eigen::Vector3f The world coordinate of the vertex.
 */
Eigen::Vector3f CuttingBox::SolveVertex(float boundX, float boundY,
                                        float boundZ) {
  Eigen::Matrix3f A;
  A.row(0) = this->dirx.transpose();
  A.row(1) = this->diry.transpose();
  A.row(2) = this->dirz.transpose();

  Eigen::Vector3f B(boundX, boundY, boundZ);

  Eigen::Vector3f V = A.colPivHouseholderQr().solve(B);

  return V;
}

bool CuttingBox::ValidBox() {
  if (this->Coord.empty() || this->Coord[0].empty() ||
      this->Coord[0][0].empty()) {
    return false;
  }

  const int x_max = this->Coord.size() - 1;
  const int y_max = this->Coord[0].size() - 1;
  const int z_max = this->Coord[0][0].size() - 1;

  const Eigen::Vector3f lower_bound_P = this->Coord[0][0][0];

  const Eigen::Vector3f upper_bound_P = this->Coord[x_max][y_max][z_max];

  std::map<int, Eigen::Vector3f> points = this->GetBoxVertices();

  for (int i = 0; i < 8; i++) {
    const Eigen::Vector3f &p = points[i];

    if (p[0] > upper_bound_P[0] || p[1] > upper_bound_P[1] ||
        p[2] > upper_bound_P[2]) {
      return false;
    }

    if (p[0] < lower_bound_P[0] || p[1] < lower_bound_P[1] ||
        p[2] < lower_bound_P[2]) {
      return false;
    }
  }

  return true;
}
