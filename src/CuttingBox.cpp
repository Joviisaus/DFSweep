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
  // this->TunePosition();
  this->TuneBoxBoundaryByConstraint();
  // this->EnergyUpdate();
}

void CuttingBox::PrecomputeForbiddenPoints() {
  int D1 = this->Field.size();
  int D2 = this->Field[0].size();
  int D3 = this->Field[0][0].size();

  ForbiddenBoundaryPoints.resize(
      D1, std::vector<std::vector<bool>>(D2, std::vector<bool>(D3, false)));

  const float FIELD_THRESHOLD = STEP_SIZE;

#ifdef ENABLE_OMP
#pragma omp parallel for collapse(3)
#endif
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (this->FieldLabel[x][y][z] < -1)
          continue;
        if (!(this->primes[this->FieldLabel[x][y][z]].isPlane) &&
            std::abs(this->Field[x][y][z]) < FIELD_THRESHOLD) {
          this->ForbiddenBoundaryPoints[x][y][z] = true;
        }
      }
    }
  }
}
bool CuttingBox::CheckForbiddenPointsInNewRegion(int bound_index,
                                                 float original_value,
                                                 float new_value) {
  const float PROJECTION_TOLERANCE = STEP_SIZE * 1e-2f;
  const Eigen::Vector3f *dir;
  int axis = bound_index / 2; // 0=X, 1=Y, 2=Z
  if (axis == 0)
    dir = &this->dirx;
  else if (axis == 1)
    dir = &this->diry;
  else
    dir = &this->dirz;
  int D1 = this->ForbiddenBoundaryPoints.size();
  int D2 = this->ForbiddenBoundaryPoints[0].size();
  int D3 = this->ForbiddenBoundaryPoints[0][0].size();
  bool ValidGrow = false;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (this->ForbiddenBoundaryPoints[x][y][z]) {
          const Eigen::Vector3f &P = this->Coord[x][y][z];
          float projection = P.dot(*dir);
          if (std::abs(projection - new_value) < PROJECTION_TOLERANCE) {
            bool is_on_new_boundary = false;
            float py = P.dot(this->diry);
            float pz = P.dot(this->dirz);
            if (axis == 0) { // 移动 X 边界
              float p_y = P.dot(this->diry);
              float p_z = P.dot(this->dirz);
              if (p_y >= this->MinY - PROJECTION_TOLERANCE &&
                  p_y <= this->MaxY + PROJECTION_TOLERANCE &&
                  p_z >= this->MinZ - PROJECTION_TOLERANCE &&
                  p_z <= this->MaxZ + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            } else if (axis == 1) { // 移动 Y 边界
              float p_x = P.dot(this->dirx);
              float p_z = P.dot(this->dirz);
              if (p_x >= this->MinX - PROJECTION_TOLERANCE &&
                  p_x <= this->MaxX + PROJECTION_TOLERANCE &&
                  p_z >= this->MinZ - PROJECTION_TOLERANCE &&
                  p_z <= this->MaxZ + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            } else if (axis == 2) { // 移动 Z 边界
              float p_x = P.dot(this->dirx);
              float p_y = P.dot(this->diry);
              if (p_x >= this->MinX - PROJECTION_TOLERANCE &&
                  p_x <= this->MaxX + PROJECTION_TOLERANCE &&
                  p_y >= this->MinY - PROJECTION_TOLERANCE &&
                  p_y <= this->MaxY + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            }
            if (is_on_new_boundary) {
              return false;
            }
          }
        }
        if (!ValidGrow && this->Field[x][y][z] >= -4 * STEP_SIZE) {
          const Eigen::Vector3f &P = this->Coord[x][y][z];
          float projection = P.dot(*dir);
          if (std::abs(projection - new_value) < PROJECTION_TOLERANCE) {
            bool is_on_new_boundary = false;
            float py = P.dot(this->diry);
            float pz = P.dot(this->dirz);
            if (axis == 0) { // 移动 X 边界
              float p_y = P.dot(this->diry);
              float p_z = P.dot(this->dirz);
              if (p_y >= this->MinY - PROJECTION_TOLERANCE &&
                  p_y <= this->MaxY + PROJECTION_TOLERANCE &&
                  p_z >= this->MinZ - PROJECTION_TOLERANCE &&
                  p_z <= this->MaxZ + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            } else if (axis == 1) { // 移动 Y 边界
              float p_x = P.dot(this->dirx);
              float p_z = P.dot(this->dirz);
              if (p_x >= this->MinX - PROJECTION_TOLERANCE &&
                  p_x <= this->MaxX + PROJECTION_TOLERANCE &&
                  p_z >= this->MinZ - PROJECTION_TOLERANCE &&
                  p_z <= this->MaxZ + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            } else if (axis == 2) { // 移动 Z 边界
              float p_x = P.dot(this->dirx);
              float p_y = P.dot(this->diry);
              if (p_x >= this->MinX - PROJECTION_TOLERANCE &&
                  p_x <= this->MaxX + PROJECTION_TOLERANCE &&
                  p_y >= this->MinY - PROJECTION_TOLERANCE &&
                  p_y <= this->MaxY + PROJECTION_TOLERANCE) {
                is_on_new_boundary = true;
              }
            }
            if (is_on_new_boundary) {
              ValidGrow = true;
            }
          }
        }
      }
    }
  }

  return ValidGrow;
}

void CuttingBox::DirCompute() {
    // 1. 初始化并校验基础方向 dirx
    this->dirx = Eigen::Vector3f::Zero();
    this->diry = Eigen::Vector3f::Zero();
    this->dirz = Eigen::Vector3f::Zero();

    if (this->id < 0 || this->id >= (int)this->SweepDir.size()) {
        std::cerr << "Error: Invalid ID for SweepDir array." << std::endl;
        return;
    }

    // dirx 保持不变
    this->dirx = this->SweepDir[this->id].normalized();

    // 2. 寻找与 dirx 最“垂直”的向量作为 diry
    float min_dot = 1.0f; // 点积越接近 0 越垂直
    int best_y_index = -1;
    const float COLINEAR_TOLERANCE = 0.999f; // 用于检查是否共线

    for (int i = 0; i < (int)this->SweepDir.size(); ++i) {
        if (i == this->id) continue;

        Eigen::Vector3f v_current = this->SweepDir[i].normalized();
        // 计算点积的绝对值，越接近 0 说明夹角越接近 90 度
        float dot_product = std::abs(this->dirx.dot(v_current));

        // 必须确保不与 dirx 共线（夹角不能太小或接近180度）
        if (dot_product < COLINEAR_TOLERANCE) {
            if (dot_product < min_dot) {
                min_dot = dot_product;
                best_y_index = i;
            }
        }
    }

    // 3. 处理查找结果
    if (best_y_index != -1) {
        Eigen::Vector3f v_y = this->SweepDir[best_y_index].normalized();

        // --- 施密特正交化确保严格垂直 ---
        // diry = v_y - (dirx · v_y) * dirx
        this->diry = (v_y - this->dirx.dot(v_y) * this->dirx).normalized();
        
        // dirz 由前两个方向叉乘得到
        this->dirz = this->dirx.cross(this->diry).normalized();

    } else {
        // 如果没有找到任何不共线的向量
        std::cout << "Single Sweep Volume: No suitable perpendicular vector found." << std::endl;
        this->diry = Eigen::Vector3f::Zero();
        this->dirz = Eigen::Vector3f::Zero();
    }
}

void CuttingBox::TuneBoxBoundaryByConstraint() {

  float *bounds[] = {&this->MinX, &this->MaxX, &this->MinY,
                     &this->MaxY, &this->MinZ, &this->MaxZ};

  auto search_for_valid_move = [&](int bound_index, int dir_multiplier,
                                   float max_search_range) -> bool {
    float *current_param = bounds[bound_index];
    float old_value = *current_param;

    float search_moves[] = {
        dir_multiplier * 5.0f * STEP_SIZE, dir_multiplier * 3.0f * STEP_SIZE,
        dir_multiplier * 1.0f * STEP_SIZE, dir_multiplier * 0.5f * STEP_SIZE};

    for (float move : search_moves) {
      float new_value = old_value + move;

      *current_param = new_value;

      if (this->MinX >= this->MaxX || this->MinY >= this->MaxY ||
          this->MinZ >= this->MaxZ) {
        *current_param = old_value; // 恢复
        continue;
      }

      if (!ValidBox()) {
        *current_param = old_value; // 恢复
        continue;
      }

      if (!CheckForbiddenPointsInNewRegion(bound_index, old_value, new_value)) {
        *current_param = old_value; // 恢复
        continue;
      }

      return true;
    }

    *current_param = old_value;
    return false;
  };

  for (int i = 0; i < 6; ++i) {

    int expansion_dir = (i % 2 != 0) ? 1 : -1;

    if (search_for_valid_move(i, expansion_dir,
                              /*max_search_range=*/100 * STEP_SIZE)) {
      continue; // 找到有效的扩张，继续下一个边界
    }

    int contraction_dir = -expansion_dir;

    if (search_for_valid_move(i, contraction_dir,
                              /*max_search_range=*/100 * STEP_SIZE)) {
      continue; // 找到有效的收缩，继续下一个边界
    }
  }

  // Tuning 完成，边界 MinX/MaxX 等被更新为找到的第一个满足约束的有效位置
}
void CuttingBox::TunePosition() {
  // Update box tune stratagy: In each iteration, sequentially test all 6
  // boundaries (i). For each boundary, apply the FIRST valid move that reduces
  // energy, and then move immediately to the next boundary (i+1).

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
    bool improvement_found = false;

    for (int i = 0; i < 6; ++i) {
      float *current_param = bounds[i];

      float moves[] = {-STEP_SIZE, STEP_SIZE, -3 * STEP_SIZE, 3 * STEP_SIZE};
      bool accepted_move_for_i = false;

      for (float move : moves) {

        float value_before_test = *current_param;

        float new_value = value_before_test + move;
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
          if (!CheckForbiddenPointsInNewRegion(i, value_before_test,
                                               new_value)) {
            valid_move = false;
          }
        }
        if (valid_move) {
          float new_energy = ComputeTotalEnergy();

          if (new_energy < current_energy) {
            current_energy = new_energy;
            improvement_found = true;
            accepted_move_for_i = true;

            break;
          } else {
            *current_param = value_before_test;
          }
        } else {
          *current_param = value_before_test;
        }
      }
    }

    if (!improvement_found) {
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
#pragma omp parallel for collapse(3) reduction(+ : total_energy)
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
        if (this->Energy[0][id][x][y][z] < MinEnergy &&
            this->FieldLabel[x][y][z] >= 0) {
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
