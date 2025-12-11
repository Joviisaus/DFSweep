#include "CuttingBox.h"

CuttingBox::CuttingBox(
    std::vector<Eigen::Vector3f> SweepDir,
    std::vector<std::vector<std::vector<std::vector<float>>>> *Energy,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord,
    std::vector<std::vector<std::vector<int>>> FieldLabel,

    std::vector<std::vector<std::vector<float>>> Field, int id) {
  this->SweepDir = SweepDir;
  this->Coord = Coord;
  this->Energy = Energy;
  this->FieldLabel = FieldLabel;
  this->Field = Field;
  this->id = id;
  STEP_SIZE = (this->Coord[0][0][0] - this->Coord[0][0][1]).norm();
  this->DirCompute();
  this->PositionInit();
  this->TunePosition();
  this->EnergyUpdate();
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
  // TODO:
  // Update box tune stratagy to enlarge the cutting box

  float *bounds[] = {&this->MinX, &this->MaxX, &this->MinY,
                     &this->MaxY, &this->MinZ, &this->MaxZ};

  const char *names[] = {"MinX", "MaxX", "MinY", "MaxY", "MinZ", "MaxZ"};

  float current_energy = ComputeTotalEnergy();
  float previous_energy = std::numeric_limits<float>::max();
  int iteration = 0;

  std::cout << "Starting Energy:" << current_energy;

  while (iteration < MAX_ITERATIONS &&
         std::abs(previous_energy - current_energy) > ENERGY_TOLERANCE) {

    previous_energy = current_energy;
    bool improvement_found = false;

    for (int i = 0; i < 6; ++i) {
      float *current_param = bounds[i];
      float original_value = *current_param;
      float best_energy_in_step = current_energy;
      float best_value = original_value;

      float moves[] = {-STEP_SIZE, STEP_SIZE};

      for (float move : moves) {
        *current_param = original_value + move;

        bool valid_move = true;
        if (i == 0 || i == 1) { // MinX, MaxX
          if (this->MinX >= this->MaxX || !ValidBox())
            valid_move = false;
        } else if (i == 2 || i == 3) { // MinY, MaxY
          if (this->MinY >= this->MaxY || !ValidBox())
            valid_move = false;
        } else if (i == 4 || i == 5) { // MinZ, MaxZ
          if (this->MinZ >= this->MaxZ || !ValidBox())
            valid_move = false;
        }

        if (!valid_move) {
          *current_param = original_value;
          continue;
        }

        float new_energy = ComputeTotalEnergy();

        if (new_energy < best_energy_in_step) {
          best_energy_in_step = new_energy;
          best_value = *current_param;
          improvement_found = true;
        }

        *current_param = original_value;
      }

      if (improvement_found && best_value != original_value) {
        *current_param = best_value;
        current_energy = best_energy_in_step;
      } else {
        *current_param = original_value;
      }
    }

    iteration++;

    if (current_energy >= previous_energy && iteration > 0) {
      break;
    }
  }
  std::cout << " Tuned Energy: " << current_energy
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
        // 当前采样点的三维坐标
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
        if (this->Energy[0][id][x][y][z] > -1e-1)
          continue;
        if (this->FieldLabel[x][y][z] != MinEnergyLabel ||
            this->Field[x][y][z] < 0)
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
  // 1. 构造系数矩阵 A (方向向量作为行向量)
  Eigen::Matrix3f A;
  A.row(0) = this->dirx.transpose();
  A.row(1) = this->diry.transpose();
  A.row(2) = this->dirz.transpose();

  // 2. 构造常数向量 B (边界值)
  Eigen::Vector3f B(boundX, boundY, boundZ);

  // 3. 求解线性方程组 A * V = B
  // 使用 Eigen 的 .solve() 方法，它会自动处理求逆
  // V = A.inverse() * B
  // 注意：如果 A 是奇异矩阵（方向共面），求解会失败，但我们已假设方向不共面。
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
