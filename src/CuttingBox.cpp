#include "CuttingBox.h"

CuttingBox::CuttingBox(
    std::vector<Eigen::Vector3f> SweepDir,
    std::vector<std::vector<std::vector<std::vector<float>>>> Energy,
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
}

/**
 * @brief Computes MatchingEnergy for all sweep directions relative to the
 * region defined by this->id, and selects three non-coplanar directions:
 * 1. dirx is fixed to SweepDir[this->id].
 * 2. diry and dirz are selected from the remaining directions based on minimum
 * MatchingEnergy and non-coplanarity relative to the chosen axes.
 */
void CuttingBox::DirCompute() {
  // 1. 初始化方向向量
  this->dirx = Eigen::Vector3f::Zero();
  this->diry = Eigen::Vector3f::Zero();
  this->dirz = Eigen::Vector3f::Zero();

  // 检查 this->id 是否在范围内
  if (this->id < 0 || this->id >= this->SweepDir.size()) {
    std::cerr << "Error: Invalid ID for SweepDir array." << std::endl;
    std::cout << "Single Sweep Volume;" << std::endl;
    return;
  }

  // 强制设置 dirx 为 this->SweepDir[this->id]
  this->dirx = this->SweepDir[this->id];

  // 2. 计算 MatchingEnergy (排除 this->id 自身)
  // 存储 (Energy, Original Index) 对，用于排序
  std::vector<std::pair<float, int>> energy_index_pairs;
  energy_index_pairs.reserve(this->SweepDir.size());

  for (int index = 0; index < this->SweepDir.size(); index++) {
    float ME = 0;
    if (index == this->id) {
      // 将自身的能量设为最大值，确保它不会被选中作为 diry 或 dirz
      ME = FLT_MAX;
    } else {
      // MatchingEnergy calculation logic:
      // Sum Energy[index] where Energy[id] is less than a negative threshold.
      for (size_t x = 0; x < this->Energy[index].size(); x++) {
        for (size_t y = 0; y < this->Energy[index][0].size(); y++) {
          for (size_t z = 0; z < this->Energy[index][0][0].size(); z++) {
            if (this->Energy[id][x][y][z] < -2e-2) {
              ME += this->Energy[index][x][y][z];
            }
          }
        }
      }
    }

    // 只需要存储用于选择 diry, dirz 的方向
    if (index != this->id) {
      energy_index_pairs.push_back({ME, index});
    }
  }

  // 3. 排序 (按 MatchingEnergy 升序排列，即找最小能量)
  std::sort(energy_index_pairs.begin(), energy_index_pairs.end(),
            [](const std::pair<float, int> &a, const std::pair<float, int> &b) {
              return a.first < b.first;
            });

  // 4. 选择 diry 和 dirz
  std::vector<int> final_indices; // 存储 diry 和 dirz 的索引
  final_indices.reserve(2);
  const float PLANE_TOLERANCE = 1e-4f;

  // 临时存储 diry 的索引和向量
  Eigen::Vector3f temp_diry = Eigen::Vector3f::Zero();

  for (const auto &pair : energy_index_pairs) {
    int current_index = pair.second;
    const Eigen::Vector3f &v_current = this->SweepDir[current_index];

    if (final_indices.empty()) {
      // 寻找 diry: 必须与 dirx 不平行
      // 检查非平行性：叉积的模长 > 容差
      if ((this->dirx.cross(v_current)).norm() > PLANE_TOLERANCE) {
        final_indices.push_back(current_index);
        temp_diry = v_current;
      }
    } else if (final_indices.size() == 1) {
      // 寻找 dirz: 必须与 dirx 和 diry 不共面

      // 检查非共面性：混合积的绝对值 > 容差
      // 混合积: dirx . (temp_diry x v_current)
      float mixed_product =
          std::abs(this->dirx.dot(temp_diry.cross(v_current)));

      if (mixed_product > PLANE_TOLERANCE) {
        final_indices.push_back(current_index);
        this->diry = temp_diry; // 确认 diry
        this->dirz = v_current; // 确认 dirz

        // 找到了两个，停止搜索
        break;
      }
    }
  }

  // 5. Post-processing and Final Assignment

  if (final_indices.size() == 2) {
    // Case 1: 找到了 diry 和 dirz
    // dirx, diry, dirz 已在上一步赋值

  } else if (final_indices.size() == 1) {
    // Case 2: 只找到了 diry (temp_diry), dirz 需要合成
    this->diry = temp_diry; // 确认 diry

    // 计算叉积合成 dirz
    Eigen::Vector3f cross_dir = this->dirx.cross(this->diry);

    // 检查叉积是否有效（理论上非平行向量的叉积不会是零向量）
    if (cross_dir.norm() < PLANE_TOLERANCE) {
      // 理论上不会发生，除非 dirx 和 diry 叉积失败，
      // 此时回退到 Single Sweep 逻辑
      final_indices.clear(); // 强制进入 Single Sweep 状态
    } else {
      this->dirz = cross_dir;
    }
  }

  if (final_indices.size() < 1) {
    // Case 3: 只找到了 dirx，或一个方向都没找到（极少情况）
    std::cout << "Single Sweep Volume;" << std::endl;

    // 清除 diry, dirz，仅保留 dirx (由步骤 1 设置)
    this->dirx = this->SweepDir[this->id]; // 确保 dirx 仍然是固定的初始方向
    this->diry = Eigen::Vector3f::Zero();
    this->dirz = Eigen::Vector3f::Zero();
    return;
  }

  // 6. 最终归一化: 确保 dirx, diry, dirz 都是单位向量
  this->dirx.normalize();
  this->diry.normalize();
  this->dirz.normalize();
}

void CuttingBox::TunePosition() {

  // 存储六个边界的指针，方便在循环中访问
  float *bounds[] = {&this->MinX, &this->MaxX, &this->MinY,
                     &this->MaxY, &this->MinZ, &this->MaxZ};

  // 存储六个边界的名称 (仅用于调试/理解)
  const char *names[] = {"MinX", "MaxX", "MinY", "MaxY", "MinZ", "MaxZ"};

  // 初始能量
  float current_energy = ComputeTotalEnergy();
  float previous_energy = std::numeric_limits<float>::max();
  int iteration = 0;

  std::cout << "Starting Energy:" << current_energy;

  while (iteration < MAX_ITERATIONS &&
         std::abs(previous_energy - current_energy) > ENERGY_TOLERANCE) {

    previous_energy = current_energy;
    bool improvement_found = false;

    // 遍历所有六个边界参数 (MinX, MaxX, MinY, MaxY, MinZ, MaxZ)
    for (int i = 0; i < 6; ++i) {
      float *current_param = bounds[i];
      float original_value = *current_param;
      float best_energy_in_step = current_energy;
      float best_value = original_value;

      // 尝试微调的方向：-STEP_SIZE (收缩) 和 +STEP_SIZE (扩张)
      // 注意：Min bounds 应该向负无穷方向移动，Max bounds
      // 应该向正无穷方向移动， 但为简化，我们尝试两个方向。更严谨的做法是：
      // 对于 MinX, MinY, MinZ: 尝试 -STEP_SIZE (收缩) 和 +STEP_SIZE (扩张)
      // 对于 MaxX, MaxY, MaxZ: 尝试 +STEP_SIZE (收缩) 和 -STEP_SIZE (扩张)
      // 以下使用一个更通用的尝试方式：
      float moves[] = {-STEP_SIZE, STEP_SIZE};

      for (float move : moves) {
        // 1. 尝试移动参数
        *current_param = original_value + move;

        // 2. 边界约束检查 (例如：MinX 必须小于 MaxX)
        // 这里的约束检查必须在外部处理，例如：
        // 如果是 MinX (i=0) 被修改，需要确保 MinX < MaxX
        // 如果是 MaxX (i=1) 被修改，需要确保 MinX < MaxX

        // 简单起见，我们只检查 X, Y, Z 对
        bool valid_move = true;
        if (i == 0 || i == 1) { // MinX, MaxX
          if (this->MinX >= this->MaxX)
            valid_move = false;
        } else if (i == 2 || i == 3) { // MinY, MaxY
          if (this->MinY >= this->MaxY)
            valid_move = false;
        } else if (i == 4 || i == 5) { // MinZ, MaxZ
          if (this->MinZ >= this->MaxZ)
            valid_move = false;
        }

        if (!valid_move) {
          // 恢复原始值并跳过这次尝试
          *current_param = original_value;
          continue;
        }

        // 3. 计算新能量
        float new_energy = ComputeTotalEnergy();

        // 4. 判断是否改善
        if (new_energy < best_energy_in_step) {
          best_energy_in_step = new_energy;
          best_value = *current_param;
          improvement_found = true;
        }

        // 5. 恢复原始值，为下一个 move 尝试做准备
        *current_param = original_value;
      }

      // 6. 更新参数到本步骤中找到的最佳值
      if (improvement_found && best_value != original_value) {
        *current_param = best_value;
        current_energy = best_energy_in_step;
      } else {
        // 确保参数回到原始值
        *current_param = original_value;
      }
    } // End of 6 parameters loop

    iteration++;

    // 如果一轮下来没有任何改善，则停止
    if (current_energy >= previous_energy && iteration > 0) {
      break;
    }

  } // End of while loop
  //
  std::cout << " Tuned Energy: " << current_energy << std::endl;
}
float CuttingBox::ComputeTotalEnergy() {
  float total_energy = 0.0f;

  // 遍历所有采样点
  size_t D1 = this->Energy[this->id].size();
  size_t D2 = this->Energy[this->id][0].size();
  size_t D3 = this->Energy[this->id][0][0].size();

  for (size_t x = 0; x < D1; ++x) {
    for (size_t y = 0; y < D2; ++y) {
      for (size_t z = 0; z < D3; ++z) {
        // 当前采样点的三维坐标
        const Eigen::Vector3f &P = this->Coord[x][y][z];

        // 将点 P 投影到 dirx, diry, dirz 轴上
        // 注意：由于 dirx, diry, dirz 可能是非正交的，
        // 严格来说，这里应该是计算在非正交基下的坐标，
        // 但根据你的 MinX 定义方式，我们暂且使用投影 (Dot Product)。
        float projX = P.dot(this->dirx);
        float projY = P.dot(this->diry);
        float projZ = P.dot(this->dirz);

        // 判断点 P 是否在六面体内部
        bool is_inside = (projX >= this->MinX && projX <= this->MaxX) &&
                         (projY >= this->MinY && projY <= this->MaxY) &&
                         (projZ >= this->MinZ && projZ <= this->MaxZ);

        if (is_inside) {
          total_energy += this->Energy[this->id][x][y][z];
        }
      }
    }
  }

  return total_energy;
}
void CuttingBox::PositionInit() {
  this->MinX = 99999;
  this->MinY = 99999;
  this->MinZ = 99999;
  this->MaxX = -99999;
  this->MaxY = -99999;
  this->MaxZ = -99999;
  int MinEnergyLabel = -1;
  float MinEnergy = FLT_MAX;
  for (size_t x = 0; x < this->Energy[id].size(); x++) {
    for (size_t y = 0; y < this->Energy[id][0].size(); y++) {
      for (size_t z = 0; z < this->Energy[id][0][0].size(); z++) {
        if (this->Energy[id][x][y][z] < MinEnergy) {
          MinEnergy = this->Energy[id][x][y][z];
          MinEnergyLabel = this->FieldLabel[x][y][z];
        }
      }
    }
  }
  for (size_t x = 0; x < this->Energy[id].size(); x++) {
    for (size_t y = 0; y < this->Energy[id][0].size(); y++) {
      for (size_t z = 0; z < this->Energy[id][0][0].size(); z++) {
        if (this->Energy[id][x][y][z] > -1e-1)
          continue;
        // if (this->FieldLabel[x][y][z] != MinEnergyLabel ||
        //     this->Field[x][y][z] < 0)
        //   continue;
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
 * @brief Calculates the world coordinates of the eight vertices of the bounding
 * box.
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
