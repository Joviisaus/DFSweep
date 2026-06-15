#include "CylinderCuttingBox.h"
#include "CuttingBox.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {
void BuildCrossFrame(const Eigen::Vector3f &axisDir, Eigen::Vector3f &crossY,
                     Eigen::Vector3f &crossZ) {
  Eigen::Vector3f arbitrary =
      (std::abs(axisDir.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axisDir.dot(arbitrary) * axisDir).normalized();
  crossZ = axisDir.cross(crossY).normalized();
}
} // namespace

CylinderCuttingBox::CylinderCuttingBox(
    const Eigen::Vector3f &axis_, const Eigen::Vector3f &origin_, float rInner,
    float rOuter, const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<std::vector<float>>> &Field,
    const std::vector<std::vector<std::vector<int>>> &FieldLabel,
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField,
    const std::set<int> &memberPrimeIds, float energyAngleThreshold)
    : axis(axis_.normalized()), origin(origin_), Coord(Coord), Field(Field),
      FieldLabel(FieldLabel), GradField(GradField),
      memberPrimeIds(memberPrimeIds), threshold(energyAngleThreshold) {
  BuildCrossFrame(axis, crossY, crossZ);
  D1 = static_cast<int>(Field.size());
  D2 = D1 > 0 ? static_cast<int>(Field[0].size()) : 0;
  D3 = D2 > 0 ? static_cast<int>(Field[0][0].size()) : 0;
  stepSize = (D1 > 0 && D2 > 0 && D3 > 1)
                 ? (Coord[0][0][0] - Coord[0][0][1]).norm()
                 : 0.05f;

  // 半径范围 = 内壁 -> 外壁（“同一中轴、不同半径”切割指标）
  const float rMargin = 0.5f * stepSize;
  MinR = std::max(0.0f, rInner - rMargin);
  MaxR = rOuter + rMargin;

  PositionInit();
  TuneAxialByConstraint();

  std::cout << "[CylinderCuttingBox] R=[" << MinR << ", " << MaxR << "] Ax=["
            << MinAx << ", " << MaxAx << "]\n";
}

float CylinderCuttingBox::AxialOf(const Eigen::Vector3f &p) const {
  return (p - origin).dot(axis);
}

float CylinderCuttingBox::RadialOf(const Eigen::Vector3f &p) const {
  Eigen::Vector3f rel = p - origin;
  return (rel - rel.dot(axis) * axis).norm();
}

float CylinderCuttingBox::RadialEnergy(int x, int y, int z) const {
  const Eigen::Vector3f &p = Coord[x][y][z];
  Eigen::Vector3f rel = p - origin;
  Eigen::Vector3f radial = rel - rel.dot(axis) * axis;
  if (radial.norm() < 1e-8f) {
    return static_cast<float>(M_PI) / 4.0f;
  }
  Eigen::Vector3f rdir = radial.normalized();
  const Eigen::Vector3f &g = GradField[x][y][z];
  if (g.norm() < 1e-8f) {
    return static_cast<float>(M_PI) / 4.0f;
  }
  float c = std::min(1.0f, std::abs(g.normalized().dot(rdir)));
  float angle = std::acos(c);
  // 径向扫掠：梯度与径向“平行”（管壁）或“垂直”（端盖）均契合
  return std::min(angle, std::abs(static_cast<float>(M_PI) / 2.0f - angle));
}

bool CylinderCuttingBox::SlabIsAnnular(float axLo, float axHi) const {
  const float rBand = 1.0f * stepSize;
  int annularSolid = 0;
  int outerSolid = 0; // 半径超出外壁的实体（说明是底座的方形截面）
  double energyAcc = 0.0;
  int energyCnt = 0;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] <= 0.0f) {
          continue; // 仅实体内部（本工程约定 Field>0 为实体）
        }
        const Eigen::Vector3f &p = Coord[x][y][z];
        float a = AxialOf(p);
        if (a < axLo || a >= axHi) {
          continue;
        }
        float r = RadialOf(p);
        if (r > MaxR + rBand) {
          outerSolid++;
        } else if (r >= MinR - rBand) {
          annularSolid++;
          energyAcc += RadialEnergy(x, y, z);
          energyCnt++;
        }
      }
    }
  }
  if (annularSolid == 0) {
    return false; // 该层没有环形实体 -> 模型在此结束
  }
  // 出现明显的环外实体 -> 进入底座的方形截面，停止生长
  if (outerSolid > annularSolid / 4 + 1) {
    return false;
  }
  // 能量约束：该层环形实体的径向扫掠契合度需足够好
  float meanEnergy =
      energyCnt > 0 ? static_cast<float>(energyAcc / energyCnt) : 1e9f;
  return meanEnergy < threshold * 1.5f;
}

void CylinderCuttingBox::PositionInit() {
  // 轴向种子 = 柱面 wall prime 的实体表面体素轴向范围
  MinAx = std::numeric_limits<float>::max();
  MaxAx = std::numeric_limits<float>::lowest();
  const float surfaceBand = 2.0f * stepSize;
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        int fl = FieldLabel[x][y][z];
        if (!memberPrimeIds.count(fl)) {
          continue;
        }
        if (std::abs(Field[x][y][z]) >= surfaceBand) {
          continue;
        }
        const Eigen::Vector3f &p = Coord[x][y][z];
        float r = RadialOf(p);
        if (r < MinR - stepSize || r > MaxR + stepSize) {
          continue; // 只用落在环壁半径附近的成员体素播种轴向
        }
        float a = AxialOf(p);
        MinAx = std::min(MinAx, a);
        MaxAx = std::max(MaxAx, a);
      }
    }
  }
  if (MinAx > MaxAx) {
    // 退化兜底：用全部实体环形体素
    for (int x = 0; x < D1; ++x) {
      for (int y = 0; y < D2; ++y) {
        for (int z = 0; z < D3; ++z) {
          if (Field[x][y][z] <= 0.0f) {
            continue;
          }
          const Eigen::Vector3f &p = Coord[x][y][z];
          float r = RadialOf(p);
          if (r < MinR - stepSize || r > MaxR + stepSize) {
            continue;
          }
          float a = AxialOf(p);
          MinAx = std::min(MinAx, a);
          MaxAx = std::max(MaxAx, a);
        }
      }
    }
  }
}

void CylinderCuttingBox::TuneAxialByConstraint() {
  if (MinAx > MaxAx) {
    return;
  }
  const float slab = stepSize;
  const int MAX_STEPS = 4 * std::max(D1, std::max(D2, D3)) + 10;

  // 向上（+axis）生长
  for (int i = 0; i < MAX_STEPS; ++i) {
    if (SlabIsAnnular(MaxAx, MaxAx + slab)) {
      MaxAx += slab;
    } else {
      break;
    }
  }
  // 向下（-axis）生长
  for (int i = 0; i < MAX_STEPS; ++i) {
    if (SlabIsAnnular(MinAx - slab, MinAx)) {
      MinAx -= slab;
    } else {
      break;
    }
  }

  const float endMargin = 0.5f * stepSize;
  MinAx -= endMargin;
  MaxAx += endMargin;
}

std::vector<std::vector<std::vector<float>>>
CylinderCuttingBox::ComputeRadialEnergyField() const {
  std::vector<std::vector<std::vector<float>>> grid(
      D1, std::vector<std::vector<float>>(
              D2, std::vector<float>(D3, EXTERIOR_SWEEP_ENERGY)));
  for (int x = 0; x < D1; ++x) {
    for (int y = 0; y < D2; ++y) {
      for (int z = 0; z < D3; ++z) {
        if (Field[x][y][z] < 0.0f) {
          continue;
        }
        grid[x][y][z] = RadialEnergy(x, y, z);
      }
    }
  }
  return grid;
}

std::map<int, Eigen::Vector3f> CylinderCuttingBox::GetBoxVertices() const {
  auto corner = [&](float a, float cy, float cz) {
    return origin + a * axis + cy * crossY + cz * crossZ;
  };
  float axVals[] = {MinAx, MaxAx};
  float yVals[] = {-MaxR, MaxR};
  float zVals[] = {-MaxR, MaxR};
  std::map<int, Eigen::Vector3f> vertices;
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
