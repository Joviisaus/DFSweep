#ifndef __CYLINDER_CUTTING_BOX_H__
#define __CYLINDER_CUTTING_BOX_H__

#include "PrimeData.h"
#include <Eigen/Eigen>
#include <map>
#include <set>
#include <vector>

/**
 * @brief 柱面扫掠的“柱面切割盒”：能量 CuttingBox 在柱坐标下的对应物。
 *
 * 与平面 CuttingBox 用 (MinX..MaxZ) 六个笛卡尔边界不同，柱面扫掠以
 * “同一中轴、不同半径”为切割指标，区域由两组边界刻画：
 *   - 半径范围 [MinR, MaxR]：内壁半径 -> 外壁半径（扫掠的源面/目标面）
 *   - 轴向范围 [MinAx, MaxAx]：沿柱轴的坐标范围
 *
 * 沿用原能量 CuttingBox 的思路：从柱面 prime 处播种，按“环形截面 +
 * 径向扫掠能量”的约束沿轴向生长，直到截面不再是干净的环（遇到底座等特征）。
 */
class CylinderCuttingBox {
public:
  CylinderCuttingBox(
      const Eigen::Vector3f &axis, const Eigen::Vector3f &origin, float rInner,
      float rOuter,
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
      const std::vector<std::vector<std::vector<float>>> &Field,
      const std::vector<std::vector<std::vector<int>>> &FieldLabel,
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField,
      const std::set<int> &memberPrimeIds, float energyAngleThreshold = 0.3f);

  float GetMinR() const { return MinR; }
  float GetMaxR() const { return MaxR; }
  float GetMinAx() const { return MinAx; }
  float GetMaxAx() const { return MaxAx; }
  Eigen::Vector3f GetAxis() const { return axis; }
  Eigen::Vector3f GetOrigin() const { return origin; }

  /// 返回包住该环形管段的定向六面体（与平面 CuttingBox 的 8 顶点约定一致）。
  std::map<int, Eigen::Vector3f> GetBoxVertices() const;

  /// 全场径向扫掠能量（实体内部体素有值，外部为 0），供 VolumeGrid 可视化。
  std::vector<std::vector<std::vector<float>>> ComputeRadialEnergyField() const;

protected:
  Eigen::Vector3f axis;
  Eigen::Vector3f origin;
  Eigen::Vector3f crossY;
  Eigen::Vector3f crossZ;
  const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord;
  const std::vector<std::vector<std::vector<float>>> &Field;
  const std::vector<std::vector<std::vector<int>>> &FieldLabel;
  const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &GradField;
  std::set<int> memberPrimeIds;
  float threshold;
  float stepSize;
  int D1, D2, D3;

  float MinR, MaxR;   // 半径范围（内壁 -> 外壁）
  float MinAx, MaxAx; // 轴向范围

  float AxialOf(const Eigen::Vector3f &p) const;
  float RadialOf(const Eigen::Vector3f &p) const;
  /// 体素的径向扫掠能量（梯度与局部径向方向的契合度，越小越好）。
  float RadialEnergy(int x, int y, int z) const;
  void PositionInit();
  void TuneAxialByConstraint();
  /// 判断某个轴向薄层是否仍是“干净环形”（径向扫掠的有效截面）。
  bool SlabIsAnnular(float axLo, float axHi) const;
};

#endif
