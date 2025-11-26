#ifndef __CUTTING_BOX_H__
#define __CUTTING_BOX_H__

#include "SweepDirFilter.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <float.h>
inline float STEP_SIZE = 0.05f;       // 每一步微调的量
const int MAX_ITERATIONS = 50;        // 最大迭代次数
const float ENERGY_TOLERANCE = 1e-4f; // 能量变化小于此值则停止

class CuttingBox {
public:
  CuttingBox(std::vector<Eigen::Vector3f> SweepDir,
             std::vector<std::vector<std::vector<std::vector<float>>>> Energy,
             std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord,
             int id);
  std::map<int, Eigen::Vector3f> GetBoxVertices();

protected:
  int id;
  std::vector<Eigen::Vector3f> SweepDir;
  std::vector<std::vector<std::vector<std::vector<float>>>> Energy;
  std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord;
  Eigen::Vector3f dirx;
  Eigen::Vector3f diry;
  Eigen::Vector3f dirz;
  float MinX;
  float MinY;
  float MinZ;
  float MaxX;
  float MaxY;
  float MaxZ;

  void DirCompute();
  void TunePosition();
  void PositionInit();
  float ComputeTotalEnergy();
  Eigen::Vector3f SolveVertex(float boundX, float boundY, float boundZ);
};
#endif
