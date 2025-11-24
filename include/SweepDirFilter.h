#ifndef __SWEEP_DIR_FILTER_H
#define __SWEEP_DIR_FILTER_H

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <iostream>

inline float RotateZero = 0.1;

class SweepDirFilter {
public:
  SweepDirFilter::SweepDirFilter(
      std::vector<Eigen::Vector3f> *SweepDir,
      std::vector<std::vector<std::vector<std::vector<float>>>>
          *SweepProjScalar,
      std::vector<std::vector<std::vector<int>>> FieldLabel);

protected:
  int Filting();
  void SweepDirFilting();
  void SweepDirCleaning();
  int xSize;
  int ySize;
  int zSize;
  int RestSize;
  std::vector<Eigen::Vector3f> *SweepDir;
  std::vector<int> MarkedSweep;
  std::vector<std::vector<std::vector<int>>> FieldLabel;
  std::vector<std::vector<std::vector<std::vector<float>>>> *SweepProjScalar;
  std::vector<std::vector<std::vector<float>>> RestField;
};
#endif
