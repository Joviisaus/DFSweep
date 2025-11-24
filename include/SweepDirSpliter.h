#ifndef __SWEEP_DIR_SPLITER_H__
#define __SWEEP_DIR_SPLITER_H__

#include "CTMesh.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <float.h>
#include <iostream>
#include <queue>
#include <unordered_set>
#ifdef ENABLE_OMP
#include <omp.h>
#endif

class SweepDirSpliter {
public:
  SweepDirSpliter(MeshLib::CTMesh *mesh, std::vector<Eigen::Vector3f> *SweepDir,
                  std::vector<std::vector<std::vector<std::vector<float>>>>
                      *SweepProjScalar,
                  std::vector<std::vector<std::vector<int>>> FieldLabel);

protected:
  MeshLib::CTMesh *mesh;
  std::vector<Eigen::Vector3f> *SweepDir;
  std::vector<std::vector<std::vector<std::vector<float>>>> *SweepProjScalar;
  std::vector<std::vector<std::vector<int>>> FieldLabel;
  std::vector<int> LabelList;
  Eigen::MatrixXi LabelTopo;
  void splitDisconnectedGroups(std::vector<std::vector<int>> &EnergyLabel,
                               const Eigen::MatrixXi &LabelTopo);
  void LabelTopoGen();
  void SweepDirSplit();
  void SweepMask(std::vector<std::vector<int>> &EnergyLabel);
};

#endif
