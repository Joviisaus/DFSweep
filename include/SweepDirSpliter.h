#ifndef __SWEEP_DIR_SPLITER_H__
#define __SWEEP_DIR_SPLITER_H__

#include "CTMesh.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <iostream>
#include <unordered_set>

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
};

#endif
