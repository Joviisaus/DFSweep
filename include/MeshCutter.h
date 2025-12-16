#ifndef __MESH_CUTTER_H__
#define __MESH_CUTTER_H__

#include "CTMesh.h"
#include "CuttingBox.h"

class MeshCutter {
public:
  MeshCutter(MeshLib::CTMesh *mesh,
             std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists);

private:
  MeshLib::CTMesh *mesh;
  CuttingBox *cb;
  std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists;

  void MeshCut();
};

#endif
