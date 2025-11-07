#ifndef __MESHVIEWER_H_
#define __MESHVIEWER_H_

#include "CTMesh.h"
#include "glm/fwd.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_grid.h"
#include <Eigen/Eigen>

class MeshViewer {
public:
  int show();
  int setMesh(MeshLib::CTMesh *mesh);
  void setGrid(
      std::vector<std::vector<std::vector<float>>> Field,
      std::vector<std::vector<std::vector<int>>> GradianceCount,
      std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjScalar,
      std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord);

protected:
  MeshLib::CTMesh *mesh;
  glm::vec3 bound_low{};
  glm::vec3 bound_high{};
  uint32_t dimX;
  uint32_t dimY;
  uint32_t dimZ;
  float *scalarVals;
  int *GradianceScalar;
  std::vector<float *> SweepProjScalars;
  std::vector<std::vector<float>> vertices;
  std::vector<std::vector<int>> faces;
};

#endif
