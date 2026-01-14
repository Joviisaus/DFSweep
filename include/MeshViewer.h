#ifndef __MESHVIEWER_H_
#define __MESHVIEWER_H_

#include "CTMesh.h"
#include "glm/fwd.hpp"
#include "polyscope/curve_network.h"
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
      std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjEnergy,
      std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists,
      std::vector<Eigen::Vector3f> SweepDir,
      std::vector<std::vector<std::vector<bool>>> ForbiddenBoundaryPoints,
      std::vector<std::vector<std::vector<float>>> GradianceDiff,
      std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord);

protected:
  MeshLib::CTMesh *mesh;
  glm::vec3 bound_low{};
  glm::vec3 bound_high{};
  std::vector<Eigen::Vector3f> SweepDir;
  std::vector<int> label;
  uint32_t dimX;
  uint32_t dimY;
  uint32_t dimZ;
  float *scalarVals;
  float *GradianceDiff;
  float *ForbiddenBoundaryPoints;
  int *GradianceScalar;
  std::vector<float *> SweepProjScalars;
  std::vector<Eigen::Vector3f> sharpPoints;
  std::vector<std::array<size_t, 2>> Curves;
  std::vector<float *> SweepProjEnergies;
  std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists;
  std::vector<std::vector<float>> vertices;
  std::vector<Eigen::Vector3f> VertColors;
  std::vector<Eigen::Vector3f> FaceColors;
  std::vector<int> FaceSweepTypes;
  std::vector<std::vector<int>> faces;
};

#endif
