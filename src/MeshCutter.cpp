
#include "MeshCutter.h"
#include "Mesh/iterators.h"

MeshCutter::MeshCutter(
    MeshLib::CTMesh *mesh,
    std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists) {
  this->mesh = mesh;
  this->CuttingHexLists = CuttingHexLists;
}

void MeshCutter::MeshCut() {
  for (MeshLib::MeshFaceIterator mfiter(this->mesh); !mfiter.end(); ++mfiter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mfiter.value());
  }
}
