#include "SweepDirSpliter.h"

SweepDirSpliter::SweepDirSpliter(
    MeshLib::CTMesh *mesh, std::vector<Eigen::Vector3f> *SweepDir,
    std::vector<std::vector<std::vector<std::vector<float>>>> *SweepProjScalar,
    std::vector<std::vector<std::vector<int>>> FieldLabel) {
  this->mesh = mesh;
  this->SweepDir = SweepDir;
  this->SweepProjScalar = SweepProjScalar;
  this->FieldLabel = FieldLabel;
  std::unordered_set<int> labelSet;
  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); mviter++) {
    MeshLib::CToolVertex *vert =
        static_cast<MeshLib::CToolVertex *>(mviter.value());
    int label = vert->label();
    labelSet.insert(label);
  }
  this->LabelList.assign(labelSet.begin(), labelSet.end());
  this->LabelTopo.resize(LabelList.size(), LabelList.size());
  this->LabelTopo.setZero();
}
