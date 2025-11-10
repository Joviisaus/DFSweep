#include "MeshViewer.h"
#include "Mesh/iterators.h"
#include "polyscope/types.h"

int MeshViewer::setMesh(MeshLib::CTMesh *mesh) {
  this->vertices.clear();
  this->faces.clear();
  int id = 1;
  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); mviter++) {
    MeshLib::CToolVertex *v =
        static_cast<MeshLib::CToolVertex *>(mviter.value());
    v->id() = id;
    id++;
    std::vector<float> point;
    point.resize(3);
    point[0] = v->point()[0];
    point[1] = v->point()[1];
    point[2] = v->point()[2];
    this->vertices.push_back(point);
  }

  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); mfiter++) {
    MeshLib::CTMesh::CFace *f =
        static_cast<MeshLib::CToolFace *>(mfiter.value());
    std::vector<int> vid;
    vid.clear();
    for (MeshLib::CTMesh::FaceVertexIterator fviter(f); !fviter.end();
         fviter++) {
      vid.push_back(fviter.value()->id() - 1);
    }
    this->faces.push_back(vid);
  }

  return 0;
}

void MeshViewer::setGrid(
    std::vector<std::vector<std::vector<float>>> Field,
    std::vector<std::vector<std::vector<int>>> GradianceCount,
    std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjScalar,
    std::vector<std::vector<std::vector<float>>> GradianceDiff,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord) {
  this->bound_low = {Coord.front().front().front()[0],
                     Coord.front().front().front()[1],
                     Coord.front().front().front()[2]};
  this->bound_high = {
      Coord.back().back().back()[0],
      Coord.back().back().back()[1],
      Coord.back().back().back()[2],
  };
  this->dimX = Coord.size();
  this->dimY = Coord.front().size();
  this->dimZ = Coord.front().front().size();
  SweepProjScalars.resize(SweepProjScalar.size());
  size_t totalSize = this->dimX * this->dimY * this->dimZ;
  this->scalarVals = new float[totalSize];
  this->GradianceScalar = new int[totalSize];
  this->GradianceDiff = new float[totalSize];
  for (int i = 0; i < SweepProjScalars.size(); i++) {
    SweepProjScalars[i] = new float[totalSize];
  }
  size_t index = 0;
  for (size_t z = 0; z < this->dimZ; ++z) {
    for (size_t y = 0; y < this->dimY; ++y) {
      for (size_t x = 0; x < this->dimX; ++x) {
        this->scalarVals[index] = Field[x][y][z];
        this->GradianceScalar[index] = GradianceCount[x][y][z];
        this->GradianceDiff[index] = GradianceDiff[x][y][z];
        for (int i = 0; i < SweepProjScalars.size(); i++) {
          SweepProjScalars[i][index] = SweepProjScalar[i][x][y][z];
        }
        index++;
      }
    }
  }
}

int MeshViewer::show() {
  polyscope::init();
  polyscope::registerSurfaceMesh("Mesh", vertices, faces);
  polyscope::VolumeGrid *psGrid = polyscope::registerVolumeGrid(
      "Field", {dimX, dimY, dimZ}, bound_low, bound_high);
  uint32_t nData = dimX * dimY * dimZ;
  polyscope::VolumeGridNodeScalarQuantity *scalarQ =
      psGrid->addNodeScalarQuantity("Distance Field",
                                    std::make_tuple(scalarVals, nData));
  psGrid->addNodeScalarQuantity("Gradiance Count",
                                std::make_tuple(GradianceScalar, nData));
  psGrid->addNodeScalarQuantity("Gradiance Diff",
                                std::make_tuple(GradianceDiff, nData));
  for (int i = 0; i < SweepProjScalars.size(); i++) {
    std::string str = "Accept for Sweep Direction " + std::to_string(i);
    psGrid->addNodeScalarQuantity(str,
                                  std::make_tuple(SweepProjScalars[i], nData));
  };
  scalarQ->setEnabled(true);
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  polyscope::view::upDir = polyscope::UpDir::NegZUp;
  polyscope::show();
  return 0;
}
