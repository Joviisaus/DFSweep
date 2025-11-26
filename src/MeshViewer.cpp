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
    std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjEnergy,
    std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists,
    std::vector<std::vector<std::vector<float>>> GradianceDiff,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord) {
  this->CuttingHexLists = CuttingHexLists;
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
  this->SweepProjScalars.resize(SweepProjScalar.size());
  this->SweepProjEnergies.resize(SweepProjEnergy.size());
  size_t totalSize = this->dimX * this->dimY * this->dimZ;
  this->scalarVals = new float[totalSize];
  this->GradianceScalar = new int[totalSize];
  this->GradianceDiff = new float[totalSize];
  for (int i = 0; i < SweepProjScalars.size(); i++) {
    SweepProjScalars[i] = new float[totalSize];
  }
  for (int i = 0; i < SweepProjEnergy.size(); i++) {
    SweepProjEnergies[i] = new float[totalSize];
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

        for (int i = 0; i < SweepProjEnergies.size(); i++) {
          SweepProjEnergies[i][index] = SweepProjEnergy[i][x][y][z];
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
  for (int i = 0; i < SweepProjEnergies.size(); i++) {
    std::string str = "Sweep Energy " + std::to_string(i);
    psGrid
        ->addNodeScalarQuantity(str,
                                std::make_tuple(SweepProjEnergies[i], nData))
        ->setColorMap("coolwarm");
  };
  scalarQ->setEnabled(true);
  // --- 绘制六面体 (Cutting Hexahedra) ---
  // 假设 CuttingHexLists 是 std::vector<std::map<int, Eigen::Vector3f>> 类型
  // 每个 map 包含 8 个角点 (索引 0 到 7)

  int hex_id = 0;
  for (const auto &hex_map : CuttingHexLists) {

    // 1. 准备顶点坐标
    // 转换为 polyscope 期望的 std::vector<glm::vec3> 格式
    std::vector<glm::vec3> hex_vertices;
    // 确保 map 中有 8 个点
    if (hex_map.size() != 8) {
      std::cerr << "Warning: Hexahedron " << hex_id
                << " does not have 8 vertices. Skipping." << std::endl;
      hex_id++;
      continue;
    }

    // Polyscope 要求顶点坐标按索引顺序排列 (0 到 7)
    for (int i = 0; i < 8; ++i) {
      // 从 Eigen::Vector3f 转换到 glm::vec3
      const Eigen::Vector3f &eigen_v = hex_map.at(i);
      hex_vertices.push_back(glm::vec3(eigen_v.x(), eigen_v.y(), eigen_v.z()));
    }

    // 2. 定义六面体的 12 条边
    // 边的索引对 (连接 hex_vertices 中的索引)
    // 假设标准的六面体顶点索引顺序如下：
    /*
        4-------7 (z max)
       /|      /|
      5-------6 |
      | 0-----|-3 (z min)
      |/      |/
      1-------2
    */
    std::vector<std::array<size_t, 2>> hex_edges = {
        // Z-edges (MinZ <-> MaxZ): (X, Y 固定)
        {0, 1}, // MinX, MinY: (000 <-> 001)
        {2, 3}, // MinX, MaxY: (010 <-> 011)
        {4, 5}, // MaxX, MinY: (100 <-> 101)
        {6, 7}, // MaxX, MaxY: (110 <-> 111)

        // Y-edges (MinY <-> MaxY): (X, Z 固定)
        {0, 2}, // MinX, MinZ: (000 <-> 010)
        {1, 3}, // MinX, MaxZ: (001 <-> 011)
        {4, 6}, // MaxX, MinZ: (100 <-> 110)
        {5, 7}, // MaxX, MaxZ: (101 <-> 111)

        // X-edges (MinX <-> MaxX): (Y, Z 固定)
        {0, 4}, // MinY, MinZ: (000 <-> 100)
        {1, 5}, // MinY, MaxZ: (001 <-> 101)
        {2, 6}, // MaxY, MinZ: (010 <-> 110)
        {3, 7}  // MaxY, MaxZ: (011 <-> 111)
    };
    // 3. 注册为曲线网络
    std::string name = "Cutting Hex " + std::to_string(hex_id);
    polyscope::registerCurveNetwork(name, hex_vertices, hex_edges)
        ->setColor({0.8f, 0.1f, 0.1f}) // 设置颜色，例如红色
        ->setRadius(0.005);            // 设置线的粗细

    hex_id++;
  }
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  polyscope::view::upDir = polyscope::UpDir::NegZUp;
  polyscope::show();
  return 0;
}
