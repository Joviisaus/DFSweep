#include "MeshViewer.h"
#include "Mesh/iterators.h"
#include "polyscope/slice_plane.h"
#include "polyscope/types.h"
#include <cmath>
#include <limits>

int MeshViewer::setMesh(MeshLib::CTMesh *mesh) {
  this->vertices.clear();
  this->faces.clear();
  this->VertColors.clear();
  this->FaceColors.clear();
  this->VertBlockColors.clear();
  this->FaceBlockColors.clear();
  this->VertSweepBlock.clear();
  this->FaceSweepBlock.clear();
  this->label.clear();
  this->FaceSweepTypes.clear();
  this->FaceBottomHighlight.clear();

  auto blockColorAt = [&](int blockIdx) -> Eigen::Vector3f {
    if (blockIdx >= 0 &&
        blockIdx < static_cast<int>(this->blockColors.size())) {
      return this->blockColors[static_cast<size_t>(blockIdx)];
    }
    return Eigen::Vector3f(0.5f, 0.5f, 0.5f);
  };

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
    int blockIdx = v->cflabel();
    this->VertSweepBlock.push_back(blockIdx);
    Eigen::Vector3f vc = blockColorAt(blockIdx);
    this->VertColors.push_back(vc);
    this->VertBlockColors.push_back(vc);
    this->label.push_back(v->label());
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
    int blockIdx = f->sweeplabel();
    this->FaceSweepBlock.push_back(blockIdx);
    Eigen::Vector3f fc = blockColorAt(blockIdx);
    this->FaceColors.push_back(fc);
    this->FaceBlockColors.push_back(fc);
    int ft = f->sweepFaceType();
    this->FaceSweepTypes.push_back(ft);
    // 底面(5)=黄，内柱面源壁(6)=青，外柱面目标壁(7)=橙
    if (ft == 5) {
      this->FaceBottomHighlight.emplace_back(1.0f, 0.92f, 0.15f);
    } else if (ft == 6) {
      this->FaceBottomHighlight.emplace_back(0.15f, 0.95f, 1.0f);
    } else if (ft == 7) {
      this->FaceBottomHighlight.emplace_back(1.0f, 0.55f, 0.15f);
    } else {
      this->FaceBottomHighlight.push_back(fc * 0.35f);
    }
  }

  return 0;
}

void MeshViewer::setGrid(
    std::vector<std::vector<std::vector<float>>> Field,
    std::vector<std::vector<std::vector<int>>> GradianceCount,
    std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjScalar,
    std::vector<std::vector<std::vector<std::vector<float>>>> SweepProjEnergy,
    std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists,
    std::vector<Eigen::Vector3f> SweepDir,
    std::vector<std::vector<std::vector<bool>>> ForbiddenBoundaryPoints,
    std::vector<std::vector<std::vector<float>>> GradianceDiff,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>> Coord,
    const std::vector<bool> &hexIsNonPlanar,
    const std::vector<Eigen::Vector3f> &blockColors,
    const std::vector<int> &displayHexIndices,
    const std::vector<std::string> &sweepEnergyNames) {
  this->CuttingHexLists = CuttingHexLists;
  this->SweepDir = SweepDir;
  this->hexIsNonPlanar = hexIsNonPlanar;
  this->blockColors = blockColors;
  this->displayHexIndices = displayHexIndices;
  this->sweepEnergyNames = sweepEnergyNames;
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
  this->ForbiddenBoundaryPoints = new float[totalSize];
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
        bool forbidden = false;
        if (!ForbiddenBoundaryPoints.empty() &&
            x < ForbiddenBoundaryPoints.size() &&
            y < ForbiddenBoundaryPoints[x].size() &&
            z < ForbiddenBoundaryPoints[x][y].size()) {
          forbidden = ForbiddenBoundaryPoints[x][y][z];
        }
        this->ForbiddenBoundaryPoints[index] = forbidden ? 1.0f : 0.0f;
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

void MeshViewer::setSweepHexMeshes(const std::vector<SweepHexMesh> &meshes) {
  this->sweepHexMeshes = meshes;
}

int MeshViewer::show() {
  polyscope::init();
  auto mesh = polyscope::registerSurfaceMesh("Mesh", vertices, faces);
  mesh->addVertexColorQuantity("Sweep Block Color", this->VertColors)
      ->setEnabled(true);
  mesh->addVertexScalarQuantity("Sweep Block", this->VertSweepBlock);
  mesh->addVertexScalarQuantity("Label", this->label);
  mesh->addFaceColorQuantity("Sweep Block Color", this->FaceColors)
      ->setEnabled(true);
  mesh->addFaceScalarQuantity("Sweep Block", this->FaceSweepBlock);
  mesh->addFaceScalarQuantity("FaceSweepTypes", this->FaceSweepTypes);
  mesh->addFaceColorQuantity("Source Faces (surface)", this->FaceBottomHighlight)
      ->setEnabled(true);
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
  psGrid->addNodeScalarQuantity(
      "Valid Points", std::make_tuple(ForbiddenBoundaryPoints, nData));
  for (int i = 0; i < SweepProjScalars.size(); i++) {
    std::string str = "Accept for Sweep Direction " + std::to_string(i);
    psGrid->addNodeScalarQuantity(str,
                                  std::make_tuple(SweepProjScalars[i], nData));
  };
  for (int i = 0; i < SweepProjEnergies.size(); i++) {
    std::string str = (i < static_cast<int>(this->sweepEnergyNames.size()) &&
                       !this->sweepEnergyNames[static_cast<size_t>(i)].empty())
                          ? this->sweepEnergyNames[static_cast<size_t>(i)]
                          : "Sweep Energy " + std::to_string(i);
    auto *eq = psGrid
                   ->addNodeScalarQuantity(
                       str, std::make_tuple(SweepProjEnergies[i], nData))
                   ->setColorMap("coolwarm");
    if (i == 0) {
      eq->setEnabled(true);
    }
  };
  scalarQ->setEnabled(false);
  // --- 绘制六面体 (Cutting Hexahedra) ---
  // 假设 CuttingHexLists 是 std::vector<std::map<int, Eigen::Vector3f>> 类型
  // 每个 map 包含 8 个角点 (索引 0 到 7)

  auto blockColor = [&](int idx) -> glm::vec3 {
    if (idx >= 0 && idx < static_cast<int>(blockColors.size())) {
      const Eigen::Vector3f &c = blockColors[static_cast<size_t>(idx)];
      return glm::vec3(c.x(), c.y(), c.z());
    }
    if (idx >= 0 && idx < static_cast<int>(SweepDir.size())) {
      const Eigen::Vector3f &d = SweepDir[static_cast<size_t>(idx)];
      return glm::vec3(std::abs(d.x()), std::abs(d.y()), std::abs(d.z()));
    }
    return glm::vec3(0.8f, 0.8f, 0.8f);
  };

  std::vector<std::array<size_t, 2>> hex_edges = {
      {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3},
      {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7},
  };

  std::vector<int> hexIndices;
  if (!displayHexIndices.empty()) {
    hexIndices = displayHexIndices;
  } else {
    for (int i = 0; i < static_cast<int>(CuttingHexLists.size()); ++i) {
      hexIndices.push_back(i);
    }
  }

  int cylinderRegionId = 0;
  for (int hex_id : hexIndices) {
    if (hex_id < 0 || hex_id >= static_cast<int>(CuttingHexLists.size())) {
      continue;
    }
    const auto &hex_map = CuttingHexLists[static_cast<size_t>(hex_id)];
    if (hex_map.size() != 8) {
      std::cerr << "Warning: Hexahedron " << hex_id
                << " does not have 8 vertices. Skipping." << std::endl;
      continue;
    }

    std::vector<glm::vec3> hex_vertices;
    hex_vertices.reserve(8);
    for (int i = 0; i < 8; ++i) {
      const Eigen::Vector3f &eigen_v = hex_map.at(i);
      hex_vertices.push_back(glm::vec3(eigen_v.x(), eigen_v.y(), eigen_v.z()));
    }

    bool isNonPlanar =
        hex_id < static_cast<int>(hexIsNonPlanar.size()) &&
        hexIsNonPlanar[static_cast<size_t>(hex_id)];
    std::string hexName =
        isNonPlanar
            ? "Cylinder Sweep Region " + std::to_string(cylinderRegionId++)
            : "Vertical Sweep Region";

    float minEdge = std::numeric_limits<float>::max();
    for (size_t ei = 0; ei < hex_edges.size(); ++ei) {
      const glm::vec3 &a = hex_vertices[hex_edges[ei][0]];
      const glm::vec3 &b = hex_vertices[hex_edges[ei][1]];
      minEdge = std::min(minEdge, glm::length(b - a));
    }
    if (minEdge < 1e-5f) {
      std::cerr << "Warning: " << hexName << " has degenerate edges (min="
                << minEdge << ").\n";
    }

    auto cn = polyscope::registerCurveNetwork(hexName, hex_vertices, hex_edges);
    glm::vec3 wireColor = blockColor(hex_id);
    if (isNonPlanar) {
      wireColor = glm::vec3(1.0f, 0.95f, 0.2f);
    }
    cn->setColor(wireColor);
    cn->setEnabled(true);
    if (isNonPlanar) {
      cylinderRegionId++;
    }
  }

  const std::array<glm::vec3, 4> hexVolColors = {
      glm::vec3(0.2f, 0.85f, 1.0f), glm::vec3(1.0f, 0.55f, 0.2f),
      glm::vec3(0.55f, 1.0f, 0.35f), glm::vec3(0.95f, 0.45f, 0.95f)};
  for (size_t mi = 0; mi < this->sweepHexMeshes.size(); ++mi) {
    const SweepHexMesh &hm = this->sweepHexMeshes[mi];
    if ((hm.hexes.empty() && hm.wedges.empty()) || hm.nodes.empty()) {
      continue;
    }

    std::vector<glm::vec3> glmVerts;
    glmVerts.reserve(hm.nodes.size());
    for (const auto &v : hm.nodes) {
      glmVerts.emplace_back(v.x(), v.y(), v.z());
    }

    // Polyscope 仅支持 tet/hex：wedge 拆成 3 个 tet
    std::vector<std::array<uint32_t, 4>> tetCells;
    tetCells.reserve(hm.wedges.size() * 3);
    for (const auto &w : hm.wedges) {
      const uint32_t a = static_cast<uint32_t>(w[0]);
      const uint32_t b = static_cast<uint32_t>(w[1]);
      const uint32_t c = static_cast<uint32_t>(w[2]);
      const uint32_t d = static_cast<uint32_t>(w[3]);
      const uint32_t e = static_cast<uint32_t>(w[4]);
      const uint32_t f = static_cast<uint32_t>(w[5]);
      tetCells.push_back({a, b, c, d});
      tetCells.push_back({b, c, d, e});
      tetCells.push_back({c, d, e, f});
    }

    std::vector<std::array<uint32_t, 8>> hexCells;
    hexCells.reserve(hm.hexes.size());
    for (const auto &h : hm.hexes) {
      std::array<uint32_t, 8> cell{};
      for (int k = 0; k < 8; ++k) {
        cell[static_cast<size_t>(k)] =
            static_cast<uint32_t>(h[static_cast<size_t>(k)]);
      }
      hexCells.push_back(cell);
    }

    std::string volName = hm.name.empty()
                              ? ("Sweep Hex Volume " + std::to_string(mi))
                              : hm.name;
    polyscope::VolumeMesh *vol = nullptr;
    if (!hexCells.empty() && !tetCells.empty()) {
      vol = polyscope::registerTetHexMesh(volName, glmVerts, tetCells, hexCells);
    } else if (!hexCells.empty()) {
      vol = polyscope::registerHexMesh(volName, glmVerts, hexCells);
    } else {
      vol = polyscope::registerTetMesh(volName, glmVerts, tetCells);
    }
    if (vol) {
      vol->setColor(hexVolColors[mi % hexVolColors.size()]);
      vol->setEnabled(true);
    }

    if (!hm.imprintEdges.empty()) {
      std::vector<glm::vec3> impVerts;
      std::vector<std::array<size_t, 2>> impEdges;
      for (const auto &seg : hm.imprintEdges) {
        size_t base = impVerts.size();
        impVerts.emplace_back(seg[0].x(), seg[0].y(), seg[0].z());
        impVerts.emplace_back(seg[1].x(), seg[1].y(), seg[1].z());
        impEdges.push_back({base, base + 1});
      }
      std::string impName = "Imprint Features " + std::to_string(mi);
      auto impNet = polyscope::registerCurveNetwork(impName, impVerts, impEdges);
      impNet->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
      impNet->setEnabled(true);
    }
  }

  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

  polyscope::view::upDir = polyscope::UpDir::NegZUp;
  polyscope::show();
  return 0;
}
