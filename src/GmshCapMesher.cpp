#include "GmshCapMesher.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace {
constexpr float kPi = 3.14159265358979323846f;

Eigen::Vector3f FromCapUV(float u, float v, float ax,
                          const SweepCapFrame &frame) {
  return frame.origin + ax * frame.axis + u * frame.crossY + v * frame.crossZ;
}

int PointId(int i, int j, int nx) { return 1 + j * (nx + 1) + i; }

void BuildCrossFrame(const Eigen::Vector3f &axisDir, Eigen::Vector3f &crossY,
                     Eigen::Vector3f &crossZ) {
  Eigen::Vector3f axis = axisDir;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f arbitrary =
      (std::abs(axis.dot(Eigen::Vector3f::UnitX())) < 0.9f)
          ? Eigen::Vector3f::UnitX()
          : Eigen::Vector3f::UnitY();
  crossY = (arbitrary - axis.dot(arbitrary) * axis).normalized();
  crossZ = axis.cross(crossY);
  if (crossZ.norm() < 1e-8f) {
    crossZ = Eigen::Vector3f::UnitZ();
  } else {
    crossZ.normalize();
  }
}

int LineIdH(int i, int j, int nx) {
  return 1 + j * nx + i;
}
int LineIdV(int i, int j, int nx, int ny) {
  return 1 + (ny + 1) * nx + i * ny + j;
}

bool RunGmsh(const std::string &geoPath, const std::string &mshPath,
             int dimension, const std::string &gmshExe) {
  std::ostringstream cmd;
  cmd << '"' << gmshExe << "\" \"" << geoPath << "\" -" << dimension
      << " -format msh2 -o \"" << mshPath << "\"";
  std::cout << "[GmshCapMesher] " << cmd.str() << "\n";
  int ret = std::system(cmd.str().c_str());
  if (!std::filesystem::exists(mshPath)) {
    std::cerr << "[GmshCapMesher] gmsh failed exit=" << ret << "\n";
    return false;
  }
  if (ret != 0) {
    std::cerr << "[GmshCapMesher] gmsh exit=" << ret
              << " (continuing with existing msh)\n";
  }
  return true;
}

bool ParseMsh2Quads(const std::string &mshPath,
                    std::vector<Eigen::Vector3f> &nodes,
                    std::vector<std::array<int, 4>> &quads) {
  std::ifstream in(mshPath);
  if (!in.is_open()) {
    return false;
  }
  nodes.clear();
  quads.clear();
  std::map<int, int> tagToIndex;
  std::string line;
  while (std::getline(in, line)) {
    if (line == "$Nodes") {
      int n = 0;
      in >> n;
      for (int i = 0; i < n; ++i) {
        int tag = 0;
        double x = 0, y = 0, z = 0;
        in >> tag >> x >> y >> z;
        tagToIndex[tag] = static_cast<int>(nodes.size());
        nodes.emplace_back(static_cast<float>(x), static_cast<float>(y),
                           static_cast<float>(z));
      }
      continue;
    }
    if (line == "$Elements") {
      int n = 0;
      in >> n;
      for (int i = 0; i < n; ++i) {
        int tag = 0, type = 0, ntags = 0;
        in >> tag >> type >> ntags;
        for (int t = 0; t < ntags; ++t) {
          int dummy = 0;
          in >> dummy;
        }
        if (type == 3) {
          std::array<int, 4> q{};
          for (int k = 0; k < 4; ++k) {
            int nodeTag = 0;
            in >> nodeTag;
            auto it = tagToIndex.find(nodeTag);
            if (it == tagToIndex.end()) {
              return false;
            }
            q[static_cast<size_t>(k)] = it->second;
          }
          quads.push_back(q);
        } else if (type == 10) {
          std::array<int, 4> q{};
          const int cornerTags[4] = {0, 1, 2, 3};
          for (int k : cornerTags) {
            int nodeTag = 0;
            in >> nodeTag;
            auto it = tagToIndex.find(nodeTag);
            if (it == tagToIndex.end()) {
              return false;
            }
            q[static_cast<size_t>(k)] = it->second;
          }
          for (int k = 4; k < 9; ++k) {
            int dummy = 0;
            in >> dummy;
          }
          quads.push_back(q);
        } else {
          int nNodes = (type == 1) ? 2 : (type == 2 ? 3 : 4);
          if (type == 15) {
            nNodes = 4;
          }
          for (int k = 0; k < nNodes; ++k) {
            int dummy = 0;
            in >> dummy;
          }
        }
      }
      break;
    }
  }
  return !nodes.empty() && !quads.empty();
}

bool WriteTranslationalTransfiniteGeo(const ImprintedCapMesh &cap,
                                      const SweepCapFrame &frame,
                                      float holeRadius,
                                      const std::string &geoPath) {
  const auto &u = cap.uSplits;
  const auto &v = cap.vSplits;
  if (u.size() < 2 || v.size() < 2) {
    return false;
  }
  const int nx = static_cast<int>(u.size()) - 1;
  const int ny = static_cast<int>(v.size()) - 1;
  const float ax = frame.axMin;
  const float lc =
      std::max(1e-3f, 0.25f * std::min(u.back() - u.front(), v.back() - v.front()) /
                            static_cast<float>(std::max(nx, ny)));

  std::ofstream geo(geoPath);
  if (!geo.is_open()) {
    return false;
  }

  geo << "Mesh.RecombineAll = 1;\n";
  geo << "Mesh.Smoothing = 0;\n";

  geo << "Mesh.Algorithm = 1;\n";
  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i <= nx; ++i) {
      int pid = PointId(i, j, nx);
      Eigen::Vector3f p = FromCapUV(u[static_cast<size_t>(i)],
                                    v[static_cast<size_t>(j)], ax, frame);
      geo << "Point(" << pid << ") = {" << p.x() << ", " << p.y() << ", "
          << p.z() << ", " << lc << "};\n";
    }
  }

  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int lid = LineIdH(i, j, nx);
      geo << "Line(" << lid << ") = {" << PointId(i, j, nx) << ", "
          << PointId(i + 1, j, nx) << "};\n";
    }
  }
  const int vLineBase = (ny + 1) * nx + 1;
  for (int i = 0; i <= nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      int lid = vLineBase + i * ny + j;
      geo << "Line(" << lid << ") = {" << PointId(i, j, nx) << ", "
          << PointId(i, j + 1, nx) << "};\n";
    }
  }

  geo << "Line Loop(1) = {";
  for (int i = 0; i < nx; ++i) {
    if (i > 0) {
      geo << ", ";
    }
    geo << LineIdH(i, 0, nx);
  }
  for (int j = 0; j < ny; ++j) {
    geo << ", " << (vLineBase + nx * ny + j);
  }
  for (int i = nx - 1; i >= 0; --i) {
    geo << ", " << -(LineIdH(i, ny, nx));
  }
  for (int j = ny - 1; j >= 0; --j) {
    geo << ", " << -(vLineBase + j);
  }
  geo << "};\n";
  geo << "Plane Surface(1) = {1};\n";

  geo << "Transfinite Curve {";
  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (i > 0 || j > 0) {
        geo << ", ";
      }
      geo << LineIdH(i, j, nx);
    }
  }
  geo << "} = 2;\n";
  geo << "Transfinite Curve {";
  for (int i = 0; i <= nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (i > 0 || j > 0) {
        geo << ", ";
      }
      geo << (vLineBase + i * ny + j);
    }
  }
  geo << "} = 2;\n";
  geo << "Transfinite Surface {1} = {" << PointId(0, 0, nx) << ", "
      << PointId(nx, 0, nx) << ", " << PointId(nx, ny, nx) << ", "
      << PointId(0, ny, nx) << "};\n";
  geo << "Recombine Surface {1};\n";
  geo << "Mesh 2;\n";
  return true;
}

bool WriteCylindricalInnerGeo(const SweepBlockRegion &block, int nTheta,
                              int nAxial, const std::string &geoPath) {
  Eigen::Vector3f axis = block.sweepAxis;
  if (axis.norm() < 1e-8f) {
    axis = Eigen::Vector3f::UnitY();
  } else {
    axis.normalize();
  }
  Eigen::Vector3f origin = block.sweepOrigin;
  if (!origin.allFinite()) {
    return false;
  }
  Eigen::Vector3f crossY, crossZ;
  BuildCrossFrame(axis, crossY, crossZ);

  const float r = std::max(1e-4f, block.radialInner);
  const float a0 = block.axialLower;
  const float a1 = block.axialUpper;
  if (!std::isfinite(r) || !std::isfinite(a0) || !std::isfinite(a1) ||
      a1 <= a0 + 1e-6f) {
    std::cerr << "[GmshCapMesher] invalid cylindrical params r=" << r
              << " ax=[" << a0 << "," << a1 << "] origin=("
              << origin.x() << "," << origin.y() << "," << origin.z() << ")\n";
    return false;
  }
  const float lc = std::max(1e-3f, (a1 - a0) / std::max(1, nAxial));

  std::ofstream geo(geoPath);
  if (!geo.is_open()) {
    return false;
  }

  geo << "Mesh.RecombineAll = 1;\n";
  geo << "Mesh.Smoothing = 0;\n";
  geo << "Mesh.Algorithm = 1;\n";

  const Eigen::Vector3f originV = origin;
  const Eigen::Vector3f axisV = axis;
  const Eigen::Vector3f crossYV = crossY;
  const Eigen::Vector3f crossZV = crossZ;
  const float rV = r;
  const float a0V = a0;
  const float a1V = a1;

  auto pos = [=](int it, int ia) {
    float theta = 2.0f * kPi * static_cast<float>(it) / static_cast<float>(nTheta);
    float axPos =
        a0V + static_cast<float>(ia) / static_cast<float>(nAxial) * (a1V - a0V);
    float ct = std::cos(theta);
    float st = std::sin(theta);
    float rx = ct * crossYV.x();
    float ry = ct * crossYV.y();
    float rz = ct * crossYV.z();
    if (std::abs(st) > 1e-8f) {
      rx += st * crossZV.x();
      ry += st * crossZV.y();
      rz += st * crossZV.z();
    }
    return Eigen::Vector3f(originV.x() + axPos * axisV.x() + rV * rx,
                           originV.y() + axPos * axisV.y() + rV * ry,
                           originV.z() + axPos * axisV.z() + rV * rz);
  };

  for (int ia = 0; ia <= nAxial; ++ia) {
    for (int it = 0; it <= nTheta; ++it) {
      int pid = 1 + ia * (nTheta + 1) + it;
      Eigen::Vector3f p = pos(it, ia);
      geo << "Point(" << pid << ") = {" << p.x() << ", " << p.y() << ", "
          << p.z() << ", " << lc << "};\n";
    }
  }

  int lineId = 1;
  for (int ia = 0; ia <= nAxial; ++ia) {
    for (int it = 0; it < nTheta; ++it) {
      int p0 = 1 + ia * (nTheta + 1) + it;
      int p1 = 1 + ia * (nTheta + 1) + it + 1;
      geo << "Line(" << lineId << ") = {" << p0 << ", " << p1 << "};\n";
      lineId++;
    }
  }
  const int hCount = (nAxial + 1) * nTheta;
  const int vBase = hCount + 1;
  for (int it = 0; it <= nTheta; ++it) {
    for (int ia = 0; ia < nAxial; ++ia) {
      int p0 = 1 + ia * (nTheta + 1) + it;
      int p1 = 1 + (ia + 1) * (nTheta + 1) + it;
      geo << "Line(" << lineId << ") = {" << p0 << ", " << p1 << "};\n";
      lineId++;
    }
  }

  geo << "Line Loop(1) = {";
  for (int it = 0; it < nTheta; ++it) {
    if (it > 0) {
      geo << ", ";
    }
    geo << (1 + it);
  }
  for (int ia = 0; ia < nAxial; ++ia) {
    geo << ", " << (vBase + nTheta * nAxial + ia);
  }
  for (int it = nTheta - 1; it >= 0; --it) {
    geo << ", " << -(1 + nAxial * nTheta + it);
  }
  for (int ia = nAxial - 1; ia >= 0; --ia) {
    geo << ", " << -(vBase + ia);
  }
  geo << "};\n";

  geo << "Surface(1) = {1};\n";

  geo << "Transfinite Curve {";
  for (int i = 1; i <= hCount; ++i) {
    if (i > 1) {
      geo << ", ";
    }
    geo << i;
  }
  geo << "} = 2;\n";
  geo << "Transfinite Curve {";
  for (int i = vBase; i < vBase + (nTheta + 1) * nAxial; ++i) {
    if (i > vBase) {
      geo << ", ";
    }
    geo << i;
  }
  geo << "} = 2;\n";

  const int c00 = 1;
  const int c10 = 1 + nTheta;
  const int c11 = 1 + nAxial * (nTheta + 1) + nTheta;
  const int c01 = 1 + nAxial * (nTheta + 1);
  geo << "Transfinite Surface {1} = {" << c00 << ", " << c10 << ", " << c11
      << ", " << c01 << "};\n";
  geo << "Recombine Surface {1};\n";
  geo << "Mesh 2;\n";
  return true;
}
} // namespace

std::string GmshCapMesher::FindExecutable() {
#ifdef GMSH_EXECUTABLE
  if (std::filesystem::exists(GMSH_EXECUTABLE)) {
    return GMSH_EXECUTABLE;
  }
#endif
  const char *env = std::getenv("GMSH_EXECUTABLE");
  if (env && std::filesystem::exists(env)) {
    return env;
  }
  const char *candidates[] = {"/opt/homebrew/bin/gmsh", "/usr/local/bin/gmsh",
                              "gmsh"};
  for (const char *c : candidates) {
    if (std::filesystem::exists(c)) {
      return c;
    }
  }
  return "";
}

bool GmshCapMesher::IsAvailable() { return !FindExecutable().empty(); }

bool GmshCapMesher::MeshTranslationalCap(ImprintedCapMesh &cap,
                                         const SweepCapFrame &frame,
                                         float holeRadius,
                                         const std::string &workDir) {
  if (holeRadius > 1e-4f) {
    // 带孔截面：压印分割已足够，用本地 structured quads 更贴合特征
    return false;
  }
  const std::string gmshExe = FindExecutable();
  if (gmshExe.empty()) {
    std::cerr << "[GmshCapMesher] gmsh not found. Install: brew install gmsh\n";
    return false;
  }
  std::filesystem::create_directories(workDir);
  const std::string geoPath = workDir + "/cap_trans.geo";
  const std::string mshPath = workDir + "/cap_trans.msh";
  if (!WriteTranslationalTransfiniteGeo(cap, frame, holeRadius, geoPath)) {
    return false;
  }
  if (!RunGmsh(geoPath, mshPath, 2, gmshExe)) {
    return false;
  }
  std::vector<Eigen::Vector3f> nodes;
  std::vector<std::array<int, 4>> quads;
  if (!ParseMsh2Quads(mshPath, nodes, quads) || quads.size() < 4) {
    std::cerr << "[GmshCapMesher] failed to parse quad mesh from " << mshPath
              << "\n";
    return false;
  }
  cap.capNodes = std::move(nodes);
  cap.quads = std::move(quads);
  std::cout << "[GmshCapMesher] translational cap: " << cap.quads.size()
            << " quads, " << cap.capNodes.size() << " nodes\n";
  return true;
}

bool GmshCapMesher::MeshCylindricalInnerCap(const SweepBlockRegion &block,
                                            int nTheta, int nAxial,
                                            ImprintedCapMesh &cap,
                                            const std::string &workDir) {
  const std::string gmshExe = FindExecutable();
  if (gmshExe.empty()) {
    return false;
  }
  std::filesystem::create_directories(workDir);
  const std::string geoPath = workDir + "/cap_cyl.geo";
  const std::string mshPath = workDir + "/cap_cyl.msh";
  if (!WriteCylindricalInnerGeo(block, nTheta, nAxial, geoPath)) {
    return false;
  }
  if (!RunGmsh(geoPath, mshPath, 2, gmshExe)) {
    return false;
  }
  std::vector<Eigen::Vector3f> nodes;
  std::vector<std::array<int, 4>> quads;
  if (!ParseMsh2Quads(mshPath, nodes, quads)) {
    return false;
  }
  cap.capNodes = std::move(nodes);
  cap.quads = std::move(quads);
  cap.sweepLayers = 1;
  std::cout << "[GmshCapMesher] cylindrical inner cap: " << cap.quads.size()
            << " quads\n";
  return true;
}
