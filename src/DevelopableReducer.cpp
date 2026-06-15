#include "DevelopableReducer.h"
#include "Mesh/iterators.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

Eigen::Matrix3d Adjugate(const Eigen::Matrix3d &H) {
  Eigen::Matrix3d adj;
  adj(0, 0) = H(1, 1) * H(2, 2) - H(1, 2) * H(2, 1);
  adj(0, 1) = H(0, 2) * H(2, 1) - H(0, 1) * H(2, 2);
  adj(0, 2) = H(0, 1) * H(1, 2) - H(0, 2) * H(1, 1);
  adj(1, 0) = H(1, 2) * H(2, 0) - H(1, 0) * H(2, 2);
  adj(1, 1) = H(0, 0) * H(2, 2) - H(0, 2) * H(2, 0);
  adj(1, 2) = H(0, 2) * H(1, 0) - H(0, 0) * H(1, 2);
  adj(2, 0) = H(1, 0) * H(2, 1) - H(1, 1) * H(2, 0);
  adj(2, 1) = H(0, 1) * H(2, 0) - H(0, 0) * H(2, 1);
  adj(2, 2) = H(0, 0) * H(1, 1) - H(0, 1) * H(1, 0);
  return adj;
}

std::vector<Eigen::Vector3d>
CollectPatchPoints(const PrimeData &prime, MeshLib::CTMesh *mesh) {
  std::vector<Eigen::Vector3d> pts;
  if (!mesh) {
    return pts;
  }
  for (MeshLib::MeshVertexIterator it(mesh); !it.end(); ++it) {
    auto *v = static_cast<MeshLib::CToolVertex *>(it.value());
    if (v->label() != prime.id || v->FeaturePoint()) {
      continue;
    }
    auto p = v->point();
    pts.emplace_back(p[0], p[1], p[2]);
  }
  return pts;
}

int PrimeIndexById(const std::vector<PrimeData> &primes, int id) {
  for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
    if (primes[i].id == id) {
      return i;
    }
  }
  return -1;
}

double FitResidualRms(const PrimeData &prime,
                      const std::vector<Eigen::Vector3d> &points) {
  if (points.empty()) {
    return std::numeric_limits<double>::max();
  }
  double sum = 0.0;
  for (const auto &p : points) {
    double f = DevelopableReducer::EvalF(prime.params, p);
    sum += f * f;
  }
  return std::sqrt(sum / points.size());
}

std::vector<double> BuildCylinderParams(const Eigen::Vector3d &center,
                                        Eigen::Vector3d axis, double radius) {
  axis.normalize();
  const double nx = axis.x(), ny = axis.y(), nz = axis.z();
  const double cx = center.x(), cy = center.y(), cz = center.z();

  std::vector<double> p(10, 0.0);
  p[7] = 1.0 - nx * nx;
  p[8] = 1.0 - ny * ny;
  p[9] = 1.0 - nz * nz;
  p[4] = -2.0 * nx * ny;
  p[5] = -2.0 * nx * nz;
  p[6] = -2.0 * ny * nz;

  // F = |u|² - (u·n)² - R², u = x - C
  p[1] = -2.0 * cx * p[7] - cy * p[4] - cz * p[5];
  p[2] = -cx * p[4] - 2.0 * cy * p[8] - cz * p[6];
  p[3] = -cx * p[5] - cy * p[6] - 2.0 * cz * p[9];

  double quadAtC = cx * cx * p[7] + cy * cy * p[8] + cz * cz * p[9] +
                   cx * cy * p[4] + cx * cz * p[5] + cy * cz * p[6];
  double linAtC = p[1] * cx + p[2] * cy + p[3] * cz;
  p[0] = -(quadAtC + linAtC) - radius * radius;

  return p;
}

Eigen::Matrix3d QuadraticPart(const std::vector<double> &params) {
  Eigen::Matrix3d Q;
  Q << 2 * params[7], params[4], params[5], params[4], 2 * params[8],
      params[6], params[5], params[6], 2 * params[9];
  return Q;
}

} // namespace

double DevelopableReducer::EvalF(const std::vector<double> &params,
                                   const Eigen::Vector3d &p) {
  if (params.size() < 10) {
    return 0.0;
  }
  double x = p.x(), y = p.y(), z = p.z();
  return params[0] + params[1] * x + params[2] * y + params[3] * z +
         params[4] * x * y + params[5] * x * z + params[6] * y * z +
         params[7] * x * x + params[8] * y * y + params[9] * z * z;
}

Eigen::Vector3d DevelopableReducer::EvalGrad(const std::vector<double> &params,
                                             const Eigen::Vector3d &p) {
  double x = p.x(), y = p.y(), z = p.z();
  return Eigen::Vector3d(
      params[1] + params[4] * y + params[5] * z + 2 * params[7] * x,
      params[2] + params[4] * x + params[6] * z + 2 * params[8] * y,
      params[3] + params[5] * x + params[6] * y + 2 * params[9] * z);
}

Eigen::Matrix3d DevelopableReducer::EvalHessian(const std::vector<double> &params) {
  Eigen::Matrix3d H;
  H << 2 * params[7], params[4], params[5], params[4], 2 * params[8], params[6],
      params[5], params[6], 2 * params[9];
  return H;
}

double DevelopableReducer::GaussianCurvature(const std::vector<double> &params,
                                             const Eigen::Vector3d &p) {
  Eigen::Vector3d g = EvalGrad(params, p);
  double gn = g.norm();
  if (gn < 1e-12) {
    return 0.0;
  }
  Eigen::Matrix3d H = EvalHessian(params);
  Eigen::Matrix3d adjH = Adjugate(H);
  double gn4 = gn * gn * gn * gn;
  double term = g.dot(adjH * g) - gn * gn * H.determinant();
  return term / gn4;
}

double DevelopableReducer::MeanAbsGaussianCurvature(const PrimeData &prime,
                                                    MeshLib::CTMesh *mesh) {
  auto pts = CollectPatchPoints(prime, mesh);
  if (pts.empty()) {
    return 0.0;
  }
  double sum = 0.0;
  int count = 0;
  for (const auto &p : pts) {
    double k = GaussianCurvature(prime.params, p);
    if (std::isfinite(k)) {
      sum += std::abs(k);
      count++;
    }
  }
  return count > 0 ? sum / count : 0.0;
}

PrimeData DevelopableReducer::FitPlane(
    const std::vector<Eigen::Vector3d> &points) {
  PrimeData out;
  out.params.assign(10, 0.0);
  out.isPlane = true;
  out.rank = 1;
  out.residual = 0.0;

  if (points.size() < 3) {
    return out;
  }

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const auto &p : points) {
    centroid += p;
  }
  centroid /= static_cast<double>(points.size());

  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto &p : points) {
    Eigen::Vector3d d = p - centroid;
    cov += d * d.transpose();
  }
  cov /= static_cast<double>(points.size());

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  Eigen::Vector3d normal = es.eigenvectors().col(0);
  if (normal.norm() < 1e-12) {
    normal = Eigen::Vector3d::UnitZ();
  } else {
    normal.normalize();
  }

  out.params[1] = normal.x();
  out.params[2] = normal.y();
  out.params[3] = normal.z();
  out.params[0] = -normal.dot(centroid);
  return out;
}

PrimeData DevelopableReducer::FitCylinder(
    const std::vector<Eigen::Vector3d> &points) {
  PrimeData out;
  out.params.assign(10, 0.0);
  out.isPlane = false;
  out.rank = 2;
  out.residual = 0.0;

  if (points.size() < 4) {
    return FitPlane(points);
  }

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const auto &p : points) {
    centroid += p;
  }
  centroid /= static_cast<double>(points.size());

  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto &p : points) {
    Eigen::Vector3d d = p - centroid;
    cov += d * d.transpose();
  }
  cov /= static_cast<double>(points.size());

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  Eigen::Vector3d axis = es.eigenvectors().col(2);
  if (axis.norm() < 1e-12) {
    axis = Eigen::Vector3d::UnitZ();
  } else {
    axis.normalize();
  }

  double r2 = 0.0;
  for (const auto &p : points) {
    Eigen::Vector3d u = p - centroid;
    Eigen::Vector3d perp = u - u.dot(axis) * axis;
    r2 += perp.squaredNorm();
  }
  double radius = std::sqrt(std::max(r2 / points.size(), 1e-12));

  out.params = BuildCylinderParams(centroid, axis, radius);
  return out;
}

PrimeData DevelopableReducer::ReduceToDevelopable(const PrimeData &prime,
                                                  MeshLib::CTMesh *mesh) {
  return ReduceToDevelopable(prime, mesh, Config());
}

PrimeData DevelopableReducer::ReduceToDevelopable(const PrimeData &prime,
                                                  MeshLib::CTMesh *mesh,
                                                  const Config &cfg) {
  PrimeData result = prime;
  auto points = CollectPatchPoints(prime, mesh);
  if (points.size() < 3) {
    return result;
  }

  if (prime.isPlane) {
    return result;
  }

  double meanK = MeanAbsGaussianCurvature(prime, mesh);
  Eigen::Matrix3d Q = QuadraticPart(prime.params);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(Q);
  Eigen::Vector3d ev = es.eigenvalues();
  double evMax = std::max({std::abs(ev[0]), std::abs(ev[1]), std::abs(ev[2])});
  double evMin = std::min({std::abs(ev[0]), std::abs(ev[1]), std::abs(ev[2])});
  double axisRatio = evMax > 1e-12 ? evMin / evMax : 0.0;

  bool coneLike = axisRatio < cfg.coneAxisRatio;
  bool needsReduction = meanK > cfg.curvatureThreshold;

  if (!needsReduction && !coneLike) {
    return result;
  }

  PrimeData plane = FitPlane(points);
  PrimeData cylinder = FitCylinder(points);

  double rmsPlane = FitResidualRms(plane, points);
  double rmsCylinder = FitResidualRms(cylinder, points);

  bool pickPlane = rmsPlane <= rmsCylinder;

  // 锥面 / 狭长三角：展宽较小 → 平面；否则柱面
  if (coneLike) {
    double planarity = evMin / (evMax + 1e-12);
    if (planarity < cfg.planePlanarityRatio) {
      pickPlane = true;
    } else {
      pickPlane = false;
    }
  }

  result = pickPlane ? plane : cylinder;
  result.id = prime.id;
  result.rank = pickPlane ? 1 : 2;
  result.residual =
      pickPlane ? rmsPlane : rmsCylinder;

  std::cout << "[DevelopableReducer] prime id=" << prime.id << " mean|K|="
            << meanK << " -> " << (pickPlane ? "plane" : "cylinder")
            << " rms=" << result.residual << "\n";

  return result;
}

void DevelopableReducer::ReduceAll(std::vector<PrimeData> &primes,
                                   MeshLib::CTMesh *mesh) {
  ReduceAll(primes, mesh, Config());
}

void DevelopableReducer::ReduceAll(std::vector<PrimeData> &primes,
                                   MeshLib::CTMesh *mesh, const Config &cfg) {
  if (!mesh || primes.empty()) {
    return;
  }

  for (auto &prime : primes) {
    prime = ReduceToDevelopable(prime, mesh, cfg);
  }

  UpdateVertexNormals(mesh, primes);
}

void DevelopableReducer::UpdateVertexNormals(
    MeshLib::CTMesh *mesh, const std::vector<PrimeData> &primes) {
  if (!mesh) {
    return;
  }

  for (MeshLib::MeshVertexIterator it(mesh); !it.end(); ++it) {
    auto *v = static_cast<MeshLib::CToolVertex *>(it.value());
    if (v->FeaturePoint()) {
      continue;
    }

    int label = v->label();
    int idx = PrimeIndexById(primes, label);
    if (idx < 0 || idx >= static_cast<int>(primes.size())) {
      continue;
    }

    const auto &params = primes[idx].params;
    if (params.size() < 10) {
      continue;
    }

    auto vertPoint = v->point();
    Eigen::Vector3d p(vertPoint[0], vertPoint[1], vertPoint[2]);
    Eigen::Vector3d n = EvalGrad(params, p);
    if (n.norm() < 1e-12) {
      continue;
    }
    n.normalize();

    Eigen::Vector3f former(v->normal()[0], v->normal()[1], v->normal()[2]);
    if (former.dot(n.cast<float>()) < 0) {
      n = -n;
    }
    v->normal()[0] = static_cast<float>(n.x());
    v->normal()[1] = static_cast<float>(n.y());
    v->normal()[2] = static_cast<float>(n.z());
  }
}
