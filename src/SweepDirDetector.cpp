#include "SweepDirDetector.h"

SweepDirDetector::SweepDirDetector(MeshLib::CTMesh *mesh,
                                   std::vector<PrimeData> primes) {
  if (primes.size() == 0) {
    std::cout << "Initial Patch Required" << std::endl;
    return;
  }
  this->mesh = mesh;
  this->primes = primes;
  PlaneFeatureLineClustering();
  SweepDirCleaning();
}

void SweepDirDetector::PlaneFeatureLineClustering() {
  this->PlaneFeatureLine.clear();
  for (MeshLib::MeshEdgeIterator meiter(this->mesh); !meiter.end(); meiter++) {
    MeshLib::CToolEdge *edge =
        static_cast<MeshLib::CToolEdge *>(meiter.value());
    int label1 =
        static_cast<MeshLib::CToolFace *>(edge->halfedge(0)->face())->label();

    int label2 =
        static_cast<MeshLib::CToolFace *>(edge->halfedge(1)->face())->label();
    if (label1 != label2 && primes[label1].isPlane && primes[label2].isPlane) {
      Eigen::Vector3f Normal1 =
          Eigen::Vector3f(primes[label1].params[1], primes[label1].params[2],
                          primes[label1].params[3]);
      Eigen::Vector3f Normal2 =
          Eigen::Vector3f(primes[label2].params[1], primes[label2].params[2],
                          primes[label2].params[3]);
      this->PlaneFeatureLine.push_back(Normal1.cross(Normal2));
    }
  }

  for (int i = 0; i < this->primes.size(); i++) {
    for (int j = i + 1; j < this->primes.size(); j++) {
      if (primes[i].isPlane && primes[j].isPlane) {
        Eigen::Vector3f Normal1 = Eigen::Vector3f(
            primes[i].params[1], primes[i].params[2], primes[i].params[3]);
        Eigen::Vector3f Normal2 = Eigen::Vector3f(
            primes[j].params[1], primes[j].params[2], primes[j].params[3]);
        if (acos(abs(Normal1.dot(Normal2))) < 1e-16) {
          this->PlaneFeatureLine.push_back(Normal1);
        }
      }
    }
  }
}

void SweepDirDetector::SweepDirCleaning() {
  SweepDir.clear();

  if (PlaneFeatureLine.empty()) {
    return;
  }

  std::vector<Eigen::Vector3f> uniqueDirections;
  std::vector<bool> isUsed(PlaneFeatureLine.size(), false);

  const float angleTolerance = 0.087f;
  const float cosTolerance = std::cos(angleTolerance);

  for (size_t i = 0; i < PlaneFeatureLine.size(); ++i) {
    if (isUsed[i]) {
      continue;
    }

    Eigen::Vector3f normalizedDir = PlaneFeatureLine[i].normalized();
    uniqueDirections.push_back(normalizedDir);
    isUsed[i] = true;

    for (size_t j = i + 1; j < PlaneFeatureLine.size(); ++j) {
      if (isUsed[j]) {
        continue;
      }

      Eigen::Vector3f dir1 = PlaneFeatureLine[i].normalized();
      Eigen::Vector3f dir2 = PlaneFeatureLine[j].normalized();

      float dotProduct = std::abs(dir1.dot(dir2));

      if (dotProduct > cosTolerance) {
        isUsed[j] = true;
      }
    }
  }

  SweepDir = std::move(uniqueDirections);
}
