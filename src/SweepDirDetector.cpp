#include "SweepDirDetector.h"

SweepDirDetector::SweepDirDetector(MeshLib::CTMesh *mesh,
                                   std::vector<PrimeData> primes) {
  if (primes.size() == 0) {
    std::cout << "Initial Patch Required" << std::endl;
    return;
  }
}
