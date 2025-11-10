#ifndef __SWEEPDIRDECTOR_H__
#define __SWEEPDIRDECTOR_H__

#include "CTMesh.h"
#include "PrimeData.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>

class SweepDirDetector {
public:
  std::vector<Eigen::Vector3f> GetSweepDir() { return this->SweepDir; };
  SweepDirDetector(MeshLib::CTMesh *mesh, std::vector<PrimeData> primes);

protected:
  MeshLib::CTMesh *mesh;
  std::vector<Eigen::Vector3f> SweepDir;
  std::vector<PrimeData> primes;
};

#endif
