#include "SweepDirFilter.h"

SweepDirFilter::SweepDirFilter(
    std::vector<Eigen::Vector3f> *SweepDir,
    std::vector<std::vector<std::vector<std::vector<float>>>> *SweepProjScalar,

    std::vector<std::vector<std::vector<int>>> FieldLabel) {
  this->SweepProjScalar = SweepProjScalar;
  this->SweepDir = SweepDir;
  this->MarkedSweep.resize(SweepDir->size());
  this->RestField.resize(this->SweepProjScalar[0][0].size());
  this->xSize = this->SweepProjScalar[0][0].size();
  this->ySize = this->SweepProjScalar[0][0][0].size();
  this->zSize = this->SweepProjScalar[0][0][0][0].size();
  this->FieldLabel = FieldLabel;
  this->RestField.resize(this->xSize);
  for (size_t i = 0; i < this->xSize; ++i) {
    this->RestField[i].resize(this->ySize);
    for (size_t j = 0; j < this->ySize; ++j) {
      this->RestField[i][j].resize(this->zSize, 1.0f);
    }
  }

  this->RestSize = xSize * ySize * zSize;
  for (auto MS : this->MarkedSweep) {
    MS = 0;
  }
  this->SweepDirFilting();
  this->SweepDirCleaning();
}

int SweepDirFilter::Filting() {
  std::vector<int> FitCount;
  FitCount.resize(SweepDir->size());

  // Update Scalar for each Sweep Direction.
  for (int i = 0; i < this->SweepDir->size(); i++) {
    if (this->MarkedSweep[i] == 1) {
      FitCount[i] = INT_MAX;
      continue;
    }
    int FitSize = 0;
    for (int x = 0; x < this->xSize; x++) {
      for (int y = 0; y < this->ySize; y++) {
        for (int z = 0; z < this->zSize; z++) {
          if (this->SweepProjScalar[0][i][x][y][z] > RotateZero &&
              this->RestField[x][y][z] == 1)
            FitSize++;
        }
      }
    }
    FitCount[i] = FitSize;
  }

  int MinDirID =
      std::min_element(FitCount.begin(), FitCount.end()) - FitCount.begin();
  this->MarkedSweep[MinDirID] = 1;

  this->RestSize = 0;
  for (int x = 0; x < this->xSize; x++) {
    for (int y = 0; y < this->ySize; y++) {
      for (int z = 0; z < this->zSize; z++) {
        if (this->SweepProjScalar[0][MinDirID][x][y][z] < RotateZero) {
          RestField[x][y][z] = 0;
        }
        if (this->FieldLabel[x][y][z] < 0)
          RestField[x][y][z] = 0;
        if (RestField[x][y][z] == 1)
          RestSize++;
      }
    }
  }
  return MinDirID;
}

void SweepDirFilter::SweepDirFilting() {
  while (this->RestSize > 0) {
    int RestSize_Former = this->RestSize;
    int Current_Poped_Filter_id = Filting();
    std::cout << "Rest Size: " << this->RestSize << std::endl;
    if (RestSize_Former - RestSize < 1e-6 * xSize * ySize * zSize) {
      this->MarkedSweep[Current_Poped_Filter_id] = 0;
      break;
    }
    std::cout << "Poped Filter id: " << Current_Poped_Filter_id << std::endl;
    bool hasZero;
    hasZero = false;
    for (auto MS : this->MarkedSweep) {
      std::cout << MS << " ";
      if (MS == 0)
        hasZero = true;
    }

    if (!hasZero) {
      break;
    }
  }

  std::cout << "Sweep Dir After Filting" << std::endl;
  for (int i = 0; i < this->MarkedSweep.size(); i++) {
    if (this->MarkedSweep[i] == 1)
      std::cout << std::endl << SweepDir[0][i] << std::endl;
  }
}

void SweepDirFilter::SweepDirCleaning() {
  std::vector<int> newMarkedSweep;
  std::vector<Eigen::Vector3f> newSweepDir;
  std::vector<std::vector<std::vector<std::vector<float>>>> newSweepProjScalar;

  newMarkedSweep.reserve(MarkedSweep.size());
  newSweepDir.reserve(MarkedSweep.size());
  newSweepProjScalar.reserve(MarkedSweep.size());

  for (size_t i = 0; i < MarkedSweep.size(); ++i) {
    if (MarkedSweep[i] == 1) {

      newMarkedSweep.push_back(MarkedSweep[i]);
      newSweepDir.push_back((*SweepDir)[i]);
      newSweepProjScalar.push_back((*SweepProjScalar)[i]);
    }
  }

  MarkedSweep = std::move(newMarkedSweep);
  (*SweepDir) = std::move(newSweepDir);
  (*SweepProjScalar) = std::move(newSweepProjScalar);
}
