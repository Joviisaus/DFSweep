#include "SweepDirSpliter.h"
#include "SweepDirFilter.h"

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

  // this->LabelTopoGen();
  this->SweepDirSplit();
}

void SweepDirSpliter::LabelTopoGen() {
  for (MeshLib::MeshVertexIterator mviter(this->mesh); !mviter.end();
       mviter++) {
    MeshLib::CToolVertex *vertex =
        static_cast<MeshLib::CToolVertex *>(mviter.value());
    int VertOrignLabel = vertex->label();
    for (MeshLib::CTMesh::VertexVertexIterator vviter(vertex); !vviter.end();
         ++vviter) {
      MeshLib::CToolVertex *vv =
          static_cast<MeshLib::CToolVertex *>(vviter.value());
      int label = vv->label();
      if (VertOrignLabel != label) {
        this->LabelTopo.coeffRef(label, VertOrignLabel) = 1;
        this->LabelTopo.coeffRef(VertOrignLabel, label) = 1;
      }
    }
  }
  for (int x = 0; x < this->SweepProjScalar[0][0].size(); x++) {
    for (int y = 0; y < this->SweepProjScalar[0][0][0].size(); y++) {
      for (int z = 0; z < this->SweepProjScalar[0][0][0][0].size(); z++) {
        int centerLabel = FieldLabel[x][y][z];
        if (centerLabel < 0)
          continue; // 无效标签跳过

        std::set<int> neighLabels;
        neighLabels.insert(centerLabel);

        const int dir[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                               {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};

        for (int k = 0; k < 6; k++) {
          int nx = x + dir[k][0];
          int ny = y + dir[k][1];
          int nz = z + dir[k][2];

          if (nx < 0 || ny < 0 || nz < 0 || nx >= FieldLabel.size() ||
              ny >= FieldLabel[0].size() || nz >= FieldLabel[0][0].size())
            continue;

          int neighLabel = FieldLabel[nx][ny][nz];
          if (neighLabel >= 0) {
            neighLabels.insert(neighLabel);
          }
        }

        if (neighLabels.size() > 1) {
          for (int l1 : neighLabels) {
            for (int l2 : neighLabels) {
              if (l1 != l2) {
                this->LabelTopo.coeffRef(l1, l2) = 1;
                this->LabelTopo.coeffRef(l2, l1) = 1;
              }
            }
          }
        }
      }
    }
  }
}

void SweepDirSpliter::SweepDirSplit() {
  std::vector<std::vector<int>> EnergyLabel;
  EnergyLabel.clear();
  EnergyLabel.resize(this->SweepDir->size());
  std::vector<std::unordered_set<int>> labelSet;
  labelSet.resize(this->SweepDir->size());
  for (int x = 0; x < this->SweepProjScalar[0][0].size(); x++) {
    for (int y = 0; y < this->SweepProjScalar[0][0][0].size(); y++) {
      for (int z = 0; z < this->SweepProjScalar[0][0][0][0].size(); z++) {
        float MaxEnergy = FLT_MIN;
        int MinEnergyDirLabel = -1;
        float MinEnergy = FLT_MAX;
        for (int SweepLabel = 0; SweepLabel < this->SweepDir->size();
             SweepLabel++) {
          float Energy = this->SweepProjScalar[0][SweepLabel][x][y][z];
          if (Energy < MinEnergy) {
            MinEnergy = Energy;
            MinEnergyDirLabel = SweepLabel;
          }
          if (Energy > MaxEnergy)
            MaxEnergy = Energy;
        }

        if (MaxEnergy > RotateZero)
          labelSet[MinEnergyDirLabel].insert(FieldLabel[x][y][z]);
      }
    }
  }
  for (int SweepLabel = 0; SweepLabel < this->SweepDir->size(); SweepLabel++) {
    EnergyLabel[SweepLabel].assign(labelSet[SweepLabel].begin(),
                                   labelSet[SweepLabel].end());
  }
  this->splitDisconnectedGroups(EnergyLabel, this->LabelTopo);
  this->SweepMask(EnergyLabel);
}

void SweepDirSpliter::splitDisconnectedGroups(
    std::vector<std::vector<int>> &EnergyLabel,
    const Eigen::MatrixXi &LabelTopo) {
  std::vector<std::vector<int>> newGroups;
  std::vector<Eigen::Vector3f> newSweepDirs;
  std::vector<std::vector<std::vector<std::vector<float>>>> newSweepProjScalar;

  for (int g = 0; g < (int)EnergyLabel.size(); ++g) {
    const std::vector<int> &group = EnergyLabel[g];
    if (group.empty())
      continue;

    Eigen::Vector3f sweepDirCopy = (*this->SweepDir)[g];
    std::vector<std::vector<std::vector<float>>> sweepProjCopy =
        (*this->SweepProjScalar)[g];

    std::vector<char> visited(group.size(), 0);

    for (int i = 0; i < (int)group.size(); ++i) {
      if (visited[i])
        continue;

      std::vector<int> connectedGroup;
      std::queue<int> q;
      q.push(i);
      visited[i] = 1;

      while (!q.empty()) {
        int idx = q.front();
        q.pop();
        connectedGroup.push_back(group[idx]);

        for (int j = 0; j < (int)group.size(); ++j) {
          if (!visited[j]) {
            int A = group[idx];
            int B = group[j];
            if (LabelTopo(A, B) != 0) {
              visited[j] = 1;
              q.push(j);
            }
          }
        }
      }
      newGroups.push_back(std::move(connectedGroup));
      newSweepDirs.push_back(sweepDirCopy);
      newSweepProjScalar.push_back(sweepProjCopy);
    }
  }

  EnergyLabel = std::move(newGroups);
  *this->SweepDir = std::move(newSweepDirs);
  *this->SweepProjScalar = std::move(newSweepProjScalar);
}

void SweepDirSpliter::SweepMask(std::vector<std::vector<int>> &EnergyLabel) {
  int dimX = this->SweepProjScalar[0][0].size();
  int dimY = this->SweepProjScalar[0][0][0].size();
  int dimZ = this->SweepProjScalar[0][0][0][0].size();
  int numDirs = this->SweepDir[0].size();

  std::vector<std::unordered_set<int>> dirLabelSets(numDirs);
  for (int dir = 0; dir < numDirs; dir++) {
    dirLabelSets[dir].insert(EnergyLabel[dir].begin(), EnergyLabel[dir].end());
  }

  std::vector<std::vector<std::vector<int>>> maxDirCache(
      dimX, std::vector<std::vector<int>>(dimY, std::vector<int>(dimZ, 0)));

#ifdef ENABLE_OMP

#pragma omp parallel for collapse(2)
  for (int x = 0; x < dimX; x++) {
    for (int y = 0; y < dimY; y++) {
      for (int z = 0; z < dimZ; z++) {
        int maxDir = 0;
        float maxVal = this->SweepProjScalar[0][0][x][y][z];

        for (int dir = 1; dir < numDirs; dir++) {
          float currentVal = this->SweepProjScalar[0][dir][x][y][z];
          if (currentVal > maxVal) {
            maxVal = currentVal;
            maxDir = dir;
          }
        }
        maxDirCache[x][y][z] = maxDir;
      }
    }
  }

// 并行替换值
#pragma omp parallel for collapse(2)
  for (int x = 0; x < dimX; x++) {
    for (int y = 0; y < dimY; y++) {
      for (int z = 0; z < dimZ; z++) {
        int label = this->FieldLabel[x][y][z];
        int maxDir = maxDirCache[x][y][z];
        float maxVal = this->SweepProjScalar[0][maxDir][x][y][z];

        for (int dir = 0; dir < numDirs; dir++) {
          if (dirLabelSets[dir].find(label) == dirLabelSets[dir].end()) {
            this->SweepProjScalar[0][dir][x][y][z] = maxVal;
          }
        }
      }
    }
  }
#else
  // 串行版本
  for (int x = 0; x < dimX; x++) {
    for (int y = 0; y < dimY; y++) {
      for (int z = 0; z < dimZ; z++) {
        int maxDir = 0;
        float maxVal = this->SweepProjScalar[0][0][x][y][z];

        for (int dir = 1; dir < numDirs; dir++) {
          float currentVal = this->SweepProjScalar[0][dir][x][y][z];
          if (currentVal > maxVal) {
            maxVal = currentVal;
            maxDir = dir;
          }
        }
        maxDirCache[x][y][z] = maxDir;
      }
    }
  }

  for (int x = 0; x < dimX; x++) {
    for (int y = 0; y < dimY; y++) {
      for (int z = 0; z < dimZ; z++) {
        int label = this->FieldLabel[x][y][z];
        int maxDir = maxDirCache[x][y][z];
        float maxVal = this->SweepProjScalar[0][maxDir][x][y][z];

        for (int dir = 0; dir < numDirs; dir++) {
          if (dirLabelSets[dir].find(label) == dirLabelSets[dir].end()) {
            this->SweepProjScalar[0][dir][x][y][z] = maxVal;
          }
        }
      }
    }
  }
#endif
}
