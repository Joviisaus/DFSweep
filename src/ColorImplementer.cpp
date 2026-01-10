#include "ColorImplementer.h"
#include "CTMesh.h"
#include "Mesh/iterators.h"
#include <unordered_map>

Implementer::Implementer(MeshLib::CTMesh *mesh) {
  this->mesh = mesh;
  Implement();
}

void Implementer::Implement() {
  for (MeshLib::MeshFaceIterator mviter(this->mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mviter.value());
    if (f->rgb().norm() < 0.1)
      f->visited() = false;
    else
      f->visited() = true;
  }

  std::map<int, int> colorMap;
  colorMap.clear();
  for (MeshLib::MeshFaceIterator mviter(this->mesh); !mviter.end(); ++mviter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mviter.value());

    if (!f->visited()) {
      colorMap.clear();

      DFS(f, &colorMap);

      if (!colorMap.empty()) {
        Paint(colorMap);
      }
    }
    colorMap.clear();
  }
}

void Implementer::DFS(MeshLib::CToolFace *startFace,
                      std::map<int, int> *colorMap) {
  if (!startFace || startFace->visited())
    return;

  std::stack<MeshLib::CToolFace *> s;
  s.push(startFace);
  startFace->visited() = true; // 入栈即标记

  while (!s.empty()) {
    MeshLib::CToolFace *curr = s.top();
    s.pop();

    for (MeshLib::CTMesh::FaceHalfedgeIterator fhiter(curr); !fhiter.end();
         ++fhiter) {
      MeshLib::CToolFace *f =
          static_cast<MeshLib::CToolFace *>(fhiter.value()->he_sym()->face());

      if (!f)
        continue;

      if (f->visited()) {
        if (f->rgb().norm() > 0.3) {
          (*colorMap)[f->sweeplabel()]++; // 注意 colorMap 是指针，写法已修正
        }
        continue;
      }

      // 关键：在入栈前检查并标记
      f->visited() = true;
      s.push(f);
    }
  }
}
void Implementer::Paint(std::map<int, int> colorMap) {
  int paintcolor = std::max_element(colorMap.begin(), colorMap.end(),
                                    [](const std::pair<int, int> &a,
                                       const std::pair<int, int> &b) {
                                      return a.second < b.second;
                                    })
                       ->first;
  for (MeshLib::MeshFaceIterator mfiter(this->mesh); !mfiter.end(); ++mfiter) {
    MeshLib::CToolFace *f = static_cast<MeshLib::CToolFace *>(mfiter.value());
    if (f->sweeplabel() == -1 && f->visited() == true)
      f->sweeplabel() = paintcolor;
  }
}
