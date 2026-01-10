#ifndef __COLOR_IMPLEMENTER_H__
#define __COLOR_IMPLEMENTER_H__

#include "CTMesh.h"
#include <algorithm>
#include <stack>

class Implementer {
public:
  Implementer(MeshLib::CTMesh *mesh);

private:
  MeshLib::CTMesh *mesh;
  void Implement();
  void DFS(MeshLib::CToolFace *face, std::map<int, int> *colorMap);
  void Paint(std::map<int, int> colorMap);
};

#endif
