#ifndef __MESH_CUTTER_H__
#define __MESH_CUTTER_H__

#include "CTMesh.h"
#include "CuttingBox.h"

class MeshCutter {
public:
  MeshCutter(MeshLib::CTMesh *mesh,
             std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists);

private:
  MeshLib::CTMesh *mesh;
  std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists;
  void MeshCut();
  bool intersectSegmentPlane(const Eigen::Vector3f &A, const Eigen::Vector3f &B,
                             const Eigen::Vector3f &p0,
                             const Eigen::Vector3f &normal,
                             Eigen::Vector3f &intersect);
  bool IntersectTriangleWithQuad(const std::vector<Eigen::Vector3f> &tri,
                                 const std::vector<Eigen::Vector3f> &quad,
                                 std::vector<Eigen::Vector3f> &outPoints);
  void ManualSplitEdge(MeshLib::CToolEdge *e, Eigen::Vector3f pos, int &vId,
                       int &fId, int cfid);

  void ManualSplitFace(MeshLib::CToolFace *f, MeshLib::CToolVertex *v1,
                       MeshLib::CToolVertex *v2, int &fId);
  void TriangulateAndCreateFaces(std::vector<MeshLib::CToolVertex *> &poly,
                                 int &fId);
  void SafeDeleteFace(MeshLib::CToolFace *f);
};
#endif
