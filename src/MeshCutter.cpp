
#include "MeshCutter.h"
#include "CTMesh.h"
#include "Mesh/iterators.h"

MeshCutter::MeshCutter(
    MeshLib::CTMesh *mesh,
    std::vector<std::map<int, Eigen::Vector3f>> CuttingHexLists) {
  this->mesh = mesh;
  this->CuttingHexLists = CuttingHexLists;
  this->MeshCut();
}

// --- 基础几何辅助函数 ---
bool isPointInPolygon(const Eigen::Vector3f &P,
                      const std::vector<Eigen::Vector3f> &poly) {
  if (poly.size() < 3)
    return false;
  Eigen::Vector3f normal = (poly[1] - poly[0]).cross(poly[2] - poly[0]);
  for (size_t i = 0; i < poly.size(); ++i) {
    Eigen::Vector3f edge = poly[(i + 1) % poly.size()] - poly[i];
    Eigen::Vector3f toP = P - poly[i];
    if (normal.dot(edge.cross(toP)) < -1e-4)
      return false;
  }
  return true;
}

void MeshCutter::MeshCut() {
  const int faceIndices[6][4] = {{0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
                                 {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5}};

  int nextVertexId = -1;
  int nextFaceId = -1;
  for (MeshLib::MeshVertexIterator mviter(mesh); !mviter.end(); mviter++)
    if (mviter.value()->id() > nextVertexId)
      nextVertexId = mviter.value()->id();
  for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); mfiter++) {
    if (mfiter.value()->id() > nextFaceId)
      nextFaceId = mfiter.value()->id();
  }

  for (auto &boxVertices : this->CuttingHexLists) {
    for (int f_idx = 0; f_idx < 6; ++f_idx) {
      // 获取六面体当前面片
      std::vector<Eigen::Vector3f> hexFace(4);
      for (int i = 0; i < 4; ++i)
        hexFace[i] = boxVertices.at(faceIndices[f_idx][i]);
      Eigen::Vector3f nQuad =
          (hexFace[1] - hexFace[0]).cross(hexFace[2] - hexFace[0]).normalized();

      // --- 第一步：手动分裂所有相交边 ---
      std::vector<MeshLib::CToolEdge *> edges;
      for (MeshLib::MeshEdgeIterator eeiter(mesh); !eeiter.end(); ++eeiter)
        edges.push_back(static_cast<MeshLib::CToolEdge *>(eeiter.value()));

      for (auto *e : edges) {
        // 如果边已被删除（可能在之前的分裂中处理过），跳过
        if (e == NULL)
          continue;

        MeshLib::CToolVertex *v1 =
            static_cast<MeshLib::CToolVertex *>(e->halfedge(0)->source());
        MeshLib::CToolVertex *v2 =
            static_cast<MeshLib::CToolVertex *>(e->halfedge(0)->target());
        Eigen::Vector3f p1(v1->point()[0], v1->point()[1], v1->point()[2]);
        Eigen::Vector3f p2(v2->point()[0], v2->point()[1], v2->point()[2]);

        Eigen::Vector3f intersect;
        if (intersectSegmentPlane(p1, p2, hexFace[0], nQuad, intersect)) {
          if (isPointInPolygon(intersect, hexFace)) {
            if ((intersect - p1).norm() > 1e-4 &&
                (intersect - p2).norm() > 1e-4) {
              ManualSplitEdge(e, intersect, nextVertexId, nextFaceId);
            }
          }
        }
      }

      // --- 第二步：在面内连接新产生的顶点 ---
      std::vector<MeshLib::CToolFace *> faces;
      for (MeshLib::MeshFaceIterator mfiter(mesh); !mfiter.end(); ++mfiter)
        faces.push_back(static_cast<MeshLib::CToolFace *>(mfiter.value()));

      for (auto *f : faces) {
        std::vector<MeshLib::CToolVertex *> onPlaneVerts;
        // 遍历面片所有顶点，找出落在当前切割平面上的点
        for (MeshLib::CTMesh::FaceVertexIterator fviter(f); !fviter.end();
             ++fviter) {
          MeshLib::CToolVertex *v =
              static_cast<MeshLib::CToolVertex *>(fviter.value());
          Eigen::Vector3f p(v->point()[0], v->point()[1], v->point()[2]);
          if (std::abs(nQuad.dot(p - hexFace[0])) < 1e-4) {
            onPlaneVerts.push_back(v);
          }
        }

        // 如果面内有两个点都在平面上，说明需要切开这个面
        if (onPlaneVerts.size() == 2) {
          ManualSplitFace(f, onPlaneVerts[0], onPlaneVerts[1], nextFaceId);
        }
      }
    }
  }
}

void MeshCutter::ManualSplitEdge(MeshLib::CToolEdge *e, Eigen::Vector3f pos,
                                 int &vId, int &fId) {
  // 1. 提取邻接面和端点
  MeshLib::CToolHalfEdge *he0 =
      static_cast<MeshLib::CToolHalfEdge *>(e->halfedge(0));
  MeshLib::CToolHalfEdge *he1 =
      (e->halfedge(0) != NULL)
          ? static_cast<MeshLib::CToolHalfEdge *>(e->halfedge(1))
          : nullptr;

  MeshLib::CToolFace *f0 = static_cast<MeshLib::CToolFace *>(he0->face());
  MeshLib::CToolFace *f1 =
      (he1) ? static_cast<MeshLib::CToolFace *>(he1->face()) : nullptr;

  MeshLib::CToolVertex *v_src =
      static_cast<MeshLib::CToolVertex *>(he0->source());
  MeshLib::CToolVertex *v_tgt =
      static_cast<MeshLib::CToolVertex *>(he0->target());

  std::vector<MeshLib::CToolVertex *> loop0;
  std::vector<MeshLib::CToolVertex *> loop1;
  std::vector<MeshLib::CToolVertex *> loop2;
  std::vector<MeshLib::CToolVertex *> loop3;

  // SafeDeleteFace(f0);
  // SafeDeleteFace(f1);

  loop0.clear();
  loop1.clear();
  loop2.clear();
  loop3.clear();
  // 4. 创建新顶点并重建面
  MeshLib::CToolVertex *nv =
      static_cast<MeshLib::CToolVertex *>(this->mesh->createVertex(++vId));
  nv->point() = CPoint(pos.x(), pos.y(), pos.z());

  loop0.push_back(v_src);
  loop0.push_back(
      static_cast<MeshLib::CToolVertex *>(he1->he_next()->target()));
  loop0.push_back(nv);

  loop1.push_back(nv);
  loop1.push_back(
      static_cast<MeshLib::CToolVertex *>(he1->he_next()->target()));
  loop1.push_back(v_tgt);

  loop2.push_back(nv);
  loop2.push_back(
      static_cast<MeshLib::CToolVertex *>(he0->he_prev()->source()));
  loop2.push_back(v_src);

  loop3.push_back(v_tgt);
  loop3.push_back(
      static_cast<MeshLib::CToolVertex *>(he0->he_prev()->source()));
  loop3.push_back(nv);

  this->mesh->deleteFace(f0);
  this->mesh->deleteFace(f1);

  this->mesh->edges().remove_if([](MeshLib::CEdge *e) {
    // 先判空避免野指针/空指针访问，再判断目标条件
    return e != nullptr && e->halfedge(0) == NULL;
  });

  v_src->edges().remove_if([](MeshLib::CEdge *e) {
    // 先判空避免野指针/空指针访问，再判断目标条件
    return e != nullptr && e->halfedge(0) == NULL;
  });
  v_tgt->edges().remove_if([](MeshLib::CEdge *e) {
    // 先判空避免野指针/空指针访问，再判断目标条件
    return e != nullptr && e->halfedge(0) == NULL;
  });
  if (!loop0.empty())
    this->mesh->createFace(loop0, ++fId);
  if (!loop1.empty())
    this->mesh->createFace(loop1, ++fId);
  if (!loop2.empty())
    this->mesh->createFace(loop2, ++fId);
  if (!loop3.empty())
    this->mesh->createFace(loop3, ++fId);
}
void MeshCutter::ManualSplitFace(MeshLib::CToolFace *f,
                                 MeshLib::CToolVertex *v1,
                                 MeshLib::CToolVertex *v2, int &fId) {
  // 1. 提取当前面片的所有顶点循环 (可能是三角形、四边形或多边形)
  std::vector<MeshLib::CToolVertex *> loop;
  MeshLib::CToolHalfEdge *he =
      static_cast<MeshLib::CToolHalfEdge *>(f->halfedge());
  MeshLib::CToolHalfEdge *start = he;
  do {
    loop.push_back(static_cast<MeshLib::CToolVertex *>(he->source()));
    he = static_cast<MeshLib::CToolHalfEdge *>(he->he_next());
  } while (he != start);

  // 2. 在循环中定位 v1 和 v2
  auto it1 = std::find(loop.begin(), loop.end(), v1);
  auto it2 = std::find(loop.begin(), loop.end(), v2);
  if (it1 == loop.end() || it2 == loop.end())
    return;

  // 确保顺序，方便切分
  if (std::distance(loop.begin(), it1) > std::distance(loop.begin(), it2))
    std::swap(it1, it2);

  // 3. 将原多边形切分为两个子多边形
  std::vector<MeshLib::CToolVertex *> poly1, poly2;

  // 子多边形 1: v1 -> ... -> v2
  for (auto it = it1; it <= it2; ++it)
    poly1.push_back(*it);

  // 子多边形 2: v2 -> ... -> 结尾 -> 开头 -> v1
  for (auto it = it2; it != loop.end(); ++it)
    poly2.push_back(*it);
  for (auto it = loop.begin(); it <= it1; ++it)
    poly2.push_back(*it);

  // 4. 删除原面，并对结果进行三角化处理
  mesh->deleteFace(f);

  // 核心：调用三角化辅助函数
  TriangulateAndCreateFaces(poly1, fId);
  TriangulateAndCreateFaces(poly2, fId);
}

bool MeshCutter::intersectSegmentPlane(const Eigen::Vector3f &A,
                                       const Eigen::Vector3f &B,
                                       const Eigen::Vector3f &p0,
                                       const Eigen::Vector3f &normal,
                                       Eigen::Vector3f &intersect) {
  Eigen::Vector3f dir = B - A;
  float denom = normal.dot(dir);
  if (std::abs(denom) < 1e-7)
    return false; // 平行

  float t = normal.dot(p0 - A) / denom;
  if (t >= 0.0f && t <= 1.0f) {
    intersect = A + t * dir;
    return true;
  }
  return false;
}

// 主求交函数
bool MeshCutter::IntersectTriangleWithQuad(
    const std::vector<Eigen::Vector3f> &tri,
    const std::vector<Eigen::Vector3f> &quad,
    std::vector<Eigen::Vector3f> &outPoints) {
  outPoints.clear();

  // 1. 计算两个面的法线
  Eigen::Vector3f nTri = (tri[1] - tri[0]).cross(tri[2] - tri[0]).normalized();
  Eigen::Vector3f nQuad =
      (quad[1] - quad[0]).cross(quad[2] - quad[0]).normalized();

  // 2. 检查三角形的边与四边形平面的交点
  for (int i = 0; i < 3; ++i) {
    Eigen::Vector3f p;
    if (intersectSegmentPlane(tri[i], tri[(i + 1) % 3], quad[0], nQuad, p)) {
      if (isPointInPolygon(p, quad))
        outPoints.push_back(p);
    }
  }

  // 3. 检查四边形的边与三角形平面的交点
  for (int i = 0; i < 4; ++i) {
    Eigen::Vector3f p;
    if (intersectSegmentPlane(quad[i], quad[(i + 1) % 4], tri[0], nTri, p)) {
      if (isPointInPolygon(p, tri))
        outPoints.push_back(p);
    }
  }

  // 去重（处理顶点正好落在边上的情况）
  for (size_t i = 0; i < outPoints.size(); ++i) {
    for (size_t j = i + 1; j < outPoints.size(); ++j) {
      if ((outPoints[i] - outPoints[j]).norm() < 1e-5) {
        outPoints.erase(outPoints.begin() + j);
        --j;
      }
    }
  }

  return outPoints.size() >= 2;
}
void MeshCutter::TriangulateAndCreateFaces(
    std::vector<MeshLib::CToolVertex *> &poly, int &fId) {
  if (poly.size() < 3)
    return;

  if (poly.size() == 3) {
    // 已经是三角形
    this->mesh->createFace(poly, fId++);
  } else if (poly.size() == 4) {
    // 四边形拆分：选择较短的对角线
    Eigen::Vector3f v0(poly[0]->point()[0], poly[0]->point()[1],
                       poly[0]->point()[2]);
    Eigen::Vector3f v1(poly[1]->point()[0], poly[1]->point()[1],
                       poly[1]->point()[2]);
    Eigen::Vector3f v2(poly[2]->point()[0], poly[2]->point()[1],
                       poly[2]->point()[2]);
    Eigen::Vector3f v3(poly[3]->point()[0], poly[3]->point()[1],
                       poly[3]->point()[2]);

    if ((v0 - v2).squaredNorm() < (v1 - v3).squaredNorm()) {
      std::vector<MeshLib::CToolVertex *> t1 = {poly[0], poly[1], poly[2]};
      std::vector<MeshLib::CToolVertex *> t2 = {poly[0], poly[2], poly[3]};
      this->mesh->createFace(t1, fId++);
      this->mesh->createFace(t2, fId++);
    } else {
      std::vector<MeshLib::CToolVertex *> t1 = {poly[0], poly[1], poly[3]};
      std::vector<MeshLib::CToolVertex *> t2 = {poly[1], poly[2], poly[3]};
      this->mesh->createFace(t1, fId++);
      this->mesh->createFace(t2, fId++);
    }
  } else {
    // 多边形扇形拆分
    for (size_t i = 1; i < poly.size() - 1; ++i) {
      std::vector<MeshLib::CToolVertex *> tri = {poly[0], poly[i], poly[i + 1]};
      this->mesh->createFace(tri, fId++);
    }
  }
}

void MeshCutter::SafeDeleteFace(MeshLib::CToolFace *f) {
  if (!f)
    return;

  // 1. 记录该面关联的所有边和顶点
  std::vector<MeshLib::CToolEdge *> edgesInFace;
  std::vector<MeshLib::CToolVertex *> vertsInFace;

  MeshLib::CToolHalfEdge *start =
      static_cast<MeshLib::CToolHalfEdge *>(f->halfedge());
  MeshLib::CToolHalfEdge *curr = start;
  do {
    edgesInFace.push_back(static_cast<MeshLib::CToolEdge *>(curr->edge()));
    vertsInFace.push_back(static_cast<MeshLib::CToolVertex *>(curr->source()));
    curr = static_cast<MeshLib::CToolHalfEdge *>(curr->he_next());
  } while (curr != start);

  // 2. 执行基础的 deleteFace
  // 注意：在标准库中，deleteFace 会移除该面对应的半边，
  // 但可能不会自动从顶点或边的数据结构中完全“摘除”不再使用的边
  this->mesh->deleteFace(f);

  // this->mesh->edges().remove_if(
  //     [](MeshLib::CEdge *e) { return e != nullptr && e->halfedge(0) == NULL;
  //     });
}
