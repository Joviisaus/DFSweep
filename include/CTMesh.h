#ifndef _CTMESH_H
#define _CTMESH_H

#include <map>
#include <vector>

#include "Geometry/Point.h"
#include "Geometry/Point2.h"
#include "Mesh/BaseMesh.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Vertex.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"

using namespace std;

namespace MeshLib {
class CToolVertex : public CVertex {
public:
  CToolVertex(){};
  ~CToolVertex(){};
  void _to_string();
  void _from_string();

  CPoint &rgb() { return m_rgb; };
  CPoint &normal() { return m_normal; };
  int &label() { return m_label; };
  bool &FeaturePoint() { return m_FeaturePoint; };
  bool &marked() { return m_marked; };

protected:
  CPoint m_rgb;
  CPoint m_normal;
  int m_label;
  bool m_FeaturePoint;
  bool m_marked;
};

inline void CToolVertex::_from_string() {
  CParser parser(m_string);

  for (std::list<CToken *>::iterator iter = parser.tokens().begin();
       iter != parser.tokens().end(); ++iter) {
    CToken *token = *iter;
    if (token->m_key == "Label") {
      int l;
      sscanf(token->m_value.c_str(), "(%d)", &l);
      this->m_label = l;
    }
  }
}

inline void CToolVertex::_to_string() {
  std::string a;
  m_string = a;

  if (1) {
    CParser parser3(m_string);
    parser3._removeToken("rgb");
    parser3._toString(m_string);
    std::stringstream iss3;
    iss3 << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ")";
    if (m_string.length() > 0) {
      m_string += " ";
    }
    m_string += iss3.str();
  }
  if (1) {
    CParser parser4(m_string);
    parser4._removeToken("normal");
    parser4._toString(m_string);
    std::stringstream iss4;
    iss4 << "normal=(" << m_normal[0] << " " << m_normal[1] << " "
         << m_normal[2] << ")";
    if (m_string.length() > 0) {
      m_string += " ";
    }
    m_string += iss4.str();
  }
  if (1) {
    CParser parser5(m_string);
    parser5._removeToken("Label");
    parser5._toString(m_string);
    std::stringstream iss5;
    iss5 << "Label=(" << m_label << ")";
    if (m_string.length() > 0) {
      m_string += " ";
    }
    m_string += iss5.str();
  }
}

class CToolEdge : public CEdge {
public:
  CToolEdge(){};
  ~CToolEdge(){};
  void _to_string();
  void _from_string();
  bool &sharp() { return this->m_sharp; };
  double &angel() { return this->m_angel; };

protected:
  bool m_sharp;
  double m_angel;
};

inline void CToolEdge::_to_string() {
  std::string a;
  m_string = a;

  CParser parser(m_string);
  parser._removeToken("sharp");
  parser._toString(m_string);
  stringstream iss;
  iss << "sharp";
  if (m_string.length() > 0) {
    m_string += " ";
  }
  if (m_sharp == true) {
    m_string += iss.str();
  }
}

inline void CToolEdge::_from_string() {
  CParser parser(m_string);

  for (std::list<CToken *>::iterator iter = parser.tokens().begin();
       iter != parser.tokens().end(); ++iter) {
    CToken *token = *iter;
    if (token->m_key == "sharp") {
      this->m_sharp = true;
    }
  }
}

class CToolFace : public CFace {
public:
  CToolFace(){};
  ~CToolFace(){};
  void _to_string();
  void _from_string();
  bool &visited() { return m_visited; };
  CPoint &normal() { return m_normal; };
  CPoint &rgb() { return m_rgb; };
  double &area() { return m_area; };
  int &label() { return m_label; };
  int &sweeplabel() { return m_sweeplabel; }

protected:
  bool m_visited;
  CPoint m_normal;
  double m_area;
  CPoint m_rgb;
  int m_label;
  int m_sweeplabel;
};

inline void CToolFace::_from_string() {
  CParser parser(m_string);

  for (std::list<CToken *>::iterator iter = parser.tokens().begin();
       iter != parser.tokens().end(); ++iter) {
    CToken *token = *iter;
    if (token->m_key == "label") {
      int l;
      sscanf(token->m_value.c_str(), "(%d)", &l);
      this->m_label = l;
    }
  }
}

inline void CToolFace::_to_string() {
  std::string a;
  m_string = a;

  CParser parser(m_string);
  parser._removeToken("label");
  parser._toString(m_string);
  stringstream iss;
  iss << "label=(" << m_label << ") ";
  if (m_string.length() > 0) {
    m_string += " ";
  }
  m_string += iss.str();

  CParser parser2(m_string);
  parser._removeToken("rgb");
  parser._toString(m_string);
  stringstream iss2;
  iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ") ";
  if (m_string.length() > 0) {
    m_string += " ";
  }
  m_string += iss.str();
}

class CToolHalfEdge : public CHalfEdge {
public:
  CToolHalfEdge(){};
  ~CToolHalfEdge(){};
  double &angle() { return m_angle; };
  void _to_string();

protected:
  double m_angle;
};

inline void CToolHalfEdge::_to_string() {
  // std::string a;
  // m_string = a;
}

template <typename V, typename E, typename F, typename H>
class CToolMesh : public CBaseMesh<V, E, F, H> {
public:
  typedef V CVertex;
  typedef E CEdge;
  typedef F CFace;
  typedef H CHalfEdge;

  typedef CBoundary<V, E, F, H> CBoundary;
  typedef CLoop<V, E, F, H> CLoop;

#define RIDGE 0;
#define VALLEY 1;
#define SADDLE 2;
#define FLAT 3;

  typedef MeshVertexIterator<V, E, F, H> MeshVertexIterator;
  typedef MeshEdgeIterator<V, E, F, H> MeshEdgeIterator;
  typedef VertexVertexIterator<V, E, F, H> VertexVertexIterator;
  typedef VertexEdgeIterator<V, E, F, H> VertexEdgeIterator;
  typedef MeshFaceIterator<V, E, F, H> MeshFaceIterator;
  typedef FaceVertexIterator<V, E, F, H> FaceVertexIterator;
  typedef VertexFaceIterator<V, E, F, H> VertexFaceIterator;
  typedef FaceHalfedgeIterator<V, E, F, H> FaceHalfedgeIterator;
  typedef VertexOutHalfedgeIterator<V, E, F, H> VertexOutHalfedgeIterator;
  typedef VertexInHalfedgeIterator<V, E, F, H> VertexInHalfedgeIterator;
  typedef FaceEdgeIterator<V, E, F, H> FaceEdgeIterator;
};

typedef CToolMesh<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> CTMesh;

} // namespace MeshLib
#endif
