#ifndef __GMSH_CAP_MESHER_H__
#define __GMSH_CAP_MESHER_H__

#include "SweepFaceImprinter.h"
#include "SweepBlock.h"
#include <string>

/// 调用 Gmsh 对压印后的截面做 Quasi-structured 四边形剖分
class GmshCapMesher {
public:
  static bool IsAvailable();
  static std::string FindExecutable();

  /// 平移扫掠体：用压印 u/v 分割做 Transfinite/Recombine 四边形网格
  static bool MeshTranslationalCap(ImprintedCapMesh &cap,
                                   const SweepCapFrame &frame,
                                   float holeRadius,
                                   const std::string &workDir);

  /// 柱面径向扫掠体：在内柱面 (θ×轴向) 上生成四边形，再由外部沿径向扫掠
  static bool MeshCylindricalInnerCap(const SweepBlockRegion &block, int nTheta,
                                    int nAxial, ImprintedCapMesh &cap,
                                    const std::string &workDir);
};

#endif
