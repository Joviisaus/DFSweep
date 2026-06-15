#ifndef __DEVELOPABLE_REDUCER_H__
#define __DEVELOPABLE_REDUCER_H__

#include "CTMesh.h"
#include "PrimeData.h"
#include <Eigen/Eigen>
#include <vector>

/**
 * @brief 将高斯曲率 K≠0 的二次解析片简化为 K=0 的可展曲面（平面 / 柱面）。
 *
 * 隐式二次曲面 F(x)=0，在采样点上估计高斯曲率；若 |K| 超过阈值，
 * 对片上的网格顶点做平面 / 柱面最小二乘拟合，取残差更小者并写回
 * PrimeData 的 10 元参数。
 */
class DevelopableReducer {
public:
  struct Config {
    double curvatureThreshold;
    double planePlanarityRatio;
    double coneAxisRatio;
    Config()
        : curvatureThreshold(1e-4), planePlanarityRatio(0.02),
          coneAxisRatio(0.08) {}
  };

  static double EvalF(const std::vector<double> &params, const Eigen::Vector3d &p);
  static Eigen::Vector3d EvalGrad(const std::vector<double> &params,
                                   const Eigen::Vector3d &p);
  static Eigen::Matrix3d EvalHessian(const std::vector<double> &params);

  /** 隐式曲面 F=0 在 p 处的高斯曲率 */
  static double GaussianCurvature(const std::vector<double> &params,
                                  const Eigen::Vector3d &p);

  /** 在 mesh 上属于 prime.id 的顶点处平均 |K| */
  static double MeanAbsGaussianCurvature(const PrimeData &prime,
                                         MeshLib::CTMesh *mesh);

  static PrimeData FitPlane(const std::vector<Eigen::Vector3d> &points);
  static PrimeData FitCylinder(const std::vector<Eigen::Vector3d> &points);

  /** 单 patch 简化为可展曲面；若已可展则原样返回 */
  static PrimeData ReduceToDevelopable(const PrimeData &prime,
                                       MeshLib::CTMesh *mesh);
  static PrimeData ReduceToDevelopable(const PrimeData &prime,
                                       MeshLib::CTMesh *mesh,
                                       const Config &cfg);

  static void ReduceAll(std::vector<PrimeData> &primes, MeshLib::CTMesh *mesh);
  static void ReduceAll(std::vector<PrimeData> &primes, MeshLib::CTMesh *mesh,
                        const Config &cfg);

  static void UpdateVertexNormals(MeshLib::CTMesh *mesh,
                                  const std::vector<PrimeData> &primes);
};

#endif
