#ifndef METAL_NEAREST_POINT_H
#define METAL_NEAREST_POINT_H

#include <Eigen/Dense>
#include <vector>

void ComputeNearestPointsMetal(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<float>> &PointList,
    std::vector<std::vector<std::vector<int>>> &NearestIndex);

void ComputeNearestPointsCPU(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<float>> &PointList,
    std::vector<std::vector<std::vector<int>>> &NearestIndex);

#endif
