#ifndef DISTANCE_FIELD_CUDA_CUH
#define DISTANCE_FIELD_CUDA_CUH
#include <Eigen/Dense>
#include <vector>

void ComputeNearestPointsCUDA(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>>& Coord,
    const std::vector<std::vector<float>>& PointList,
	std::vector<std::vector<std::vector<Eigen::Vector3f>>>& NearestPoint);
#endif
