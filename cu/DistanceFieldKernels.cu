#include "DistanceFieldCUDA.cuh"
#include <cmath>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>


// 基础内核：每个线程处理一个体素，对所有点进行搜索
__global__ void ComputeDistanceFieldBruteForceKernel(const float *points,
                                                     int num_points,
                                                     const float *coords,
                                                     int num_voxels,
                                                     float *distances) {

  int voxel_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (voxel_idx >= num_voxels)
    return;

  // 获取当前体素坐标
  float voxel_x = coords[voxel_idx * 3];
  float voxel_y = coords[voxel_idx * 3 + 1];
  float voxel_z = coords[voxel_idx * 3 + 2];

  float min_distance_sq = FLT_MAX;

  // 遍历所有点，找到最小距离
  for (int i = 0; i < num_points; ++i) {
    float point_x = points[i * 3];
    float point_y = points[i * 3 + 1];
    float point_z = points[i * 3 + 2];

    float dx = point_x - voxel_x;
    float dy = point_y - voxel_y;
    float dz = point_z - voxel_z;

    float distance_sq = dx * dx + dy * dy + dz * dz;

    if (distance_sq < min_distance_sq) {
      min_distance_sq = distance_sq;
    }
  }

  distances[voxel_idx] = sqrtf(min_distance_sq);
}

// 优化的内核：使用共享内存缓存点数据
__global__ void ComputeDistanceFieldSharedKernel(const float *points,
                                                 int num_points,
                                                 const float *coords,
                                                 int num_voxels,
                                                 float *distances) {

  extern __shared__ float shared_points[];

  int voxel_idx = blockIdx.x * blockDim.x + threadIdx.x;
  int tid = threadIdx.x;

  // 将点数据加载到共享内存（分批处理）
  int points_per_block = blockDim.x;
  int num_batches = (num_points + points_per_block - 1) / points_per_block;

  float voxel_x = 0, voxel_y = 0, voxel_z = 0;

  if (voxel_idx < num_voxels) {
    voxel_x = coords[voxel_idx * 3];
    voxel_y = coords[voxel_idx * 3 + 1];
    voxel_z = coords[voxel_idx * 3 + 2];
  }

  float min_distance_sq = FLT_MAX;

  // 分批处理点数据
  for (int batch = 0; batch < num_batches; ++batch) {
    int point_offset = batch * points_per_block;
    int local_point_idx = point_offset + tid;

    // 加载点到共享内存
    if (local_point_idx < num_points) {
      shared_points[tid * 3] = points[local_point_idx * 3];
      shared_points[tid * 3 + 1] = points[local_point_idx * 3 + 1];
      shared_points[tid * 3 + 2] = points[local_point_idx * 3 + 2];
    }
    __syncthreads();

    // 处理当前批次中的点
    int points_in_batch = min(points_per_block, num_points - point_offset);

    for (int i = 0; i < points_in_batch; ++i) {
      if (voxel_idx >= num_voxels)
        continue;

      float point_x = shared_points[i * 3];
      float point_y = shared_points[i * 3 + 1];
      float point_z = shared_points[i * 3 + 2];

      float dx = point_x - voxel_x;
      float dy = point_y - voxel_y;
      float dz = point_z - voxel_z;

      float distance_sq = dx * dx + dy * dy + dz * dz;

      if (distance_sq < min_distance_sq) {
        min_distance_sq = distance_sq;
      }
    }
    __syncthreads();
  }

  if (voxel_idx < num_voxels) {
    distances[voxel_idx] = sqrtf(min_distance_sq);
  }
}

// 使用空间分割的优化内核（需要预先构建空间结构）
__global__ void ComputeDistanceFieldSpatialKernel(
    const float *points, int num_points, const float *coords, int num_voxels,
    const int *point_indices, const int *grid_structure, float *distances,
    int grid_size, float cell_size) {

  // 这里可以实现基于空间分割的优化算法
  // 需要预先在CPU上构建空间栅格结构
  // 实现相对复杂，这里提供框架
}
