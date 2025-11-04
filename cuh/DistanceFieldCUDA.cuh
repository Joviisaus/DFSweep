#ifndef DISTANCE_FIELD_CUDA_CUH
#define DISTANCE_FIELD_CUDA_CUH

#include <Eigen/Dense>
#include <cuda_runtime.h>
#include <vector>

class DistanceFieldCUDA {
public:
  // 构造函数和析构函数
  DistanceFieldCUDA();
  ~DistanceFieldCUDA();

  // 初始化CUDA设备
  static bool InitializeCUDA();

  // 计算距离场的主要接口
  bool ComputeDistanceField(
      const std::vector<std::vector<float>> &pointList,
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &coord,
      std::vector<std::vector<std::vector<float>>> &field);

  // 获取设备信息
  static void PrintDeviceInfo();

private:
  // CUDA内存管理
  bool AllocateDeviceMemory(int totalPoints, int totalVoxels);
  void FreeDeviceMemory();

  // 数据拷贝
  bool CopyDataToDevice(
      const std::vector<std::vector<float>> &pointList,
      const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &coord);
  bool CopyResultToHost(std::vector<std::vector<std::vector<float>>> &field);

  // 内核启动参数
  void SetupKernelParameters(int totalVoxels, dim3 &gridDims, dim3 &blockDims);

  // 设备指针
  float *d_points_ = nullptr;    // 点云数据 [x1, y1, z1, x2, y2, z2, ...]
  float *d_coords_ = nullptr;    // 体素坐标 [x1, y1, z1, x2, y2, z2, ...]
  float *d_distances_ = nullptr; // 输出距离场

  int num_points_ = 0;
  int num_voxels_ = 0;

  // 缓存用于性能优化
  float *h_points_ = nullptr;
  float *h_coords_ = nullptr;
};

// CUDA错误检查宏
#define CUDA_CHECK(call)                                                       \
  {                                                                            \
    cudaError_t cuda_status = call;                                            \
    if (cuda_status != cudaSuccess) {                                          \
      std::cerr << "CUDA Error: " << cudaGetErrorString(cuda_status) << " at " \
                << __FILE__ << ":" << __LINE__ << std::endl;                   \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

#endif
