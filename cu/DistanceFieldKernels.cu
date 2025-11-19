// DistanceFieldKernels.cu
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <Eigen/Core>
#include <vector>
#include <cstdio>
#include <cstdlib>

// ============= CUDA 错误检查 =============
inline void checkCuda(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        printf("CUDA Error [%s]: %s\n", msg, cudaGetErrorString(err));
        std::exit(EXIT_FAILURE);
    }
}

// ============= Vec3 结构（三维向量）============
struct Vec3 {
    float x, y, z;

    __host__ __device__ Vec3() : x(0), y(0), z(0) {}
    __host__ __device__ Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    __host__ __device__ Vec3(const Eigen::Vector3f& e) : x(e.x()), y(e.y()), z(e.z()) {}

    __host__ __device__ Vec3 operator-(const Vec3& b) const {
        return Vec3(x - b.x, y - b.y, z - b.z);
    }
    __host__ __device__ float squaredNorm() const {
        return x * x + y * y + z * z;
    }
};

// ============= CUDA Kernel =============
__global__ void NearestPointKernel(const Vec3* d_coord,
    const Vec3* d_points,
    Vec3* d_result,
    int numCoord,
    int numPoints)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numCoord) return;

    Vec3 query = d_coord[idx];

    float minDist = 1e30f;
    Vec3 nearest = Vec3();

    for (int i = 0; i < numPoints; ++i) {
        float dist = (query - d_points[i]).squaredNorm();
        if (dist < minDist) {
            minDist = dist;
            nearest = d_points[i];
        }
    }

    d_result[idx] = nearest;
}

// ============= 对外调用函数 =============
void ComputeNearestPointsCUDA(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>>& Coord,
    const std::vector<std::vector<float>>& PointList,
    std::vector<std::vector<std::vector<Eigen::Vector3f>>>& NearestPoint)
{
    int xSize = (int)Coord.size();
    if (xSize == 0) return;
    int ySize = (int)Coord[0].size();
    int zSize = (int)Coord[0][0].size();
    int total = xSize * ySize * zSize;
    int numPoints = (int)PointList.size();

    // ---- flatten Coord ----
    std::vector<Vec3> flatCoord;
    flatCoord.reserve(total);
    for (int i = 0; i < xSize; ++i)
        for (int j = 0; j < ySize; ++j)
            for (int k = 0; k < zSize; ++k)
                flatCoord.emplace_back(Coord[i][j][k]);

    // ---- convert PointList ----
    std::vector<Vec3> ptList;
    ptList.reserve(numPoints);
    for (int i = 0; i < numPoints; ++i)
        ptList.emplace_back(PointList[i][0], PointList[i][1], PointList[i][2]);

    // ---- allocate device ----
    Vec3* d_coord = nullptr, * d_points = nullptr, * d_result = nullptr;
    checkCuda(cudaMalloc(&d_coord, sizeof(Vec3) * total), "Malloc d_coord");
    checkCuda(cudaMalloc(&d_points, sizeof(Vec3) * numPoints), "Malloc d_points");
    checkCuda(cudaMalloc(&d_result, sizeof(Vec3) * total), "Malloc d_result");

    // ---- copy to GPU ----
    checkCuda(cudaMemcpy(d_coord, flatCoord.data(), sizeof(Vec3) * total, cudaMemcpyHostToDevice), "Memcpy d_coord");
    if (numPoints > 0)
        checkCuda(cudaMemcpy(d_points, ptList.data(), sizeof(Vec3) * numPoints, cudaMemcpyHostToDevice), "Memcpy d_points");

    // ---- launch kernel ----
    const int threads = 256;
    const int blocks = (total + threads - 1) / threads;
    NearestPointKernel << <blocks, threads >> > (d_coord, d_points, d_result, total, numPoints);
    checkCuda(cudaGetLastError(), "Kernel launch");
    checkCuda(cudaDeviceSynchronize(), "Sync after kernel");

    // ---- copy results back ----
    std::vector<Vec3> flatResult(total);
    checkCuda(cudaMemcpy(flatResult.data(), d_result, sizeof(Vec3) * total, cudaMemcpyDeviceToHost), "Memcpy back");

    // ---- unflatten ----
    NearestPoint.resize(xSize);
    int idx = 0;
    for (int i = 0; i < xSize; ++i) {
        NearestPoint[i].resize(ySize);
        for (int j = 0; j < ySize; ++j) {
            NearestPoint[i][j].resize(zSize);
            for (int k = 0; k < zSize; ++k) {
                Vec3 v = flatResult[idx++];
                NearestPoint[i][j][k] = Eigen::Vector3f(v.x, v.y, v.z);
            }
        }
    }

    cudaFree(d_coord);
    cudaFree(d_points);
    cudaFree(d_result);
}
