#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "MetalNearestPoint.h"
#include <iostream>
#include <vector>

struct Vec3Metal {
    float x, y, z;
};

void ComputeNearestPointsMetal(
    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &Coord,
    const std::vector<std::vector<float>> &PointList,
    std::vector<std::vector<std::vector<int>>> &NearestIndex) {

    int xSize = (int)Coord.size();
    if (xSize == 0) return;
    int ySize = (int)Coord[0].size();
    int zSize = (int)Coord[0][0].size();
    int total = xSize * ySize * zSize;
    int numPoints = (int)PointList.size();

    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            std::cerr << "Metal: No GPU device found, falling back to CPU." << std::endl;
            extern void ComputeNearestPointsCPU(
                const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &,
                const std::vector<std::vector<float>> &,
                std::vector<std::vector<std::vector<int>>> &);
            ComputeNearestPointsCPU(Coord, PointList, NearestIndex);
            return;
        }

        std::cout << "Metal GPU: " << [device.name UTF8String] << std::endl;

        NSString *libPath = [[NSBundle mainBundle] pathForResource:@"NearestPoint" ofType:@"metallib"];
        id<MTLLibrary> library = nil;
        NSError *error = nil;

        if (libPath) {
            library = [device newLibraryWithFile:libPath error:&error];
        }

        if (!library) {
            NSString *srcPath = nil;
            NSArray<NSString *> *searchPaths = @[
                [[NSProcessInfo processInfo].arguments[0] stringByDeletingLastPathComponent],
                @".",
                @".."
            ];
            for (NSString *base in searchPaths) {
                NSString *candidate = [base stringByAppendingPathComponent:@"NearestPoint.metal"];
                if ([[NSFileManager defaultManager] fileExistsAtPath:candidate]) {
                    srcPath = candidate;
                    break;
                }
                candidate = [base stringByAppendingPathComponent:@"metal/NearestPoint.metal"];
                if ([[NSFileManager defaultManager] fileExistsAtPath:candidate]) {
                    srcPath = candidate;
                    break;
                }
                candidate = [base stringByAppendingPathComponent:@"../metal/NearestPoint.metal"];
                if ([[NSFileManager defaultManager] fileExistsAtPath:candidate]) {
                    srcPath = candidate;
                    break;
                }
            }

            if (!srcPath) {
                std::cerr << "Metal: Cannot find NearestPoint.metal shader file, falling back to CPU." << std::endl;
                extern void ComputeNearestPointsCPU(
                    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &,
                    const std::vector<std::vector<float>> &,
                    std::vector<std::vector<std::vector<int>>> &);
                ComputeNearestPointsCPU(Coord, PointList, NearestIndex);
                return;
            }

            NSString *source = [NSString stringWithContentsOfFile:srcPath
                                                         encoding:NSUTF8StringEncoding
                                                            error:&error];
            if (!source) {
                std::cerr << "Metal: Failed to read shader: " << [error.localizedDescription UTF8String] << std::endl;
                extern void ComputeNearestPointsCPU(
                    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &,
                    const std::vector<std::vector<float>> &,
                    std::vector<std::vector<std::vector<int>>> &);
                ComputeNearestPointsCPU(Coord, PointList, NearestIndex);
                return;
            }

            MTLCompileOptions *opts = [[MTLCompileOptions alloc] init];
            opts.fastMathEnabled = YES;
            library = [device newLibraryWithSource:source options:opts error:&error];
            if (!library) {
                std::cerr << "Metal: Shader compile error: " << [error.localizedDescription UTF8String] << std::endl;
                extern void ComputeNearestPointsCPU(
                    const std::vector<std::vector<std::vector<Eigen::Vector3f>>> &,
                    const std::vector<std::vector<float>> &,
                    std::vector<std::vector<std::vector<int>>> &);
                ComputeNearestPointsCPU(Coord, PointList, NearestIndex);
                return;
            }
        }

        id<MTLFunction> function = [library newFunctionWithName:@"nearestPointKernel"];
        if (!function) {
            std::cerr << "Metal: Kernel function not found." << std::endl;
            return;
        }

        id<MTLComputePipelineState> pipeline =
            [device newComputePipelineStateWithFunction:function error:&error];
        if (!pipeline) {
            std::cerr << "Metal: Pipeline error: " << [error.localizedDescription UTF8String] << std::endl;
            return;
        }

        std::vector<Vec3Metal> flatCoord;
        flatCoord.reserve(total);
        for (int i = 0; i < xSize; ++i)
            for (int j = 0; j < ySize; ++j)
                for (int k = 0; k < zSize; ++k) {
                    const auto &c = Coord[i][j][k];
                    flatCoord.push_back({c.x(), c.y(), c.z()});
                }

        std::vector<Vec3Metal> ptList;
        ptList.reserve(numPoints);
        for (int i = 0; i < numPoints; ++i)
            ptList.push_back({PointList[i][0], PointList[i][1], PointList[i][2]});

        id<MTLBuffer> coordBuf = [device newBufferWithBytes:flatCoord.data()
                                                     length:sizeof(Vec3Metal) * total
                                                    options:MTLResourceStorageModeShared];
        id<MTLBuffer> pointsBuf = [device newBufferWithBytes:ptList.data()
                                                      length:sizeof(Vec3Metal) * numPoints
                                                     options:MTLResourceStorageModeShared];
        id<MTLBuffer> indexBuf = [device newBufferWithLength:sizeof(int) * total
                                                    options:MTLResourceStorageModeShared];
        id<MTLBuffer> numPtsBuf = [device newBufferWithBytes:&numPoints
                                                      length:sizeof(int)
                                                     options:MTLResourceStorageModeShared];

        id<MTLCommandQueue> queue = [device newCommandQueue];
        id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [cmdBuf computeCommandEncoder];

        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:coordBuf   offset:0 atIndex:0];
        [encoder setBuffer:pointsBuf  offset:0 atIndex:1];
        [encoder setBuffer:indexBuf   offset:0 atIndex:2];
        [encoder setBuffer:numPtsBuf  offset:0 atIndex:3];

        NSUInteger threadGroupSize = pipeline.maxTotalThreadsPerThreadgroup;
        if (threadGroupSize > 256) threadGroupSize = 256;
        MTLSize gridSize = MTLSizeMake(total, 1, 1);
        MTLSize groupSize = MTLSizeMake(threadGroupSize, 1, 1);
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
        [encoder endEncoding];

        [cmdBuf commit];
        [cmdBuf waitUntilCompleted];

        if (cmdBuf.error) {
            std::cerr << "Metal: Execution error: "
                      << [cmdBuf.error.localizedDescription UTF8String] << std::endl;
            return;
        }

        int *resultPtr = (int *)[indexBuf contents];
        NearestIndex.resize(xSize);
        int idx = 0;
        for (int i = 0; i < xSize; ++i) {
            NearestIndex[i].resize(ySize);
            for (int j = 0; j < ySize; ++j) {
                NearestIndex[i][j].resize(zSize);
                for (int k = 0; k < zSize; ++k)
                    NearestIndex[i][j][k] = resultPtr[idx++];
            }
        }

        std::cout << "Metal: Nearest point computation complete." << std::endl;
    }
}
