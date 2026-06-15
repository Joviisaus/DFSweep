#include <metal_stdlib>
using namespace metal;

struct Vec3 {
    float x, y, z;
};

kernel void nearestPointKernel(
    device const Vec3 *coords   [[buffer(0)]],
    device const Vec3 *points   [[buffer(1)]],
    device int        *indices  [[buffer(2)]],
    constant int      &numPoints [[buffer(3)]],
    uint               gid      [[thread_position_in_grid]])
{
    Vec3 query = coords[gid];
    float minDist = 1e30f;
    int bestIdx = -1;

    for (int i = 0; i < numPoints; ++i) {
        float dx = query.x - points[i].x;
        float dy = query.y - points[i].y;
        float dz = query.z - points[i].z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < minDist) {
            minDist = d2;
            bestIdx = i;
        }
    }

    indices[gid] = bestIdx;
}
