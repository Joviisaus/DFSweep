# DFSweep

A computational tool for extracting sweep volumes through medial axis and medial surface extraction using distance fields.

## Overview

DFSweep utilizes distance field computations to extract medial axes (2D) and medial surfaces (3D), enabling precise sweep volume generation. This tool is designed for applications in geometric modeling, computer-aided design, and motion planning, leveraging efficient distance field algorithms for robust performance.

## Core Methodology

1. **Distance Field Calculation**: Compute signed distance fields (SDF) for input manifold geometries, representing the minimum distance from any point to the manifold boundary.

2. **Medial Feature Extraction**:
   - Identify medial axes (2D) and medial surfaces (3D) as loci of points with multiple equidistant boundary points, using gradient analysis of the distance field.
   - Refine features using a user-specified epsilon value for gradient zero-crossing detection.

3. **Sweep Volume Generation**: Combine extracted medial features with trajectory information to construct the complete sweep volume, with optional smooth field enhancement using prime files.

## Dependencies

- C++17 or later
- Eigen3 (linear algebra)
- Polyscope (visualization)
- CLI11 (command-line parsing)
- OpenMP (optional, for parallel acceleration)
- CUDA / Metal (optional, for GPU nearest-point acceleration)
- MeshLib (included as 3rd party)

## Build Instructions

```bash
git clone https://github.com/yourusername/dfsweep.git
cd dfsweep

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4

# Install (optional)
make install
```

### Notes

- OpenMP support is enabled automatically if found.
- On macOS, Metal GPU acceleration is used when available; on Windows with MSVC, CUDA may be used.
- Eigen3 will be downloaded automatically if not found on the system.
- After adding new source files, re-run `cmake ..` in the build directory.

## Usage

The executable is built as `DFSweep` (or `./build/DFSweep` from the build tree).

### Basic workflow

```bash
# Required: input mesh (.obj or .m)
./DFSweep -i input_manifold.obj

# With analytic surface patches (prime file, .txt)
./DFSweep -i model.m -p prime.txt
```

Output mesh is written next to the input as `<stem>_output.m`. Cutting planes are exported to `PlaneLists.txt`.

### Processing pipeline (order)

1. Load mesh → `SetMesh`
2. Optional: `readPrime` from `-p`
3. Optional: **developable reduction** (`--developable`) — must run before distance field
4. `GridScalar` + `ComputeDistanceField`
5. **Sweep decomposition** (default `CuttingBox` or `-g` generalized mode)
6. Visualization (Polyscope)

### Command-line options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input` | `-i` | *(required)* | Input mesh (`.obj` or `.m`) |
| `--primefile` | `-p` | — | Prime file (`.txt`) with quadric patch parameters |
| `--epsilon` | `-e` | `1e-2` | Gradient zero threshold |
| `--SampleSize` | `-s` | `100` | Grid resolution per axis |
| `--RotateZero` | `-r` | `0.1` | Max rotation angle per sweep direction |
| `--ParallelAngel` | `-d` | `1e-2` | Parallel-plane detection tolerance |
| `--alpha` | `-a` | `0.6` | Weight in sweep energy function |
| `--generalized-sweep` | `-g` | off | Use iso-surface sweep blocks instead of `CuttingBox` |
| `--sweep-threshold` | `-t` | `0.3` | Angle threshold (radians) for generalized sweep direction constraint |
| `--developable` | — | off | Replace non-developable quadrics with plane/cylinder |
| `--developable-k` | — | `1e-4` | Mean \|Gaussian curvature\| threshold to trigger reduction |
| `--developable-out` | — | — | Export simplified prime parameters to a text file |

---

### Default sweep decomposition (`CuttingBox`)

Fits oriented bounding boxes along detected sweep directions (requires prime file for full pipeline after distance field).

```bash
./DFSweep -i model.m -p prime.txt
```

---

### Generalized sweep decomposition (`SweepBlock`)

Alternative to `CuttingBox` for **non-axis-aligned** sweeps: translational, rotational, and radial patterns.

- Grows blocks **inward from analytic surface patches** along distance-field iso-surfaces.
- Keeps voxels where the **field direction** is perpendicular or tangent to the local iso-surface normal (within `-t` radians).
- Selects a **minimum set of blocks** to cover all interior voxels (greedy set cover).
- Validates sweep-body topology (outer base, inner base, side).

**Requires** `-p` / `--primefile`.

```bash
# Enable generalized mode
./DFSweep -i model.m -p prime.txt -g

# Stricter direction constraint (smaller angle, in radians)
./DFSweep -i model.m -p prime.txt -g --sweep-threshold 0.2
```

Approximate threshold guide:

| `--sweep-threshold` | Meaning |
|---------------------|---------|
| `0.3` (default) | ~17° deviation allowed from perpendicular or tangent |
| `0.15` | ~8.6° — tighter, smaller blocks |
| `0.5` | ~29° — looser, larger blocks |

---

### Developable surface reduction (`DevelopableReducer`)

Simplifies quadric patches with **non-zero Gaussian curvature** (e.g. ellipsoid, paraboloid) to **developable** surfaces (**K = 0**): **plane** or **cylinder**. Cone-like patches (rank-2 quadrics) are mapped to a plane if nearly flat along the axis, otherwise to a cylinder.

**Requires** `-p`. Runs **before** the distance field is computed.

```bash
./DFSweep -i model.m -p prime.txt --developable

# Custom curvature threshold and export new primes
./DFSweep -i model.m -p prime.txt \
  --developable \
  --developable-k 1e-3 \
  --developable-out primes_developable.txt
```

---

### Combined example

```bash
# Developable primes → generalized sweep → output + visualization
./DFSweep -i part.m -p prime.txt \
  --developable --developable-out primes_dev.txt \
  -g --sweep-threshold 0.3
```

---

### Programmatic API (C++)

```cpp
DistanceField df(&mesh);
df.readPrime("prime.txt");

// Optional: K≠0 → plane/cylinder
df.ReducePrimesToDevelopable(1e-4, "primes_dev.txt");

df.GridScalar(100);
df.ComputeDistanceField();

// Either:
df.RunCuttingBoxPipeline();                    // CuttingBox + MeshCutter
// Or:
df.GeneralizedSweepDecomposition(0.3f);        // SweepBlock + MeshCutter

auto blocks = df.GetSweepBlocks();
auto hexes  = df.getCuttingHex();
```

### Prime file

The prime file lists analytic patches (`m_primes[id]->GetParams()` with 10 coefficients per quadric). It is required for:

- Smooth / labeled distance field evaluation
- `--developable`
- `-g` / `--generalized-sweep`

Without `-p`, only basic distance-field visualization runs; sweep cutting and the new modules are not available.
