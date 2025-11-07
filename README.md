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

- CUDA (optional, for GPU acceleration)

- MeshLib (included as 3rd party)

  ## Build Instructions

  ```bash
  git clone https://github.com/yourusername/dfsweep.git
  cd dfsweep

  # Create build directory

  mkdir build && cd build

  # Configure with CMake

  cmake .. -DCMAKE_BUILD_TYPE=Release

  # Build the project

  make -j4

  # Install (optional)

  make install
  ```

### Notes:

- OpenMP support is enabled automatically if found

- CUDA acceleration is optional and requires a CUDA-enabled compiler

- Eigen3 will be downloaded automatically if not found in the system

  ## Usage：

  ```bash
  # Basic usage with mandatory input file
  ./DFSweep -i input_manifold.obj

  # Specify epsilon value for gradient thresholding
  ./DFSweep --input complex_shape.m --epsilon 1e-1f

  # Use prime file for smooth field generation
  ./DFSweep -i base_shape.obj -p prime.txt
  ```
