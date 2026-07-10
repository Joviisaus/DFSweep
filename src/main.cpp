#include "CLI/CLI.hpp"
#include "DFContainer.h"
#include "MeshViewer.h"
#include <CLI/CLI.hpp>
#include <iostream>

int main(int argc, char **argv) {
  string input_file;
  string prime_file;
  CLI::App app{"DFSweep input"};
  app.add_option("-i,--input", input_file,
                 "Input File, Manifold in .obj and .m is Supported")
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("-e,--epsilon", epsilon,
                 "Gradiance Zero Value,default as 5e-1f");
  app.add_option("-s,--SampleSize", SampleSize,
                 "Grid Sample Size,default is 100");
  app.add_option("-r,--RotateZero", RotateZero,
                 "RotateZero is difining max rotation angel in each sweep "
                 "direction,default is 0.1");
  app.add_option("-d,--ParallelAngel", ParallelAngel,
                 "Adjust if a pair of plane is parallel,default is 1e-2");
  app.add_option("-a,--alpha", Alpha,
                 "Parament in Energy Function,default is 0.5");
  app.add_option("-p,--primefile", prime_file,
                 "Prime File if a Smooth Field is Needed")
      ->check(CLI::ExistingFile);

  bool useGeneralizedSweep = false;
  bool cutMesh = false;
  bool forceTwoSweepBodies = false;
  float sweepAngleThreshold = 0.3f;
  app.add_flag("--cut-mesh", cutMesh,
               "Cut mesh with sweep boxes (slow; off by default)");
  app.add_flag("-g,--generalized-sweep", useGeneralizedSweep,
               "Use generalized iso-surface sweep decomposition (SweepBlock) "
               "instead of CuttingBox");
  app.add_flag(
      "--two-sweep-bodies", forceTwoSweepBodies,
      "Force tube+base two-body split (CylinderPrime). Default is planar "
      "CuttingBox only — non-planar primes do not switch to cylinder sweep");
  app.add_option("-t,--sweep-threshold", sweepAngleThreshold,
                 "Angular threshold (radians) for sweep direction constraint "
                 "in generalized mode");

  bool reduceDevelopable = false;
  double developableKThreshold = 1e-4;
  std::string developableExport;
  app.add_flag(
      "--developable",
      reduceDevelopable,
      "Replace non-developable quadric patches (K!=0) with plane/cylinder (K=0)");
  app.add_option("--developable-k", developableKThreshold,
                 "Mean |Gaussian curvature| threshold to trigger reduction");
  app.add_option("--developable-out", developableExport,
                 "Export reduced prime parameters to this text file");

  CLI11_PARSE(app, argc, argv);
  MeshLib::CTMesh mesh;
  std::cout << "Loading mesh:" << input_file << std::endl;
  if (input_file.substr(input_file.find_last_of(".") + 1) == "obj") {
    mesh.read_obj(input_file.c_str());
  } else if (input_file.substr(input_file.find_last_of(".") + 1) == "m") {
    mesh.read_m(input_file.c_str());
  } else {
    std::cout << "Unsupported file format. Please use .obj or .m files."
              << std::endl;
    return 1;
  }
  std::cout << "Mesh loaded successfully." << std::endl;
  std::cout << "Mesh info:" << std::endl;
  std::cout << "Number of vertices: " << mesh.numVertices() << std::endl;
  std::cout << "Number of faces: " << mesh.numFaces() << std::endl;

#ifdef ENABLE_CUDA
  std::cout << "CUDA detected" << std::endl;
#endif
#ifdef ENABLE_OMP
  std::cout << "openmp detected" << std::endl;
#endif

  DistanceField DistanceField;
  DistanceField.SetMesh(&mesh);
  if (!prime_file.empty() &&
      prime_file.substr(prime_file.find_last_of(".") + 1) == "txt") {
    DistanceField.readPrime(prime_file);
  }
  if (reduceDevelopable) {
    if (prime_file.empty()) {
      std::cerr << "--developable requires -p/--primefile.\n";
      return 1;
    }
    DistanceField.ReducePrimesToDevelopable(developableKThreshold,
                                            developableExport);
  }
  DistanceField.GridScalar(SampleSize);
  DistanceField.ComputeDistanceField();

  if (useGeneralizedSweep) {
    if (prime_file.empty()) {
      std::cerr << "Generalized sweep requires --primefile for surface patches."
                << std::endl;
      return 1;
    }
    DistanceField.GeneralizedSweepDecomposition(sweepAngleThreshold, false);
  } else if (!prime_file.empty()) {
    // 能量已由 ComputeSweepDirectionEnergies 算好 → 有几个能量就建几个 CuttingBox。
    // --two-sweep-bodies 仅用于显式管体+底座特例。
    if (forceTwoSweepBodies) {
      DistanceField.DecomposeIntoTwoSweepBodies(sweepAngleThreshold);
    } else {
      DistanceField.RunCuttingBoxPipeline(cutMesh);
    }
  }

  DistanceField.ApplySweepVisualization();
  DistanceField.exportPlanesToFile("PlaneLists.txt");
  // DistanceField.SaveFieldToBinary("distance_field.bin");
  std::string output_str =
      (std::filesystem::path(input_file).parent_path() /
       (std::filesystem::path(input_file).stem().string() + "_output" + ".m"))
          .string();
  mesh.write_m(output_str.c_str());

  std::cout << "DistanceField Build" << std::endl;

  MeshViewer viewer;
  viewer.setGrid(DistanceField.getField(), DistanceField.getGradianceCount(),
                 DistanceField.GetSweepProjScalar(),
                 DistanceField.GetSweepProjEnergy(),
                 DistanceField.getCuttingHex(), DistanceField.getSweepDir(),
                 DistanceField.ForbiddenBoundaryPoints,
                 DistanceField.getGradianceDiff(), DistanceField.getCoord(),
                 DistanceField.getSweepBlockNonPlanar(),
                 DistanceField.getSweepBlockColors(),
                 DistanceField.getDisplayHexIndices(),
                 DistanceField.getSweepEnergyNames());
  viewer.setMesh(&mesh);
  std::cout << "Grid and Mesh Settled" << std::endl;
  viewer.show();
  return 0;
}
