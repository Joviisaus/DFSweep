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
  app.add_option("-a,--ParallelAngel", ParallelAngel,
                 "Adjust if a pair of plane is parallel,default is 1e-2");
  app.add_option("-p,--primefile", prime_file,
                 "Prime File if a Smooth Field is Needed")
      ->check(CLI::ExistingFile);

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
  if (prime_file.substr(prime_file.find_last_of(".") + 1) == "txt") {
    DistanceField.readPrime(prime_file);
  }
  DistanceField.GridScalar(SampleSize);
  DistanceField.ComputeDistanceField();
  // DistanceField.SaveFieldToBinary("distance_field.bin");

  std::cout << "DistanceField Build" << std::endl;

  MeshViewer viewer;
  viewer.setGrid(DistanceField.getField(), DistanceField.getGradianceCount(),
                 DistanceField.GetSweepProjScalar(),
                 DistanceField.getGradianceDiff(), DistanceField.getCoord());
  viewer.setMesh(&mesh);
  std::cout << "Grid and Mesh Settled" << std::endl;
  viewer.show();
  return 0;
}
