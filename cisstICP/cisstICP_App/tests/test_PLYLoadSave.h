#ifndef TEST_PLYLOADSAVE_H
#define TEST_PLYLOADSAVE_H

#include "cisstMesh.h"

void test_PLYLoadSave()
{
  cisstMesh mesh;
  mesh.LoadPLY("C://workspace//cisstICP//test_data//ProximalFemur2.ply");

  std::cout << "v: " << mesh.vertices.size() << std::endl
    << " f: " << mesh.faces.size() << std::endl
    << " fn: " << mesh.faceNormals.size() << std::endl
    << " fnbr: " << mesh.faceNeighbors.size() << std::endl
    << " vn: " << mesh.vertexNormals.size() << std::endl;

  std::cout << mesh.vertices(0) << std::endl
    << mesh.faces(0) << std::endl;

  if (mesh.faceNormals.size() > 0) {
    std::cout << mesh.faceNormals(0) << std::endl;
  }
  if (mesh.faceNeighbors.size() > 0) {
    std::cout << mesh.faceNeighbors(0) << std::endl;
  }
  if (mesh.vertexNormals.size() > 0) {
    std::cout << mesh.vertexNormals(0) << std::endl;
  }

  mesh.SavePLY("C://workspace//cisstICP//test_data//ProximalFemur_new.ply");
}

#endif // TEST_PLYLOADSAVE_H