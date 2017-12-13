// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************

#include "cisstMesh.h"
#include "utilities.h"

#include <cisstNumerical/nmrLSSolver.h>

#include <assert.h>
#undef NDEBUG       // enable assert in release mode

#include <fstream>

// compares vertex for storing in std::Map
//  NOTE: this routine is used for loading a mesh from an STL file;
//        it is used to detect multiple copies of the same vertex
//        in order to prevent storing the same vertex coordinate to
//        multiple locations in the vertices array
struct VertexCompare {
  bool operator() (const vct3 &k1, const vct3 &k2) const
  { // Return true if k1 goes before k2 in the strict weak 
    //  ordering for the map object;
    //  i.e. if k1 is not strictly less than k2, return false
    //       otherwise, return true
    //
    // Note: https://ece.uwaterloo.ca/~dwharder/aads/Relations/Weak_ordering/

    // order by x,y,z in that order
    if (k1.X() < k2.X())
      return true;
    else if (k1.X() == k2.X()) {
      if (k1.Y() < k2.Y())
        return true;
      else if (k1.Y() == k2.Y()) {
        if (k1.Z() < k2.Z())
          return true;
      }
    }
    return false;
  };
};


void cisstMesh::ResetMesh()
{
  vertices.SetSize(0);
  faces.SetSize(0);
  faceNormals.SetSize(0);
  vertexNormals.SetSize(0);
  faceNeighbors.SetSize(0);
  TriangleCov.SetSize(0);
  TriangleCovEig.SetSize(0);
}

void cisstMesh::InitializeNoiseModel()
{
  TriangleCov.SetSize(faces.size());
  TriangleCovEig.SetSize(faces.size());

  TriangleCov.SetAll(vct3x3(0.0));
  TriangleCovEig.SetAll(vct3(0.0));
}

void cisstMesh::InitializeNoiseModel(
  double noiseInPlaneVar,
  double noisePerpPlaneVar)
{
  vct3x3 M, M0;
  vctRot3 R;
  vct3 z(0.0, 0.0, 1.0);
  vct3 norm;

  if (faceNormals.size() != faces.size())
  {
    std::cout << "ERROR: must initialize face normals in order to compute mesh noise model" << std::endl;
    TriangleCov.SetSize(0);
    TriangleCovEig.SetSize(0);
    assert(0);
  }

  TriangleCov.SetSize(faces.size());
  TriangleCovEig.SetSize(faces.size());

  // set covariance eigenvalues (in order of decreasing magnitude)
  if (noiseInPlaneVar >= noisePerpPlaneVar)
  {
    TriangleCovEig.SetAll(vct3(noiseInPlaneVar, noiseInPlaneVar, noisePerpPlaneVar));
  }
  else
  {
    TriangleCovEig.SetAll(vct3(noisePerpPlaneVar, noiseInPlaneVar, noiseInPlaneVar));
  }

  // compute covariance matrices
  for (unsigned int i = 0; i < faces.size(); i++)
  {
    TriangleCov[i] = ComputePointCovariance(faceNormals[i], noisePerpPlaneVar, noiseInPlaneVar);
  }
}

void cisstMesh::SaveTriangleCovariances(std::string &filePath)
{
  std::cout << "Saving mesh covariances to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cout << "ERROR: failed to open file for saving cov: " << filePath << std::endl;
    assert(0);
  }
  unsigned int numCov = this->TriangleCov.size();
  //fs << "NUMCOV " << numCov << "\n";
  for (unsigned int i = 0; i < numCov; i++)
  {
    fs << this->TriangleCov.at(i).Row(0) << " "
      << this->TriangleCov.at(i).Row(1) << " "
      << this->TriangleCov.at(i).Row(2) << "\n";
  }
}

void cisstMesh::LoadPLY(const std::string &input_file) {
  ply_obj.read_ply_mesh(input_file,
    &vertices, &faces, &faceNormals, &faceNeighbors, &vertexNormals);
}

void cisstMesh::SavePLY(const std::string &output_file) {
  ply_obj.write_ply_mesh(output_file,
    &vertices, &faces, &faceNormals, &faceNeighbors, &vertexNormals);
}

int cisstMesh::LoadMesh(
  const vctDynamicVector<vct3> *vertices,
  const vctDynamicVector<vctInt3> *faces,
  const vctDynamicVector<vct3> *face_normals,
  const vctDynamicVector<vctInt3> *face_neighbors,
  const vctDynamicVector<vct3> *vertex_normals
  )
{
  ResetMesh();

  if (!vertices || !faces) {
    std::cout << "ERROR: vertices and faces must not be null" << std::endl;
    assert(0);
  }
  if (vertices->size() < 1 || faces->size() < 1)
  {
    std::cout << "ERROR: vertices and faces must not be empty" << std::endl;
    assert(0);
  }
  this->vertices = *vertices;
  this->faces = *faces;

  if (face_normals) {
    if (face_normals->size() != this->faces.size()) {
      std::cout << "ERROR: number of face normals does not equal number of faces" << std::endl;
      assert(0);
    }
    this->faceNormals = *face_normals;
  }

  if (face_neighbors) {
    if (face_neighbors->size() != this->faces.size()) {
      std::cout << "ERROR: number of face neighbors does not equal number of faces" << std::endl;
      assert(0);
    }
    this->faceNeighbors = *face_neighbors;
  }

  if (vertex_normals) {
    if (vertex_normals->size() != this->vertices.size()) {
      std::cout << "ERROR: number of face neighbors does not equal number of faces" << std::endl;
      assert(0);
    }
    this->vertexNormals = *vertex_normals;
  }

  InitializeNoiseModel();

  return 0;
}


//// -- Deprecated I/O --
//
//int  cisstMesh::LoadMesh(
//  const vctDynamicVector<vct3> &V,
//  const vctDynamicVector<vctInt3> &T)
//{
//  ResetMesh();
//
//  if (V.size() < 1 || T.size() < 1)
//  {
//    std::cout << "ERROR: invalid input" << std::endl;
//    assert(0);
//  }
//
//  vertices = V;
//  faces = T;
//  ComputeFaceNormalsFromVertices();
//
//  InitializeNoiseModel();
//
//  return 0;
//}
//
//int cisstMesh::LoadMesh(
//  const vctDynamicVector<vct3> &V,
//  const vctDynamicVector<vctInt3> &T,
//  const vctDynamicVector<vct3> &N)
//{
//  ResetMesh();
//
//  if (V.size() < 1 || T.size() < 1 || T.size() != N.size())
//  {
//    std::cout << "ERROR: invalid input" << std::endl;
//    assert(0);
//  }
//
//  vertices = V;
//  faces = T;
//  faceNormals = N;
//
//  InitializeNoiseModel();
//
//  return 0;
//}
//
//int cisstMesh::LoadMeshFile(const std::string &meshFilePath)
//{
//  int rv;
//
//  ResetMesh();
//  
//  rv = AddMeshFile(meshFilePath);
//
//  InitializeNoiseModel();
//
//  return rv;
//}
//
//int cisstMesh::LoadMeshFileMultiple(const std::vector<std::string> &meshFilePaths)
//{
//  int rv;
//
//  ResetMesh();
//  
//  std::vector<std::string>::const_iterator iter;
//  for (iter = meshFilePaths.begin(); iter != meshFilePaths.end(); iter++)
//  {
//    rv = AddMeshFile(*iter);
//    if (rv == -1) return -1;
//  }
//
//  InitializeNoiseModel();
//
//  return 0;
//}
//
//int cisstMesh::AddMeshFile(const std::string &meshFilePath)
//{
//  // Load mesh from ASCII file having format:
//  //
//  //  POINTS numPoints
//  //  px py pz
//  //   ...
//  //  px py pz
//  //  TRIANGLES numTriangles
//  //  vx1 vx2 vx3
//  //    ...
//  //  vx1 vx2 vx3
//  //  NORMALS numTriangles
//  //  nx ny nz
//  //    ...
//  //  nx ny nz
//  //
//  //  where vx's are indices into the points array
//
//  //std::cout << "Loading mesh file: " << meshFilePath << std::endl;
//
//  float f1, f2, f3;
//  int d1, d2, d3;
//  unsigned int itemsRead;
//  std::string line;
//
//  unsigned int vOffset = vertices.size();
//  unsigned int tOffset = faces.size();
//
//  // open file
//  std::ifstream meshFile;
//  meshFile.open(meshFilePath.c_str());
//  if (!meshFile.is_open())
//  {
//    std::cout << "ERROR: failed to open mesh file" << std::endl;
//    return -1;
//  }
//
//  // read points
//  //  (vertex coordinates)
//  unsigned int numPoints;
//  std::getline(meshFile, line);
//  itemsRead = std::sscanf(line.c_str(), "POINTS %u", &numPoints);
//  if (itemsRead != 1)
//  {
//    std::cout << "ERROR: expected POINTS header at line: " << line << std::endl;
//    return -1;
//  }
//  vct3 v;
//  vertices.resize(vOffset + numPoints);  // non destructive resize
//  unsigned int pointCount = 0;
//  while (meshFile.good() && pointCount < numPoints)
//  {
//    std::getline(meshFile, line);
//    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &f1, &f2, &f3);
//    if (itemsRead != 3)
//    {
//      std::cout << "ERROR: expected a point value at line: " << line << std::endl;
//      return -1;
//    }
//    v[0] = f1;
//    v[1] = f2;
//    v[2] = f3;
//    vertices.at(vOffset + pointCount).Assign(v);
//    pointCount++;
//  }
//  if (meshFile.bad() || meshFile.fail() || pointCount != numPoints)
//  {
//    std::cout << "ERROR: read points from mesh file failed; last line read: " << line << std::endl;
//    return -1;
//  }
//
//  // read triangles
//  //  (three values that index the points array, one for each vertex)
//  unsigned int numTriangles;
//  std::getline(meshFile, line);
//  itemsRead = std::sscanf(line.c_str(), "TRIANGLES %u", &numTriangles);
//  if (itemsRead != 1)
//  {
//    std::cout << "ERROR: expected TRIANGLES header at line: " << line << std::endl;
//    return -1;
//  }
//  vctInt3 T;
//  faces.resize(tOffset + numTriangles); // non destructive size change
//  unsigned int triangleCount = 0;
//  while (meshFile.good() && triangleCount < numTriangles)
//  {
//    std::getline(meshFile, line);
//    itemsRead = std::sscanf(line.c_str(), "%d %d %d", &d1, &d2, &d3);
//    if (itemsRead != 3)
//    {
//      std::cout << "ERROR: expeced three index values on line: " << line << std::endl;
//      return -1;
//    }
//    T.Assign(d1 + vOffset, d2 + vOffset, d3 + vOffset);
//    faces.at(tOffset + triangleCount) = T;
//    triangleCount++;
//  }
//  if (meshFile.bad() || meshFile.fail() || triangleCount != numTriangles)
//  {
//    std::cout << "ERROR: while reading triangles from mesh file; last line read: " << line << std::endl;
//    return -1;
//  }
//
//  // read triangle normals
//  unsigned int numNormals;
//  std::getline(meshFile, line);
//  itemsRead = std::sscanf(line.c_str(), "NORMALS %u", &numNormals);
//  if (numNormals != numTriangles)
//  {
//    std::cout << std::endl << "ERROR: number of triangles and normals in the mesh do not match" << std::endl << std::endl;
//    return -1;
//  }
//  if (itemsRead != 1)
//  {
//    std::cout << std::endl << " ===> WARNING: normal vectors are missing in this legacy mesh file <=== " << std::endl << std::endl;
//    return 0;
//  }
//  vct3 N;
//  faceNormals.resize(tOffset + numNormals); // non destructive size change
//  unsigned int normCount = 0;
//  while (meshFile.good() && normCount < numNormals)
//  {
//    std::getline(meshFile, line);
//    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &f1, &f2, &f3);
//    if (itemsRead != 3)
//    {
//      std::cout << "ERROR: expeced three decimal values on line: " << line << std::endl;
//      return -1;
//    }
//    N.Assign(f1, f2, f3);
//    faceNormals.at(tOffset + normCount) = N.Normalized();
//    normCount++;
//  }
//  if (meshFile.bad() || meshFile.fail() || normCount != numNormals)
//  {
//    std::cout << "ERROR: while reading normals from mesh file; last line read: " << line << std::endl;
//    return -1;
//  }
//
//  //std::cout << " Mesh Load Complete (Points: " << vertices.size()
//  //  << ", Triangles: " << faces.size() << ")" << std::endl;
//  return 0;
//}
//
//// Save current mesh to a single mesh file
//int cisstMesh::SaveMeshFile(const std::string &filePath)
//{
//  // Save ASCII file in format:
//  //
//  //  POINTS numPoints
//  //  px py pz
//  //   ...
//  //  px py pz
//  //  TRIANGLES numTriangles
//  //  vx1 vx2 vx3
//  //    ...
//  //  vx1 vx2 vx3
//  //  NORMALS numTriangles
//  //  nx ny nz
//  //    ...
//  //  nx ny nz
//  //
//  //  where vx's are indices into the points array
//
//  //std::cout << "Saving mesh to file: " << filePath << std::endl;
//  std::ofstream fs(filePath.c_str());
//  if (!fs.is_open())
//  {
//    std::cout << "ERROR: failed to open file for saving mesh: " << filePath << std::endl;
//    return -1;
//  }
//  fs << "POINTS " << vertices.size() << "\n";
//  for (unsigned int i = 0; i < vertices.size(); i++)
//  {
//    fs << vertices.at(i)[0] << " "
//      << vertices.at(i)[1] << " "
//      << vertices.at(i)[2] << "\n";
//  }
//  fs << "TRIANGLES " << faces.size() << "\n";
//  for (unsigned int i = 0; i < faces.size(); i++)
//  {
//    fs << faces.at(i)[0] << " "
//      << faces.at(i)[1] << " "
//      << faces.at(i)[2] << "\n";
//  }
//  fs << "NORMALS " << faceNormals.size() << "\n";
//  for (unsigned int i = 0; i < faceNormals.size(); i++)
//  {
//    fs << faceNormals.at(i)[0] << " "
//      << faceNormals.at(i)[1] << " "
//      << faceNormals.at(i)[2] << "\n";
//  }
//  fs.close();
//  return 0;
//}
//
//
//// Load mesh from an STL File (.stl)
//int cisstMesh::LoadMeshFromSTLFile(const std::string &stlFilePath)
//{
//  std::cout << "Building mesh from STL file: " << stlFilePath << std::endl;
//  vertices.SetSize(0);
//  faces.SetSize(0);
//  faceNormals.SetSize(0);
//
//  // Assumes STL File has the following format:
//  //
//  //  solid ascii
//  //   facet normal -0.206371 -0.756911 -0.620078
//  //    outer loop
//  //     vertex 171.49 93.7415 0.593199
//  //     vertex 170.849 93.4625 1.14687
//  //     vertex 169.704 93.7397 1.18966
//  //    endloop
//  //   endfacet
//  //    ...
//  //  endsolid
//
//  // open mesh file
//  std::ifstream stlFile;
//  stlFile.open(stlFilePath.c_str());
//  if (!stlFile.is_open())
//  {
//    std::cout << "ERROR: failed to open STL file" << std::endl;
//    return -1;
//  }
//
//  //--- Process STL File ---//
//  //  Initially, store all vertices in a Map in order to ensure uniqueness
//  //  of vertices and to assign sequential indices to each 
//  vct3 v, n;
//  vctInt3 T;
//  float f1, f2, f3;
//  unsigned int itemsRead;
//  std::vector<vctInt3> triangles;         // temp storage
//  std::vector<vct3> triangleNormals;      // temp storage
//  unsigned int vertexCount = 0;
//  unsigned int triangleCount = 0;
//
//  typedef std::map<vct3, unsigned int, VertexCompare> VertexMapType;
//  VertexMapType vertexMap;
//  std::pair<vct3, unsigned int>             mapInput;
//  std::pair<VertexMapType::iterator, bool>  mapResult;
//
//  // read mesh header
//  std::string line;
//  std::getline(stlFile, line);
//  // Note: different STL formats exist, such as
//  //       "solid ascii" or "solid vcg" => just check for "solid"
//  if (line.find("solid") == std::string::npos)
//    //if (line.compare(0,5, "solid") != 0)
//  {
//    std::cout << "ERROR: unrecognized STL format, missing \"solid ascii\" or \"solid vcg\", etc" << std::endl;
//    return -1;
//  }
//  unsigned int update;
//  while (stlFile.good())
//  {
//    // print progress update
//    update = triangleCount % 10000;
//    if (update == 0)
//    {
//      std::cout << "Triangle #: " << triangleCount / 1000 << "k" << std::endl;
//    }
//    // check for end of file (i.e. line: "endsolid")
//    std::getline(stlFile, line);
//    if (line.find("endsolid") != std::string::npos)
//      //if (line.compare(0,8, "endsolid") == 0)
//    { // no more triangles in file
//      break;
//    }
//
//    // read data for next triangle
//    // line: " facet normal -0.206371 -0.756911 -0.620078"
//    //  NOTE: %*[ \n\t] throws away leading white space
//    itemsRead = std::sscanf(line.c_str(), "%*[ \t]facet normal %f %f %f", &f1, &f2, &f3);
//    if (itemsRead != 3)
//    {
//      std::cout << "ERROR: STL file missing \"facet normal %f %f %f\"" << std::endl;
//      return -1;
//    }
//    n.Assign(f1, f2, f3);
//    // line: "  outer loop"
//    std::getline(stlFile, line);
//    if (line.find("outer loop") == std::string::npos)
//      //if (line.compare(0,12, "  outer loop") != 0)
//    {
//      std::cout << "ERROR: STL file missing \"outer loop\"" << std::endl;
//      return -1;
//    }
//    // 3 lines of: "   vertex 171.49 93.7415 0.593199"
//    for (unsigned int i = 0; i < 3; i++)
//    {
//      std::getline(stlFile, line);
//      itemsRead = std::sscanf(line.c_str(), "%*[ \t]vertex %f %f %f", &f1, &f2, &f3);
//      if (itemsRead != 3)
//      {
//        std::cout << "ERROR: STL file missing \"vertex %f %f %f\"" << std::endl;
//        return -1;
//      }
//      v.Assign(f1, f2, f3);
//      // add vertex to map, along with this vertex's index position
//      mapInput.first = v;
//      mapInput.second = vertexCount;
//      mapResult = vertexMap.insert(mapInput);
//      if (mapResult.second == true)
//      { // new vertex was added to map => increment the index position
//        //  for the next vertex
//        vertexCount++;
//        update = vertexCount % 10000;
//        if (update == 0)
//        {
//          std::cout << "Vertex #: " << vertexCount / 1000 << "k" << std::endl;
//        }
//      }
//      // store index position to triangle
//      T[i] = mapResult.first->second;
//      if (mapResult.first->second < 0)
//      {
//        std::cout << "ERROR: vertex maps to negative index" << std::endl;
//        return -1;
//      }
//    }
//    // line: "  end loop"
//    std::getline(stlFile, line);
//    if (line.find("endloop") == std::string::npos)
//      //if (line.compare(0,10, "  endloop") != 0)
//    {
//      std::cout << "ERROR: STL file missing \"  endloop\"" << std::endl;
//      return -1;
//    }
//    // line: " endfacet"
//    std::getline(stlFile, line);
//    //if (line.compare(0,9, " endfacet") != 0)
//    if (line.find("endfacet") == std::string::npos)
//    {
//      std::cout << "ERROR: STL file missing \" endfacet\"" << std::endl;
//      return -1;
//    }
//    // add mesh triangle
//    triangles.push_back(T);
//    triangleNormals.push_back(n);
//    triangleCount++;
//  } // while (stlFile.good())
//  if (!stlFile.good())
//  {
//    std::cout << "ERROR: processing of STL encountered error" << std::endl;
//    return -1;
//  }
//
//  //--- Mesh Post-Processing ---//
//
//  // move triangles to mesh from std::vector
//  faces.SetSize(triangles.size());
//  faceNormals.SetSize(triangles.size());
//  for (unsigned int i = 0; i < triangles.size(); i++)
//  {
//    faces.at(i) = triangles[i];
//    faceNormals.at(i) = triangleNormals[i];
//  }
//  // move vertices to mesh from std::map
//  vertices.SetSize(vertexMap.size());
//  VertexMapType::iterator mapIter;
//  for (mapIter = vertexMap.begin(); mapIter != vertexMap.end(); mapIter++)
//  { // move each vertex to its assigned index position in the mesh array
//    vertices.at(mapIter->second).Assign(mapIter->first);
//  }
//
//  InitializeNoiseModel();
//
//  //std::cout << "Mesh Construction Complete (Points: " << vertices.size()
//  //  << ", Triangles: " << faces.size() << ")" << std::endl;
//
//  return 0;
//}
//
//void cisstMesh::ComputeFaceNormalsFromVertices()
//{
//  int nFaces = faces.size();
//  faceNormals.SetSize(nFaces);
//  vct3 n, v1, v2, v3;
//
//  for (int i = 0; i < nFaces; i++)
//  {
//    // compute face normal by right-hand rule
//    v1 = vertices(faces(i)[0]);
//    v2 = vertices(faces(i)[1]);
//    v3 = vertices(faces(i)[2]);
//    n = vctCrossProduct(v2 - v1, v3 - v1);
//    if (n.Norm() < 1e-12) {
//      std::cout << "ERROR: the computed face normal vector is extremely small; setting normal to all zeros" << std::endl;
//      n.SetAll(0.0);
//    }
//    else {
//      n.NormalizedSelf();
//    }
//    faceNormals(i) = n;
//  }
//}
//
//
//
////--- Legacy I/O ---//
//
//void cisstMesh::ReadMeshFile(const char *fn)
//{
//  char ext[10];
//  unsigned int i, j;
//
//  ResetMesh();
//
//  i = strlen(fn) - 1;
//  while (i >= 0 && fn[i] != '.') i--;
//
//  i++;
//  for (j = i; j < strlen(fn); j++)
//  {
//    ext[j - i] = fn[j];
//  }
//  ext[j - i] = 0;
//
//  if (!strcmp(ext, "sur") || !strcmp(ext, "SUR"))
//  {
//    ReadSURMeshFile(fn);
//  }
//  else if (!strcmp(ext, "sfc") || !strcmp(ext, "SFC"))
//  {
//    ReadSFCMeshFile(fn);
//  }
//  else
//  {
//    std::cout << "ERROR: Read mesh requires .sfc, . SFC, .sur, or .SUR files" << std::endl;
//  };
//
//  InitializeNoiseModel();
//}
//
//void cisstMesh::WriteMeshFile(const char *fn)
//{
//  char ext[10];
//  unsigned int i, j;
//
//  i = strlen(fn) - 1;
//  while (i >= 0 && fn[i] != '.') i--;
//
//  i++;
//  for (j = i; j < strlen(fn); j++)
//  {
//    ext[j - i] = fn[j];
//  }
//  ext[j - i] = 0;
//
//  if (!strcmp(ext, "sur") || !strcmp(ext, "SUR"))
//  {
//    WriteSURMeshFile(fn);
//  }
//  else
//  {
//    std::cout << "SUR Meshes require .sur or .SUR files" << std::endl;
//  };
//}
//
//void cisstMesh::ReadSURMeshFile(const char *fn)
//{
//  FILE *fp;
//
//  if ((fp = fopen(fn, "r")) == NULL)
//  {
//    std::cout << "SUR mesh file not found" << std::endl;
//    return;
//  }
//
//  int vert_n, face_n;
//  float f1, f2, f3;
//
//  fscanf(fp, "%d\n", &vert_n);	// off file
//  vertices.SetSize(vert_n);
//
//  int i;
//  for (i = 0; i < vert_n; i++)
//  {
//    fscanf(fp, "%f %f %f\n", &f1, &f2, &f3);
//    vertices[i] = vct3(f1, f2, f3);
//  }
//
//  fscanf(fp, "%d\n", &face_n);	// sur file
//  faces.SetSize(face_n);
//  faceNeighbors.SetSize(face_n);
//
//  char buff[1024];
//  for (i = 0; i < face_n; i++)
//  {
//    int a, b, c; int d = 1, e = 1, f = 1;
//    //JA have to fix this, probably should be using streams...
//    fgets(buff, 1024, fp);
//    //sscanf(buff,"%d %d %d",&a,&b,&c);
//    sscanf(buff, "%d %d %d %d %d %d\n", &a, &b, &c, &d, &e, &f);
//    //fscanf(fp,"%d %d %d %f %f %f\n",&a,&b,&c,&d,&e,&f);
//    //cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<endl;
//    faces[i].Assign(a, b, c);
//    faceNeighbors[i].Assign(d, e, f);
//  };
//  fclose(fp);
//}
//
//void cisstMesh::ReadSFCMeshFile(const char *fn)
//{
//  FILE *fp;
//
//  if ((fp = fopen(fn, "r")) == NULL)
//  {
//    std::cout << "ERROR: SFC mesh file not found" << std::endl;
//    return;
//  }
//
//  int vert_n, face_n;
//  float f1, f2, f3;
//
//  fscanf(fp, "%d\n", &vert_n);	// off file
//  vertices.SetSize(vert_n);
//
//  int i;
//  for (i = 0; i < vert_n; i++)
//  {
//    fscanf(fp, "%f %f %f\n", &f1, &f2, &f3);
//    vertices[i] = vct3(f1, f2, f3);
//  }
//
//  fscanf(fp, "%d\n", &face_n);	// sur file
//  faces.SetSize(face_n);
//  faceNormals.SetSize(face_n);
//
//
//  for (i = 0; i < face_n; i++)
//  {
//    int a, b, c; int d = 1, e = 1, f = 1;
//    //JA have to fix this, probably should be using streams...
//
//    fscanf(fp, "%d %d %d\n", &a, &b, &c); d = e = f = -1;  // rht hack to get around not having neighbors
//    //fscanf(fp,"%d %d %d %f %f %f\n",&a,&b,&c,&d,&e,&f);
//    //cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<endl;
//    faces[i].Assign(a, b, c);
//    faceNeighbors[i].Assign(d, e, f);
//#if 0
//    printf("%d: ",i); FT.Print(stdout); printf("\n");
//#endif
//  };
//  fclose(fp);
//}
//
//void cisstMesh::WriteSURMeshFile(const char *fn)
//{
//  FILE *fp;
//
//  if ((fp = fopen(fn, "w")) == NULL)
//  {
//    return;
//  }
//
//  int vert_n, face_n;
//
//  vert_n = vertices.size();
//  face_n = faces.size();
//  fprintf(fp, "%d\n", vert_n);											// sur file
//
//  int i;
//
//  for (i = 0; i < vert_n; i++)
//  {
//    vct3 V = vertices[i];
//    // fprintf(fp,"%f %f %f   ; Vertex %d\n",V.x, V.y, V.z, i);;
//    fprintf(fp, "%f %f %f\n", V[0], V[1], V[2]);;
//  };
//
//  fprintf(fp, "%d\n", face_n);										// sur file
//
//  for (i = 0; i < face_n; i++)
//  {
//    vctInt3 face = faces[i];
//    vctInt3 neighbor(-1, -1, -1);
//    if (faceNeighbors.size() > 0) neighbor = faceNeighbors[i];
//    fprintf(fp, "%d %d %d %d %d %d\n",
//      face[0], face[1], face[2],
//      neighbor[0], neighbor[1], neighbor[2]);
//  };
//
//  fclose(fp);
//}
