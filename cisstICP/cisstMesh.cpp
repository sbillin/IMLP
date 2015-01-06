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
#include <cisstNumerical/nmrLSSolver.h>

#include <assert.h>
#undef NDEBUG       // enable assert in release mode

#include <fstream>

// compares vertex for storing in std::Map
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


// compute a noise covariance matrix for each triangle of mesh having different noise
//  magnitude in-plane vs. out-of-plane
void cisstMesh::ComputeTriangleNoiseModels(
  double noiseInPlaneVar,
  double noisePerpPlaneVar)
{
  vct3x3 M, M0;
  vctRot3 R;
  vct3 z(0.0, 0.0, 1.0);
  vct3 norm;

  // define eigenvalues of noise covariance
  //  set plane perpendicular noise component along z-axis
  M0.SetAll(0.0);
  M0.Element(0, 0) = noiseInPlaneVar;
  M0.Element(1, 1) = noiseInPlaneVar;
  M0.Element(2, 2) = noisePerpPlaneVar;

  // set covariance eigenvalues (in order of decreasing magnitude)
  if (noiseInPlaneVar >= noisePerpPlaneVar)
  {
    TriangleCovEig.SetAll(vct3(noiseInPlaneVar, noiseInPlaneVar, noisePerpPlaneVar));
  }
  else
  {
    TriangleCovEig.SetAll(vct3(noisePerpPlaneVar, noiseInPlaneVar, noiseInPlaneVar));
  }

  for (unsigned int i = 0; i < NumTriangles(); i++)
  {
    norm = TriangleNorm(i);

    // find rotation to align triangle plane normal with the z-axis
    vct3 xProd = vctCrossProduct(norm, z);
    if (xProd.Norm() <= 1e-6)  // protect from divide by zero
    { // norm is already oriented with z-axis
      R = vctRot3::Identity();
    }
    else
    {
      // the norm of the cross product is the same for angles of x deg & x+180 deg
      //  between two vectors => use dot product to determine the angle
      //   NOTE: the angle corresponding to the cross product axis is always > 0;
      //         acos of the dot product gives the correct form
      //   NOTE: the problem with using norm of cross product isn't that we aren't
      //         going the right direction, but rather that we don't rotate far enough
      //         if A & B are seperated by more than 90 degrees.  I.e. if angular
      //         seperation between A & B is 100 degrees, then asin(norm(AxB)) gives
      //         the same angle as if A & B are seperated by 80 degrees => we don't
      //         know if the actual angle is X or 90+X using the norm of cross product.
      vct3 ax = xProd.Normalized();
      double an = acos(vctDotProduct(norm, z));
      //double an = asin(t.Norm());
      vctAxAnRot3 R_AxAn(ax, an);
      R = vctRot3(R_AxAn);
    }

    // compute noise covariance M of this sample and its decomposition:
    //    M = U*S*V'
    // rotate to align normal with z-axis, apply noise covariance, rotate back
    TriangleCov[i] = R.Transpose()*M0*R;
  }
}

void cisstMesh::SaveTriangleCovariances(std::string &filePath)
{
  std::cout << "Saving point cloud covariances to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file for saving cov: " << filePath << std::endl;
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


// called by destructor
void cisstMesh::Reset()
{
  VertexCoordinates.SetSize(0);
  Triangles.SetSize(0);
  TriangleCov.SetSize(0);
  TriangleCovEig.SetSize(0);
}


// Build mesh from a list of vertices, triangles, and normals
int cisstMesh::LoadMesh(
  const vctDynamicVector<vct3> &V,
  const vctDynamicVector<vctInt3> &T,
  const vctDynamicVector<vct3> &N)
{
  //std::cout << "Building mesh from input vectors" << std::endl;

  if (V.size() < 1 || T.size() < 1 || T.size() != N.size())
  {
    std::cerr << "ERROR: invalid input" << std::endl;
    assert(0);
  }

  this->VertexCoordinates = V;

  unsigned int numT = T.size();
  this->Triangles.SetSize(numT);

  cisstTriangle Tri(this);   // provide pointer to this mesh  
  for (unsigned int i = 0; i < numT; i++)
  {
    Tri.Vx[0] = T[i].Element(0);
    Tri.Vx[1] = T[i].Element(1);
    Tri.Vx[2] = T[i].Element(2);
    Tri.norm.Assign(N.Element(i));
    Tri.ComputeBoundingBox();
    this->Triangles.at(i) = Tri;
  }

  std::cout << " Mesh Build Complete (Points: " << this->VertexCoordinates.size()
    << ", Triangles: " << this->Triangles.size() << ")" << std::endl;
  return 0;
}

// Build new mesh from a single mesh file
int cisstMesh::LoadMeshFile(const std::string &meshFilePath)
{
  int rv;
  VertexCoordinates.SetSize(0);
  Triangles.SetSize(0);
  rv = AddMeshFile(meshFilePath);

  TriangleCov.SetSize(Triangles.size());
  TriangleCovEig.SetSize(Triangles.size());
  TriangleCov.SetAll(vct3x3(0.0));
  TriangleCovEig.SetAll(vct3(0.0));
  return rv;
}

// Build new mesh containing data from multiple mesh files
int cisstMesh::LoadMeshFileMultiple(const std::vector<std::string> &meshFilePaths)
{
  int result;
  this->VertexCoordinates.SetSize(0);
  this->Triangles.SetSize(0);
  std::vector<std::string>::const_iterator iter;
  for (iter = meshFilePaths.begin(); iter != meshFilePaths.end(); iter++)
  {
    result = AddMeshFile(*iter);
    if (result == -1) return -1;
  }

  TriangleCov.SetSize(Triangles.size());
  TriangleCovEig.SetSize(Triangles.size());
  TriangleCov.SetAll(vct3x3(0.0));
  TriangleCovEig.SetAll(vct3(0.0));
  return 0;
}

// Load mesh file, adding it to the current mesh while preserving all
//  data currently existing in the mesh
int cisstMesh::AddMeshFile(const std::string &meshFilePath)
{
  // Load mesh from ASCII file having format:
  //
  //  POINTS numPoints
  //  px py pz
  //   ...
  //  px py pz
  //  TRIANGLES numTriangles
  //  vx1 vx2 vx3
  //    ...
  //  vx1 vx2 vx3
  //  NORMALS numTriangles
  //  nx ny nz
  //    ...
  //  nx ny nz
  //
  //  where vx's are indices into the points array

  std::cout << "Loading mesh file: " << meshFilePath << std::endl;

  float f1, f2, f3;
  int d1, d2, d3;
  unsigned int itemsRead;
  std::string line;

  unsigned int vOffset = this->VertexCoordinates.size();
  unsigned int tOffset = this->Triangles.size();

  // open file
  std::ifstream meshFile;
  meshFile.open(meshFilePath.c_str());
  if (!meshFile.is_open())
  {
    std::cerr << "ERROR: failed to open mesh file" << std::endl;
    return -1;
  }

  // read points
  //  (vertex coordinates)
  unsigned int numPoints;
  std::getline(meshFile, line);
  itemsRead = std::sscanf(line.c_str(), "POINTS %u", &numPoints);
  if (itemsRead != 1)
  {
    std::cerr << "ERROR: expected POINTS header at line: " << line << std::endl;
    return -1;
  }
  vct3 v;
  this->VertexCoordinates.resize(vOffset + numPoints);  // non destructive size change
  unsigned int pointCount = 0;
  while (meshFile.good() && pointCount < numPoints)
  {
    std::getline(meshFile, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &f1, &f2, &f3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: expected a point value at line: " << line << std::endl;
      return -1;
    }
    v[0] = f1;
    v[1] = f2;
    v[2] = f3;
    this->VertexCoordinates.at(vOffset + pointCount).Assign(v);
    pointCount++;
  }
  if (meshFile.bad() || meshFile.fail() || pointCount != numPoints)
  {
    std::cerr << "ERROR: read points from mesh file failed; last line read: " << line << std::endl;
    return -1;
  }

  // read triangles
  //  (three values that index the points array, one for each vertex)
  unsigned int numTriangles;
  std::getline(meshFile, line);
  itemsRead = std::sscanf(line.c_str(), "TRIANGLES %u", &numTriangles);
  if (itemsRead != 1)
  {
    std::cerr << "ERROR: expected TRIANGLES header at line: " << line << std::endl;
    return -1;
  }
  cisstTriangle T(this);   // provide pointer to this mesh
  this->Triangles.resize(tOffset + numTriangles); // non destructive size change
  unsigned int triangleCount = 0;
  while (meshFile.good() && triangleCount < numTriangles)
  {
    std::getline(meshFile, line);
    itemsRead = std::sscanf(line.c_str(), "%d %d %d", &d1, &d2, &d3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: expeced three index values on line: " << line << std::endl;
      return -1;
    }
    T.Vx[0] = d1 + vOffset;
    T.Vx[1] = d2 + vOffset;
    T.Vx[2] = d3 + vOffset;
    T.ComputeBoundingBox();
    this->Triangles.at(tOffset + triangleCount) = T;
    triangleCount++;
  }
  if (meshFile.bad() || meshFile.fail() || triangleCount != numTriangles)
  {
    std::cerr << "ERROR: while reading triangles from mesh file; last line read: " << line << std::endl;
    return -1;
  }

  // read triangle normals
  unsigned int numNormals;
  std::getline(meshFile, line);
  itemsRead = std::sscanf(line.c_str(), "NORMALS %u", &numNormals);
  if (numNormals != numTriangles)
  {
    std::cout << std::endl << "ERROR: number of triangles and normals in the mesh do not match" << std::endl << std::endl;
    return -1;
  }
  if (itemsRead != 1)
  {
    std::cout << std::endl << " ===> WARNING: normal vectors are missing in this legacy mesh file <=== " << std::endl << std::endl;
    return 0;
  }
  vct3 tempNorm;
  triangleCount = 0;
  while (meshFile.good() && triangleCount < numTriangles)
  {
    std::getline(meshFile, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &f1, &f2, &f3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: expeced three decimal values on line: " << line << std::endl;
      return -1;
    }
    tempNorm.Assign(f1, f2, f3);
    this->Triangles.at(tOffset + triangleCount).norm.Assign(tempNorm.Normalized());
    triangleCount++;
  }
  if (meshFile.bad() || meshFile.fail() || triangleCount != numTriangles)
  {
    std::cerr << "ERROR: while reading normals from mesh file; last line read: " << line << std::endl;
    return -1;
  }

  std::cout << " Mesh Load Complete (Points: " << this->VertexCoordinates.size()
    << ", Triangles: " << this->Triangles.size() << ")" << std::endl;
  return 0;
}

// Save current mesh to a single mesh file
int cisstMesh::SaveMeshFile(const std::string &filePath)
{
  // Save ASCII file in format:
  //
  //  POINTS numPoints
  //  px py pz
  //   ...
  //  px py pz
  //  TRIANGLES numTriangles
  //  vx1 vx2 vx3
  //    ...
  //  vx1 vx2 vx3
  //  NORMALS numTriangles
  //  nx ny nz
  //    ...
  //  nx ny nz
  //
  //  where vx's are indices into the points array

  std::cout << "Saving mesh to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file for saving mesh: " << filePath << std::endl;
    return -1;
  }
  fs << "POINTS " << this->VertexCoordinates.size() << "\n";
  for (unsigned int i = 0; i < this->VertexCoordinates.size(); i++)
  {
    fs << this->VertexCoordinates.at(i)[0] << " "
      << this->VertexCoordinates.at(i)[1] << " "
      << this->VertexCoordinates.at(i)[2] << "\n";
  }
  fs << "TRIANGLES " << this->Triangles.size() << "\n";
  for (unsigned int i = 0; i < this->Triangles.size(); i++)
  {
    fs << this->Triangles.at(i).Vx[0] << " "
      << this->Triangles.at(i).Vx[1] << " "
      << this->Triangles.at(i).Vx[2] << "\n";
  }
  fs << "NORMALS " << this->Triangles.size() << "\n";
  for (unsigned int i = 0; i < this->Triangles.size(); i++)
  {
    fs << this->Triangles.at(i).norm[0] << " "
      << this->Triangles.at(i).norm[1] << " "
      << this->Triangles.at(i).norm[2] << "\n";
  }
  fs.close();
  return 0;
}


// Load mesh from an STL File (.stl)
int cisstMesh::LoadMeshFromSTLFile(const std::string &stlFilePath)
{
  std::cout << "Building mesh from STL file: " << stlFilePath << std::endl;
  this->VertexCoordinates.SetSize(0);
  this->Triangles.SetSize(0);

  // Assumes STL File has the following format:
  //
  //  solid ascii
  //   facet normal -0.206371 -0.756911 -0.620078
  //    outer loop
  //     vertex 171.49 93.7415 0.593199
  //     vertex 170.849 93.4625 1.14687
  //     vertex 169.704 93.7397 1.18966
  //    endloop
  //   endfacet
  //    ...
  //  endsolid

  // open mesh file
  std::ifstream stlFile;
  stlFile.open(stlFilePath.c_str());
  if (!stlFile.is_open())
  {
    std::cerr << "ERROR: failed to open STL file" << std::endl;
    return -1;
  }

  //--- Process STL File ---//
  //  Initially, store all vertices in a Map in order to ensure uniqueness
  //  of vertices and to assign sequential indices to each 
  vct3 v;
  float f1, f2, f3;
  unsigned int itemsRead;
  cisstTriangle newT(this);    // provide pointer to this mesh
  std::vector<cisstTriangle> triangles; // for temp storage
  unsigned int vertexCount = 0;
  unsigned int triangleCount = 0;

  typedef std::map<vct3, unsigned int, VertexCompare>   VertexMapType;
  VertexMapType vertexMap;
  std::pair<vct3, unsigned int>             mapInput;
  std::pair<VertexMapType::iterator, bool>  mapResult;

  // read mesh header
  std::string line;
  std::getline(stlFile, line);
  // Note: different STL formats exist, such as
  //       "solid ascii" or "solid vcg" => just check for "solid"
  if (line.find("solid") == std::string::npos)
    //if (line.compare(0,5, "solid") != 0)
  {
    std::cerr << "ERROR: unrecognized STL format, missing \"solid ascii\" or \"solid vcg\", etc" << std::endl;
    return -1;
  }
  unsigned int update;
  while (stlFile.good())
  {
    // print progress update
    update = triangleCount % 10000;
    if (update == 0)
    {
      std::cout << "Triangle #: " << triangleCount / 1000 << "k" << std::endl;
    }
    // check for end of file (i.e. line: "endsolid")
    std::getline(stlFile, line);
    if (line.find("endsolid") != std::string::npos)
      //if (line.compare(0,8, "endsolid") == 0)
    { // no more triangles in file
      break;
    }

    // read data for next triangle
    // line: " facet normal -0.206371 -0.756911 -0.620078"
    //  NOTE: %*[ \n\t] throws away leading white space
    itemsRead = std::sscanf(line.c_str(), "%*[ \t]facet normal %f %f %f", &f1, &f2, &f3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: STL file missing \"facet normal %f %f %f\"" << std::endl;
      return -1;
    }
    // store normal to triangle
    newT.norm.Assign(f1, f2, f3);
    // line: "  outer loop"
    std::getline(stlFile, line);
    if (line.find("outer loop") == std::string::npos)
      //if (line.compare(0,12, "  outer loop") != 0)
    {
      std::cerr << "ERROR: STL file missing \"outer loop\"" << std::endl;
      return -1;
    }
    // 3 lines of: "   vertex 171.49 93.7415 0.593199"
    for (unsigned int i = 0; i < 3; i++)
    {
      std::getline(stlFile, line);
      itemsRead = std::sscanf(line.c_str(), "%*[ \t]vertex %f %f %f", &f1, &f2, &f3);
      if (itemsRead != 3)
      {
        std::cerr << "ERROR: STL file missing \"vertex %f %f %f\"" << std::endl;
        return -1;
      }
      v.Element(0) = f1;
      v.Element(1) = f2;
      v.Element(2) = f3;
      // add vertex to map, along with this vertex's index position
      mapInput.first = v;
      mapInput.second = vertexCount;
      mapResult = vertexMap.insert(mapInput);
      if (mapResult.second == true)
      { // new vertex was added to map => increment the index position
        //  for the next vertex
        vertexCount++;
        update = vertexCount % 10000;
        if (update == 0)
        {
          std::cout << "Vertex #: " << vertexCount / 1000 << "k" << std::endl;
        }
      }
      // store index position to triangle
      newT.Vx[i] = mapResult.first->second;
      if (mapResult.first->second < 0)
      {
        std::cerr << "ERROR: vertex maps to negative index" << std::endl;
        return -1;
      }
    }
    // line: "  end loop"
    std::getline(stlFile, line);
    if (line.find("endloop") == std::string::npos)
      //if (line.compare(0,10, "  endloop") != 0)
    {
      std::cerr << "ERROR: STL file missing \"  endloop\"" << std::endl;
      return -1;
    }
    // line: " endfacet"
    std::getline(stlFile, line);
    //if (line.compare(0,9, " endfacet") != 0)
    if (line.find("endfacet") == std::string::npos)
    {
      std::cerr << "ERROR: STL file missing \" endfacet\"" << std::endl;
      return -1;
    }
    // add mesh triangle
    triangles.push_back(newT);
    triangleCount++;
  } // while (stlFile.good())
  if (!stlFile.good())
  {
    std::cerr << "ERROR: processing of STL encountered error" << std::endl;
    return -1;
  }

  //--- Mesh Post-Processing ---//

  // move triangles to mesh from std::vector
  this->Triangles.SetSize(triangles.size());
  unsigned int i = 0;
  std::vector<cisstTriangle>::iterator iterT;
  for (iterT = triangles.begin(); iterT != triangles.end(); iterT++)
  {
    this->Triangles.at(i) = (*iterT);
    i++;
  }
  // move vertices to mesh from std::map
  this->VertexCoordinates.SetSize(vertexMap.size());
  VertexMapType::iterator mapIter;
  for (mapIter = vertexMap.begin(); mapIter != vertexMap.end(); mapIter++)
  { // move each vertex to its assigned index position in the mesh array
    this->VertexCoordinates.at(mapIter->second).Assign(mapIter->first);
  }

  // compute bounding boxes of triangles, now that vertex array is populated
  for (unsigned int i = 0; i < this->Triangles.size(); i++)
  {
    this->Triangles.at(i).ComputeBoundingBox();
  }

  //// test
  //std::cout << std::endl << "Looking for triangle norm reversals..." << std::endl;
  //vct3 v1, v2, v3, n;
  //for (unsigned int i = 0; i<this->NumTriangles(); i++)
  //{
  //  v1 = this->VertexCoords(this->TriangleVertex(i,0));
  //  v2 = this->VertexCoords(this->TriangleVertex(i,1));
  //  v3 = this->VertexCoords(this->TriangleVertex(i,2));
  //  n.CrossProductOf(v2-v1,v3-v1);
  //  if (n.DotProduct(this->TriangleNorm(i)) < 0)
  //  { 
  //    std::cout << " ---> Norm reversed for triangle: " << i << std::endl;
  //  }
  //}

  std::cout << "Mesh Construction Complete (Points: " << this->VertexCoordinates.size()
    << ", Triangles: " << this->Triangles.size() << ")" << std::endl;

  return 0;
}

//// Create mesh from a Legacy VTK File (.vtk)
//// This method needs update to compute normal vectors.
////  It has not been updated because MATLAB script has been written which
////  is more convenient for this purpose.
//int cisstMesh::LoadMeshFromLegacyVTKFile( std::string &vtkFilePath )
//{
//  std::cout << "Building mesh from VTK file: " << vtkFilePath << std::endl;
//  this->Triangles.SetSize(0);
//  this->VertexCoordinates.SetSize(0);
//
//  // open vtk file
//  std::ifstream vtkFile;
//  vtkFile.open(vtkFilePath.c_str());
//  if (!vtkFile.is_open())
//  {
//    std::cerr << "ERROR: failed to open VTK file" << std::endl;
//    return -1;
//  }
//
//  //--- Process VTK File ---//
//  //  The vtk file format already has each point stored uniquely
//  //  => no need for Map as done for STL format
//  vct3 v;
//  float f1, f2, f3, f4, f5, f6, f7, f8, f9;
//  unsigned int itemsRead;
//  unsigned int u1;
//
//  // read header
//  std::string line;
//  std::getline( vtkFile,line );
//  if (line.compare(0,22, "# vtk DataFile Version") != 0)
//  {
//    std::cerr << "ERROR: expected Legacy VTK file version at line: " << line << std::endl;
//    return -1;
//  }
//  std::getline( vtkFile, line );
//  std::getline( vtkFile, line );
//  if (line.compare(0,5, "ASCII") != 0)
//  {
//    std::cerr << "ERROR: expected ASCII format at line: " << line << std::endl;
//    return -1;
//  }
//  std::getline( vtkFile, line );
//  if (line.compare(0,16, "DATASET POLYDATA") != 0)
//  {
//    std::cerr << "ERROR: expected POLYDATA data type at line: " << line << std::endl;
//    return -1;
//  }
//
//  // read points
//  unsigned int numPoints;
//  std::getline( vtkFile, line );
//  itemsRead = std::sscanf ( line.c_str(), "POINTS %lu float", &numPoints);
//  if (itemsRead != 1)
//  {
//    std::cerr << "ERROR: expected POINTS header at line: " << line << std::endl;
//    return -1;
//  }
//  //  Note: up to 3 sets of points (9 coord values) may exist on the same line
//  this->VertexCoordinates.SetSize(numPoints);
//  unsigned int pointCount = 0;
//  while ( vtkFile.good() && pointCount < numPoints )
//  {
//    std::getline( vtkFile, line );
//    itemsRead = std::sscanf ( line.c_str(), "%f %f %f %f %f %f %f %f %f", 
//      &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9);
//    if (itemsRead <= 0 || itemsRead % 3 != 0)
//    {
//      std::cerr << "ERROR: \"" << itemsRead << 
//        "\" items read is not a multiple of 3 in line: " << line << std::endl;
//      return -1;
//    }
//    if (itemsRead >= 3)
//    {
//      v.Element(0) = f1;
//      v.Element(1) = f2;
//      v.Element(2) = f3;
//      this->VertexCoordinates.at(pointCount).Assign(v);
//      pointCount++;
//      if (pointCount >= numPoints) break;
//    }
//    if (itemsRead >= 6)
//    {
//      v.Element(0) = f4;
//      v.Element(1) = f5;
//      v.Element(2) = f6;
//      this->VertexCoordinates.at(pointCount).Assign(v);
//      pointCount++;
//      if (pointCount >= numPoints) break;
//    }
//    if (itemsRead == 9)
//    {
//      v.Element(0) = f7;
//      v.Element(1) = f8;
//      v.Element(2) = f9;
//      this->VertexCoordinates.at(pointCount).Assign(v);
//      pointCount++;
//      if (pointCount >= numPoints) break;
//    }
//  }
//  if (!vtkFile.good() || pointCount != numPoints)
//  {
//      std::cerr << "ERROR: read points from VTK file failed; last line read: " << line << std::endl;
//      return -1;
//  }
//
//  //read triangles
//  unsigned int numStrips;
//  std::getline( vtkFile, line );
//  itemsRead = std::sscanf ( line.c_str(), "TRIANGLE_STRIPS %u %u", &numStrips, &u1);
//  if (itemsRead != 2)
//  {
//    std::cerr << "ERROR: expected TRIANGLE_STRIPS header at line: " << line << std::endl;
//    return -1;
//  }
//  //  Note: any number of triangles may exit on the same line
//  //  Line Format: n v1 v2 v3 v4 ... vn
//  //    where n = number of vertices on line (index values for the point array)
//  //    each ordered set (vi, vi+1, vi+2) forms a complete triangle
//  //    => n-2 triangles exist on each line
//  // due to complexity, read each line as a stream
//  std::vector<cisstTriangle> triangles;   // temp storage
//  cisstTriangle Tprev(this), T(this);   // provide reference to this mesh
//  std::stringstream ss;
//  unsigned int numIndices;
//  unsigned int stripCount = 0;
//  while ( vtkFile.good() && stripCount < numStrips )
//  {
//    std::getline( vtkFile, line );
//    ss.str( line );
//    ss >> numIndices;
//    // first three values form first triangle
//    ss >> Tprev.Vx[0] >> Tprev.Vx[1] >> Tprev.Vx[2];
//    triangles.push_back(Tprev);
//
//    // subsequent index values in this triangle strip each complete a new triangle
//    //  start at 3, since 3 vertices already read
//    for (unsigned int i = 3; i < numIndices; i++)
//    {
//      // last 2 indices of prev triangle form first 2 indices of this triangle
//      T.Vx[0] = Tprev.Vx[1];
//      T.Vx[1] = Tprev.Vx[2];
//      ss >> T.Vx[2];
//      triangles.push_back(T);
//      Tprev = T;
//    }
//    stripCount++;
//  }
//  if (!vtkFile.good() || stripCount != numStrips)
//  {
//      std::cerr << "ERROR: read triangle strips from VTK file failed; last line read: " << line << std::endl;
//      return -1;
//  }
//
//  //--- Mesh Post-Processing ---//
//
//  // move triangles to mesh from std::vector
//  this->Triangles.SetSize(triangles.size());
//  unsigned int i = 0;
//  std::vector<cisstTriangle>::iterator iterT;
//  for (iterT = triangles.begin(); iterT != triangles.end(); iterT++)
//  {
//    Triangles.at(i) = (*iterT);
//  }
//
//  std::cout << "Mesh Construction Complete (Points: " << this->VertexCoordinates.size() 
//    << ", Triangles: " << this->Triangles.size() << ")" << std::endl;
//  return 0;
//}



//--- Legacy I/O ---//

void cisstMesh::ReadMeshFile(const char *fn)
{
  char ext[10];
  unsigned int i, j;

  i = strlen(fn) - 1;
  while (i >= 0 && fn[i] != '.') i--;

  i++;
  for (j = i; j < strlen(fn); j++)
  {
    ext[j - i] = fn[j];
  }
  ext[j - i] = 0;

  if (!strcmp(ext, "sur") || !strcmp(ext, "SUR"))
  {
    this->Reset();
    ReadSURMeshFile(fn);
  }
  else if (!strcmp(ext, "sfc") || !strcmp(ext, "SFC"))
  {
    this->Reset();
    ReadSFCMeshFile(fn);
  }
  else
  {
    cisstError("Read mesh requires .sfc, . SFC, .sur, or .SUR files");
  };

  return;
}

void cisstMesh::WriteMeshFile(const char *fn)
{
  char ext[10];
  unsigned int i, j;

  i = strlen(fn) - 1;
  while (i >= 0 && fn[i] != '.') i--;

  i++;
  for (j = i; j < strlen(fn); j++)
  {
    ext[j - i] = fn[j];
  }
  ext[j - i] = 0;

  if (!strcmp(ext, "sur") || !strcmp(ext, "SUR"))
  {
    WriteSURMeshFile(fn);
  }
  else
  {
    cisstError("SUR Meshes require .sur or .SUR files");
  };

  return;
}

void cisstMesh::ReadSURMeshFile(const char *fn)
{
  FILE *fp;

  if ((fp = fopen(fn, "r")) == NULL)
  {
    cisstError("SUR mesh file not found");
    return;
  }
  /*
    ifstream ifs(fn);
    if(!ifs) return;
    */

  int vert_n, face_n;
  float f1, f2, f3;

  // JA 3/2000 silence unused var warnings
  //Vec3 v0,v1,v2;
  //double x_b=0,y_b=0,z_b=0;

  fscanf(fp, "%d\n", &vert_n);	// off file
  VertexCoordinates.SetSize(vert_n);

  int i;
  for (i = 0; i < vert_n; i++)
  {
    fscanf(fp, "%f %f %f\n", &f1, &f2, &f3);
    VertexCoordinates[i] = vct3(f1, f2, f3);
  }

  fscanf(fp, "%d\n", &face_n);	// sur file
  Triangles.SetSize(face_n);

  char buff[1024];
  for (i = 0; i < face_n; i++)
  {
    int a, b, c; int d = 1, e = 1, f = 1;
    //JA have to fix this, probably should be using streams...
    fgets(buff, 1024, fp);
    //sscanf(buff,"%d %d %d",&a,&b,&c);
    sscanf(buff, "%d %d %d %d %d %d\n", &a, &b, &c, &d, &e, &f);
    //fscanf(fp,"%d %d %d %f %f %f\n",&a,&b,&c,&d,&e,&f);
    //cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<endl;
    //JA change pointer to reference
    cisstTriangle& FT = Triangles[i];
    FT.VertexIndex(0) = a; FT.NeighborIndex(0) = d;
    FT.VertexIndex(1) = b; FT.NeighborIndex(1) = e;
    FT.VertexIndex(2) = c; FT.NeighborIndex(2) = f;
    FT.myMesh = this;
    FT.ComputeBoundingBox();
  };
  fclose(fp);

  return;
}

void cisstMesh::ReadSFCMeshFile(const char *fn)
{
  FILE *fp;

  if ((fp = fopen(fn, "r")) == NULL)
  {
    cisstError("SFC mesh file not found");
    return;
  }
  /*
   ifstream ifs(fn);
   if(!ifs) return;
   */

  int vert_n, face_n;
  float f1, f2, f3;

  // JA 3/2000 silence unused var warnings
  //Vec3 v0,v1,v2;
  //double x_b=0,y_b=0,z_b=0;

  fscanf(fp, "%d\n", &vert_n);	// off file
  VertexCoordinates.SetSize(vert_n);

  int i;
  for (i = 0; i < vert_n; i++)
  {
    fscanf(fp, "%f %f %f\n", &f1, &f2, &f3);
    VertexCoordinates[i] = vct3(f1, f2, f3);
  }

  fscanf(fp, "%d\n", &face_n);	// sur file
  Triangles.SetSize(face_n);


  for (i = 0; i < face_n; i++)
  {
    int a, b, c; int d = 1, e = 1, f = 1;
    //JA have to fix this, probably should be using streams...

    fscanf(fp, "%d %d %d\n", &a, &b, &c); d = e = f = -1;  // rht hack to get around not having neighbors
    //fscanf(fp,"%d %d %d %f %f %f\n",&a,&b,&c,&d,&e,&f);
    //cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<endl;
    //JA change pointer to reference
    cisstTriangle& FT = Triangles[i];
    FT.VertexIndex(0) = a; FT.NeighborIndex(0) = d;
    FT.VertexIndex(1) = b; FT.NeighborIndex(1) = e;
    FT.VertexIndex(2) = c; FT.NeighborIndex(2) = f;
    FT.myMesh = this;
    FT.ComputeBoundingBox();
#if 0
    printf("%d: ",i); FT.Print(stdout); printf("\n");
#endif
  };
  fclose(fp);

  return;
}

void cisstMesh::WriteSURMeshFile(const char *fn)
{
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL)
  {
    return;
  }

  int vert_n, face_n;

  vert_n = VertexCoordinates.size();
  face_n = Triangles.size();
  fprintf(fp, "%d\n", vert_n);											// sur file

  int i;

  for (i = 0; i < vert_n; i++)
  {
    vct3 V = VertexCoordinates[i];
    // fprintf(fp,"%f %f %f   ; Vertex %d\n",V.x, V.y, V.z, i);;
    fprintf(fp, "%f %f %f\n", V[0], V[1], V[2]);;
  };

  fprintf(fp, "%d\n", face_n);										// sur file

  for (i = 0; i < face_n; i++)
  {
    cisstTriangle* T = &Triangles[i];
    fprintf(fp, "%d %d %d %d %d %d\n",
      T->VertexIndex(0), T->VertexIndex(1), T->VertexIndex(2),
      T->NeighborIndex(0), T->NeighborIndex(1), T->NeighborIndex(2));
  };

  fclose(fp);

  return;
}

//void cisstMesh::Print(FILE* fp)
//{
//  int vert_n, face_n;
//
//  vert_n = VertexCoordinates.size();
//  face_n = Triangles.size();
//  fprintf(fp, "%d\n", vert_n);											// sur file
//
//  int i;
//
//  for (i = 0; i < vert_n; i++)
//  {
//    vct3 V = VertexCoordinates[i];
//    // fprintf(fp,"%f %f %f   ; Vertex %d\n",V.x, V.y, V.z, i);;
//    fprintf(fp, "%5d: %f %f %f\n", i, V[0], V[1], V[2]);;
//  };
//
//  fprintf(fp, "%d\n", face_n);										// sur file
//
//  for (i = 0; i < face_n; i++)
//  {
//    cisstTriangle* T = &Triangles[i];
//    fprintf(fp, "%5d (%4d %4d %4d %5d %5d %5d)", i,
//      T->VertexIndex(0), T->VertexIndex(1), T->VertexIndex(2),
//      T->NeighborIndex(0), T->NeighborIndex(1), T->NeighborIndex(2));
//    for (int vx = 0; vx < 3; vx++)
//    {
//      fprintfVct3Bracketed(fp, VertexCoord(T->VertexIndex(vx)));
//    };
//    fprintf(fp, "\n");
//  };
//
//}
