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


#include "cisstPointCloudDir.h"


// Load points and point normals from file
int cisstPointCloudDir::ReadPointsFromFile(
  vctDynamicVector<vct3> &pts,
  vctDynamicVector<vct3> &norms,
  std::string &filePath)
{
  pts.SetSize(0);
  norms.SetSize(0);
  return AppendPointsFromFile(pts, norms, filePath);
}


// Append points and point normals from file
int cisstPointCloudDir::AppendPointsFromFile(
  vctDynamicVector<vct3> &pts,
  vctDynamicVector<vct3> &norms,
  std::string &filePath)
{
  // Text file format:
  //
  //  POINTS numPoints
  //  px py pz nx ny nz
  //   ...
  //  px py pz nx ny nz
  //

  unsigned int itemsRead;
  std::string line;
  float f1, f2, f3;
  float n1, n2, n3;

  unsigned int pOffset;
  pOffset = pts.size();

  std::cout << "Reading pts & normals from file: " << filePath << std::endl;
  std::ifstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filePath << std::endl;
    assert(0);
    return -1;
  }

  // read pts
  unsigned int numPoints;
  std::getline(fs, line);
  itemsRead = std::sscanf(line.c_str(), "POINTS %u", &numPoints);
  if (itemsRead != 1)
  {
    std::cerr << "ERROR: expected POINTS header at line: " << line << std::endl;
    assert(0);
    return -1;
  }
  vct3 v;
  pts.resize(pOffset + numPoints);    // non-destructive
  unsigned int pointCount = 0;
  while (fs.good() && pointCount < numPoints)
  {
    std::getline(fs, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &f1, &f2, &f3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: expected a point value at line: " << line << std::endl;
      assert(0);
      return -1;
    }
    v[0] = f1; v[1] = f2; v[2] = f3;
    pts.at(pOffset + pointCount).Assign(v);
    pointCount++;
  }
  if (fs.bad() || fs.fail() || pointCount != numPoints)
  {
    std::cerr << "ERROR: read points from file failed; last line read: " << line << std::endl;
    assert(0);
    return -1;
  }

  // read normals
  unsigned int numNormals;
  std::getline(fs, line);
  itemsRead = std::sscanf(line.c_str(), "POINT_ORIENTATIONS %u", &numNormals);
  assert(numNormals == numPoints);
  if (itemsRead != 1)
  {
    std::cerr << "ERROR: expected POINT_ORIENTATIONS header at line: " << line << std::endl;
    assert(0);
    return -1;
  }
  vct3 n;
  norms.resize(pOffset + numPoints);  // non-destructive
  unsigned int normCount = 0;
  while (fs.good() && normCount < numNormals)
  {
    std::getline(fs, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f", &n1, &n2, &n3);
    if (itemsRead != 3)
    {
      std::cerr << "ERROR: expected a normal value at line: " << line << std::endl;
      assert(0);
      return -1;
    }
    n[0] = n1; n[1] = n2; n[2] = n3;
    norms.at(pOffset + normCount).Assign(n);
    normCount++;
  }
  if (fs.bad() || fs.fail() || normCount != numNormals)
  {
    std::cerr << "ERROR: read normals from file failed; last line read: " << line << std::endl;
    assert(0);
    return -1;
  }

  fs.close();
  std::cout << " ..." << pts.size() << " sample pts & normals" << std::endl;
  return 0;
}

// Write points and point normals to file
int cisstPointCloudDir::WritePointsToFile(
  vctDynamicVector<vct3> &pts,
  vctDynamicVector<vct3> &norms,
  std::string &filePath)
{
  // Text file format:
  //
  //  POINTS numPoints
  //  px py pz nx ny nz
  //   ...
  //  px py pz nx ny nz
  //

  //std::cout << "Saving pts & normals to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filePath << std::endl;
    assert(0);
    return -1;
  }
  fs << "POINTS " << pts.size() << "\n";
  for (unsigned int i = 0; i < pts.size(); i++)
  {
    fs << pts.at(i)[0] << " " << pts.at(i)[1] << " " << pts.at(i)[2] << "\n";
  }
  fs << "POINT_ORIENTATIONS " << pts.size() << "\n";
  for (unsigned int i = 0; i < pts.size(); i++)
  {
    fs << norms.at(i)[0] << " " << norms.at(i)[1] << " " << norms.at(i)[2] << "\n";
  }
  fs.close();
  return 0;
}