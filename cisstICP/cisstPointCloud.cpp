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

#include "cisstPointCloud.h"


// Load points from file
int cisstPointCloud::ReadPointsFromFile(
  vctDynamicVector<vct3> &pts,
  std::string &filePath)
{
  pts.SetSize(0);
  return AppendPointsFromFile(pts, filePath);
}

// Append points from file
int cisstPointCloud::AppendPointsFromFile(
  vctDynamicVector<vct3> &pts,
  std::string &filePath)
{
  // Text file format:
  //
  //  POINTS numPoints
  //  px py pz
  //   ...
  //  px py pz
  //

  unsigned int itemsRead;
  std::string line;
  float f1, f2, f3;

  unsigned int pOffset;
  pOffset = pts.size();

  std::cout << "Reading pts from file: " << filePath << std::endl;
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
      std::cerr << "ERROR: expeced a point value at line: " << line << std::endl;
      assert(0);
      return -1;
    }
    v[0] = f1;
    v[1] = f2;
    v[2] = f3;
    pts.at(pOffset + pointCount).Assign(v);
    pointCount++;
  }
  if (fs.bad() || fs.fail() || pointCount != numPoints)
  {
    std::cerr << "ERROR: read pts from file failed; last line read: " << line << std::endl;
    assert(0);
    return -1;
  }
  fs.close();

  std::cout << " ..." << pts.size() << " sample pts" << std::endl;

  return 0;
}

// Write points to file
int cisstPointCloud::WritePointsToFile(
  vctDynamicVector<vct3> &pts,
  std::string &filePath)
{
  // Text file format:
  //
  //  POINTS numPoints
  //  px py pz
  //   ...
  //  px py pz
  //
  //  where vx's are indices into the pts array

  std::cout << "Saving pts to file: " << filePath << std::endl;
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
    fs << pts.at(i)[0] << " "
      << pts.at(i)[1] << " "
      << pts.at(i)[2] << "\n";
  }
  fs.close();
  return 0;
}
