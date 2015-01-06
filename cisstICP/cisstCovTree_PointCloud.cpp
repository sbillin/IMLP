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

#include "cisstCovTree_PointCloud.h"
#include "cisstCovTreeNode.h"

// Build tree from point cloud input
//  this constructor does not define the noise model for the point cloud;
//  the user must do this manually by calling the noise model methods of this class
cisstCovTree_PointCloud::cisstCovTree_PointCloud(vctDynamicVector<vct3> &pointCloud,
  int nThresh, double diagThresh)
{
#ifdef ENABLE_NOISE_MODEL
  avgVarInPlane = 0.0;

  std::cout << std::endl << "===> WARNING: noise model for covariance tree point cloud is not initialized by this constructor"
    << std::endl << std::endl;
#endif

  points = pointCloud;
  NData = pointCloud.size();
  DataIndices = new int[NData];
  for (int i = 0; i < NData; i++)
  {
    DataIndices[i] = i;
  }
  Top = new cisstCovTreeNode(DataIndices, NData, this, NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(nThresh, diagThresh);

#ifdef DebugCovTree
  fprintf(debugFile, "Point Cloud Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
#endif
}

cisstCovTree_PointCloud::~cisstCovTree_PointCloud()
{
  if (Top) delete Top;
  if (DataIndices) delete DataIndices;
}

vct3 cisstCovTree_PointCloud::DatumSortPoint(int datum) const
{ // datum sort point is the point itself
#ifdef DebugCovTree
  return points.at(datum);
#else
  return points.Element(datum);
#endif
}

void cisstCovTree_PointCloud::EnlargeBounds(const vctFrm3& F, int datum, cisstBoundingBox& BB) const
{
  BB.Include(F*points.Element(datum));
}

//void cisstCovTree_PointCloud::PrintDatum(FILE* chan, int level, int datum)
//{
//  for (int i = 0; i < level; i++)
//    fprintf(chan, " ");
//  fprintf(chan, "%5d (%f,%f,%f):", datum, points.at(datum).X(), points.at(datum).Y(), points.at(datum).Z());
//}


#ifdef ENABLE_NOISE_MODEL

// build tree from mesh, converting mesh to point cloud by choosing
//  the centerpoint from each triangle as a point in the cloud;
//  also defines measurement noise over the cloud by defining in-plane variance
//  as average square distance to three vertices and out-of-plane variance 
//  as defined by the argument.
cisstCovTree_PointCloud::cisstCovTree_PointCloud(
  cisstMesh &mesh,
  int nThresh, double diagThresh,
  double noisePerpPlaneSD)
{
  double noiseInPlaneVar, noisePerpPlaneVar;
  double sqrDist, avgSqrDist;
  vct3 c, v0, v1, v2;

  noisePerpPlaneVar = noisePerpPlaneSD*noisePerpPlaneSD;

  // build point cloud from triangle centers
  NData = mesh.NumTriangles();
  DataIndices = new int[NData];
  for (int i = 0; i < NData; i++)
  {
    DataIndices[i] = i;
  }
  points.SetSize(NData);
  for (int i = 0; i < NData; i++)
  {
    points.at(i) = mesh.Triangles.at(i).Midpoint();
  }

  // Define noise properties of point cloud
  //  use triangles to determine in-plane noise
  pointCov.SetSize(NData);
  pointCovEig.SetSize(NData);
  avgSqrDist = 0.0;
  for (int i = 0; i < NData; i++)
  {
    c = mesh.Triangles.at(i).Midpoint();

    // compute in-plane noise
    mesh.VerticesCoords(i, v0, v1, v2);
    sqrDist = (v0 - c).NormSquare();
    sqrDist += (v1 - c).NormSquare();
    sqrDist += (v2 - c).NormSquare();
    sqrDist /= 3.0;
    noiseInPlaneVar = sqrDist;

    // set noise model for this point
    pointCov[i] = ComputeTriangleCovariance(mesh.Triangles.at(i).norm,
      noiseInPlaneVar, noisePerpPlaneVar);
    // must list eigenvalues in descending order
    if (noiseInPlaneVar >= noisePerpPlaneVar)
    {
      pointCovEig[i].Element(0) = noiseInPlaneVar;
      pointCovEig[i].Element(1) = noiseInPlaneVar;
      pointCovEig[i].Element(2) = noisePerpPlaneVar;
    }
    else
    {
      pointCovEig[i].Element(0) = noisePerpPlaneVar;
      pointCovEig[i].Element(1) = noiseInPlaneVar;
      pointCovEig[i].Element(2) = noiseInPlaneVar;
    }
    avgSqrDist += sqrDist;
  }
  avgSqrDist /= NData;
  avgVarInPlane = avgSqrDist;

  std::cout << "Avg VarInPlane: " << avgVarInPlane << std::endl;

  // Build tree nodes
  Top = new cisstCovTreeNode(DataIndices, NData, this, NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(nThresh, diagThresh);

  // Define noise model of the nodes
  ComputeNodeNoiseModels();

#ifdef DebugCovTree
  fprintf(debugFile, "Point Cloud Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
  fprintf(debugFile, " Noise Model:\n  perp-plane variance = %d\n  in-plane variance = %d (avg)\n\n", noisePerpPlaneVar, avgSqrDist);
#endif
}

// build tree from mesh, converting mesh to point cloud by choosing
//  the centerpoint from each triangle as a point in the cloud;
//  also define noise model for points as defined by in-plane and perp-plane
//  standard deviations
cisstCovTree_PointCloud::cisstCovTree_PointCloud(cisstMesh &mesh,
  int nThresh, double diagThresh,
  double noisePerpPlaneSD, double noiseInPlaneSD)
{
  double noiseInPlaneVar, noisePerpPlaneVar;

  noisePerpPlaneVar = noisePerpPlaneSD*noisePerpPlaneSD;
  noiseInPlaneVar = noiseInPlaneSD*noiseInPlaneSD;

  // build point cloud from triangle centers
  NData = mesh.NumTriangles();
  DataIndices = new int[NData];
  for (int i = 0; i < NData; i++)
  {
    DataIndices[i] = i;
  }
  points.SetSize(NData);
  for (int i = 0; i < NData; i++)
  {
    points.at(i) = mesh.Triangles.at(i).Midpoint();
  }

  // Define noise properties of point cloud
  pointCov.SetSize(NData);
  pointCovEig.SetSize(NData);
  // must list eigenvalues in descending order
  if (noiseInPlaneVar >= noisePerpPlaneVar)
  {
    pointCovEig.SetAll(vct3(noiseInPlaneVar, noiseInPlaneVar, noisePerpPlaneVar));
  }
  else
  {
    pointCovEig.SetAll(vct3(noisePerpPlaneVar, noiseInPlaneVar, noiseInPlaneVar));
  }
  for (int i = 0; i < NData; i++)
  {
    // set noise model for this point
    pointCov[i] = ComputeTriangleCovariance(mesh.Triangles.at(i).norm,
      noiseInPlaneVar, noisePerpPlaneVar);
  }
  avgVarInPlane = noiseInPlaneVar;
  //std::cout << "AvgVarInPlane: " << avgVarInPlane << std::endl;

  // Build tree nodes
  Top = new cisstCovTreeNode(DataIndices, NData, this, NULL);
  NNodes = 0; NNodes++;
  treeDepth = Top->ConstructSubtree(nThresh, diagThresh);

  // Define noise model of the nodes
  ComputeNodeNoiseModels();

#ifdef DebugCovTree
  fprintf(debugFile, "Point Cloud Cov Tree built: NNodes=%d  NData=%d  TreeDepth=%d\n", NumNodes(), NumData(), TreeDepth());
  fprintf(debugFile, " Noise Model:\n  perp-plane variance = %d\n  in-plane variance = %d\n\n", noisePerpPlaneVar, noiseInPlaneVar);
#endif
}


// compute a noise covariance matrix having different noise
//  magnitude in-plane vs. out-of-plane for the given plane norm
vct3x3 cisstCovTree_PointCloud::ComputeTriangleCovariance(const vct3 &norm,
  double noiseInPlaneVar,
  double noisePerpPlaneVar)
{
  vct3x3 M, M0;
  vctRot3 R;
  vct3 z(0.0, 0.0, 1.0);

  // set eigenvalues of noise covariance
  //  set plane perpendicular noise component along z-axis
  M0.SetAll(0.0);
  M0.Element(0, 0) = noiseInPlaneVar;
  M0.Element(1, 1) = noiseInPlaneVar;
  M0.Element(2, 2) = noisePerpPlaneVar;

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
  M = R.Transpose()*M0*R;

  return M;
}

void cisstCovTree_PointCloud::SavePointCloudCov(std::string &filePath)
{
  std::cout << "Saving point cloud covariances to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file for saving cov: " << filePath << std::endl;
    assert(0);
  }
  unsigned int numCov = this->pointCov.size();
  //fs << "NUMCOV " << numCov << "\n";
  for (unsigned int i = 0; i < numCov; i++)
  {
    fs << this->pointCov.at(i).Row(0) << " "
      << this->pointCov.at(i).Row(1) << " "
      << this->pointCov.at(i).Row(2) << "\n";
  }
}


#endif  // ENABLE_NOISE_MODEL