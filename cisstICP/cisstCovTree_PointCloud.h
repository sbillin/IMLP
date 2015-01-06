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
#ifndef _cisstCovTree_PointCloud_h
#define _cisstCovTree_PointCloud_h

#include <limits>

#include "cisstCovTreeBase.h"
#include "cisstMesh.h"
#include "cisstPointCloud.h"


class cisstCovTree_PointCloud : public cisstCovTreeBase
{ 
  //
  // This class implements a covariance tree for a point cloud shape
  //


  //--- Variables ---//

public:

  vctDynamicVector<vct3>    points;   // point cloud


  //--- Methods ---//

  // constructor
  //  pointCloud - target shape from which to construct the tree
  //                (each point of the cloud becomes a datum in the tree)
  //  nThresh    - min number of datums to subdivide a node
  //  diagThresh - min physical size to subdivide a node
	cisstCovTree_PointCloud( vctDynamicVector<vct3> &pointCloud,
                          int nThresh, double diagThresh);

  // destructor
  ~cisstCovTree_PointCloud();


  //--- Base Class Virtual Methods ---//

  virtual vct3 DatumSortPoint(int datum) const;  // return sort point of this datum
  virtual void EnlargeBounds(const vctFrm3& F, int datum, cisstBoundingBox& BB) const;


#ifdef ENABLE_NOISE_MODEL

  //--- Noise Model Parameters ---//

  // Noise model
  // these may have to be set manually by user code
  vctDynamicVector<vct3x3>  pointCov;       // covariance of measurement noise
  vctDynamicVector<vct3>    pointCovEig;    // eigenvalues of covariance
  double avgVarInPlane;   // avg square distance between center of mesh triangle and its vertices

  //--- Noise Model Methods ---//

  // constructs point cloud from mesh triangle centers and
  //  defines noise model over the point cloud and nodes
  //  perp-plane noise model component is static
  //  in-plane noise model component defined dynamically over the mesh
  cisstCovTree_PointCloud(
    cisstMesh &mesh,
    int nThresh, double diagThresh,
    double noisePerpPlaneSD);

  // constructs point cloud from mesh triangle centers and
  //  defines noise model over the point cloud and nodes
  //  perp-plane and in-plane noise model components are static
  cisstCovTree_PointCloud(
    cisstMesh &mesh,
    int nThresh, double diagThresh,
    double noisePerpPlaneSD, double noiseInPlaneSD);

  vct3x3 ComputeTriangleCovariance(
    const vct3 &norm,
    double noiseInPlaneSD,
    double noisePerpPlaneSD);

  void SavePointCloud(std::string &file)
  {
    cisstPointCloud::WritePointsToFile(points, file);
  }

  void SavePointCloudCov(std::string &file);

  vct3x3& DatumCov(int datum)         // return measurement noise model for this datum
  {
    return pointCov.Element(datum);
  }
  vct3x3* DatumCovPtr(int datum)      // return measurement noise model for this datum
  {
    return &(pointCov.Element(datum));
  }

  vct3& DatumCovEig(int datum)         // return measurement noise model for this datum
  {
    return pointCovEig.Element(datum);
  }
  vct3* DatumCovEigPtr(int datum)      // return measurement noise model for this datum
  {
    return &(pointCovEig.Element(datum));
  }

#endif

};

#endif