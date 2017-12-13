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
#ifndef _algDirICP_IMLOP_Mesh_h
#define _algDirICP_IMLOP_Mesh_h

#include "algDirICP_IMLOP.h"
#include "DirPDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"


class algDirICP_IMLOP_Mesh : public algDirICP_IMLOP
{
  //
  // This class implements the vMFG algorithm for 
  //  a mesh target shape.
  //


  //--- Algorithm Parameters ---//

protected:

  DirPDTree_Mesh *pDirTree;
  TriangleClosestPointSolver TCPS;


  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_IMLOP_Mesh(
    DirPDTree_Mesh *pDirTree,
    vctDynamicVector<vct3> &samplePts,
    vctDynamicVector<vct3> &sampleNorms,
    double kinit = 0.0, double sigma2init = 1.0, double wRpos = 0.5,
    double kfactor = 1.0,
    bool dynamicallyEstParams = true) :
      algDirICP_IMLOP(pDirTree, samplePts, sampleNorms, kinit, sigma2init, wRpos, kfactor, dynamicallyEstParams),
      pDirTree(pDirTree),
      TCPS(pDirTree->mesh)
  {}

  // destructor
  virtual ~algDirICP_IMLOP_Mesh() {}


  //--- PD Tree Interface Methods ---//

  double FindClosestPointOnDatum(
    const vct3 &v, const vct3 &n,
    vct3 &closest, vct3 &closestNorm,
    int datum);

  int  DatumMightBeCloser(
    const vct3 &v, const vct3 &n,
    int datum,
    double ErrorBound);

};

#endif
