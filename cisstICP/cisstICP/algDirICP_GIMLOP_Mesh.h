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
#ifndef _algDirICP_GIMLOP_Mesh_h
#define _algDirICP_GIMLOP_Mesh_h

#include "algDirICP_GIMLOP.h"
#include "DirPDTree_Mesh.h"
#include "TriangleClosestPointSolver.h"


class algDirICP_GIMLOP_Mesh : public algDirICP_GIMLOP
{ 
  //
  // This class implements the Kent algorithm for 
  //  a mesh target shape.
  //


  //--- Algorithm Parameters ---//

protected:

  DirPDTree_Mesh *pDirTree;
  TriangleClosestPointSolver TCPS;


  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_GIMLOP_Mesh(
    DirPDTree_Mesh *pDirTree,
    vctDynamicVector<vct3> &samplePts,
    vctDynamicVector<vct3> &sampleNorms,
    const vctDynamicVector<double> &argK,
    const vctDynamicVector<double> &argE,
    const vctDynamicVector<vct3x2> &argL,
    const vctDynamicVector<vct3x3> &M,
    PARAM_EST_TYPE paramEst = PARAMS_FIXED)
    : algDirICP_GIMLOP(pDirTree, samplePts, sampleNorms, argK, argE, argL, M, paramEst),
    pDirTree(pDirTree),
	  TCPS(pDirTree->mesh)
  {};

  // destructor
  virtual ~algDirICP_GIMLOP_Mesh() {}


  //--- PD Tree Interface Methods ---//

  double FindClosestPointOnDatum(
    const vct3 &Xp, const vct3 &Xn,
    vct3 &closest, vct3 &closestNorm,
    int datum);

  int  DatumMightBeCloser(
    const vct3 &Xp, const vct3 &Xn,
    int datum,
    double ErrorBound);

};
#endif
