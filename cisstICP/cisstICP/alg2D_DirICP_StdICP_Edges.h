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
#ifndef _alg2D_DirICP_StdICP_Edges_h
#define _alg2D_DirICP_StdICP_Edges_h

#include "alg2D_DirICP_StdICP.h"
#include "alg2D_DirPDTree_CP_Edges.h"

class alg2D_DirICP_StdICP_Edges : public alg2D_DirICP_StdICP, alg2D_DirPDTree_CP_Edges
{
  //
  // This class implements the standard ICP algorithm for 
  //  a target shape composed of 2D edges.
  //

  //--- Algorithm Parameters ---//

public:
  
  //DirPDTree2D_Edges *pDirTreeEdges;
  
  // Store match lambdas so that 3D match positions
  //  can be computed following computation of 2D matches:
  //   matchPt = (edgeV2-edgeV1)*matchLambda + edgeV1
  vctDoubleVec    matchLambdas;

  //// temporary buffer storage used to help determine matchLambdas
  ////  (since current architecture cannot store the matchLambdas directly)
  //vctDoubleVec    searchLambdas;


  //--- Algorithm Methods ---//

public:

  // constructor
  alg2D_DirICP_StdICP_Edges(
    DirPDTree2D_Edges *pDirTree,
    vctDynamicVector<vct2> &samplePts,
    vctDynamicVector<vct2> &sampleNorms)
    : alg2D_DirICP_StdICP(pDirTree, samplePts, sampleNorms),
    alg2D_DirPDTree_CP_Edges(pDirTree)
    //pDirTreeEdges(pDirTree)
  {
    SetSamples(samplePts, sampleNorms);
    //matchLambdas.SetSize(samplePts.size());
    //searchLambdas.SetSize(pDirTree->EdgeList.numEdges);
  }

  // destructor
  virtual ~alg2D_DirICP_StdICP_Edges() {}

  void SetSamples(
    vctDynamicVector<vct2> &samplePts,
    vctDynamicVector<vct2> &sampleNorms)
  {
    alg2D_DirICP_StdICP::SetSamples(samplePts, sampleNorms);

    matchLambdas.SetSize(samplePts.size());
  }

  void  SamplePreMatch(unsigned int sampleIndex) {}
  void  SamplePostMatch(unsigned int sampleIndex);

};
#endif
