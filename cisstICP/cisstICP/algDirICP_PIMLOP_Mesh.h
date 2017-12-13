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
#ifndef _algDirICP_PIMLOP_Mesh_h
#define _algDirICP_PIMLOP_Mesh_h

#include "algDirICP_PIMLOP.h"
#include "algDirPDTree_vonMisesPrj_Mesh.h"

class algDirICP_PIMLOP_Mesh 
  : public algDirICP_PIMLOP, public algDirPDTree_vonMisesPrj_Mesh
{ 
  //
  // This class implements the base class algorithm for a mesh target shape.
  //

  //--- Algorithm Parameters ---//


  //--- Algorithm Methods ---//

public:

  // constructor
  algDirICP_PIMLOP_Mesh(
    DirPDTree_Mesh *pDirTree,
    const vctDynamicVector<vct3> &samplePts,
    const vctDynamicVector<vct2> &sampleNorms2d,
    const vctDynamicVector<vctRot3> &Rx_pln,
    const vctDynamicVector<double> &sample_k,
    const vctDynamicVector<vct3x3> &sample_M)
    : algDirICP_PIMLOP(pDirTree, samplePts, sampleNorms2d, Rx_pln, sample_k, sample_M),
    algDirPDTree_vonMisesPrj_Mesh(pDirTree)
  {}
 
  // destructor
  virtual ~algDirICP_PIMLOP_Mesh() {}


  // standard ICP algorithm virtual methods
  void SamplePreMatch(unsigned int sampleIndex);
};


#endif
