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
#ifndef _cisstCovTreeBase_h
#define _cisstCovTreeBase_h

#include <stdio.h>
#include <assert.h>
#include <cisstVector.h>
#include <cisstCommon.h>
#include <cisstNumerical.h>

#include "cisstBoundingBox.h"
#include "cisstCovTreeNode.h"
#include "cisstAlgorithmCovTree.h"

//#define DebugCovTree


class cisstCovTreeBase
{
  //
  // This is the base class for a covariance tree.
  //  This class defines the entry point for performing
  //  a search. The type of datum (i.e. triangle, point, etc.)
  //  is unknown to this class and must be defined within a
  //  derived class, making this class abstract.
  // This class stores a pointer to an algorithm object
  //  which implements the key search routines.
  //

  friend class cisstCovTreeNode;


  //--- Variables ---//

public:

  // reference to algorithm must exist here so that all nodes 
  //  may access it
  cisstAlgorithmCovTree *pAlgorithm;

//protected:

#ifdef DebugCovTree
  FILE *debugFile;
  FILE *debugFile2;
#endif

  int NData, NNodes, treeDepth;
  int* DataIndices;
  cisstCovTreeNode *Top;

  ~cisstCovTreeBase(){};


  //--- Methods ---//

public:

  // constructors
  cisstCovTreeBase() :
    NData(0), NNodes(0), treeDepth(0),
    DataIndices(NULL), Top(NULL), pAlgorithm(NULL)
  {
#ifdef DebugCovTree
    debugFile = fopen("debugCovTree.txt","w");
    debugFile2 = fopen("debugCovTree2.txt","w");
#endif
  };

  int FastInitializeProximalDatum(const vct3 &v, vct3 &proxPoint);

  void SetSearchAlgorithm(cisstAlgorithmCovTree *pAlg)
  {
    pAlgorithm = pAlg;
  }

  // Returns the index for the datum in the tree that has lowest match error for
  //  the given point and set the closest point values
  int FindClosestDatum(
    const vct3 &v,
    vct3 &closestPoint,
    int prevDatum,
    double &matchError,
    unsigned int &numNodesSearched);

  int NumData() const { return NData; };
  int NumNodes() const { return NNodes; };
  int TreeDepth() const { return treeDepth; };

  // debug routines
  int   ValidateClosestDatum(const vct3 &v, vct3 &closestPoint);
  int   ValidateClosestDatum_ByEuclideanDist(const vct3 &v, vct3 &closestPoint);
  int   FindTerminalNode(int datum, cisstCovTreeNode **termNode);
  void  PrintTerminalNodes(std::ofstream &fs);


  //--- Virtual Methods ---//
  //
  // These methods require a known datum type
  //
  virtual vct3  DatumSortPoint(int datum) const = 0;
  virtual void  EnlargeBounds(const vctFrm3& F, int datum, cisstBoundingBox& BB) const = 0;


#ifdef ENABLE_NOISE_MODEL
  //--- Noise Model Methods ---//
  virtual vct3x3& DatumCov(int datum) = 0;    // return measurement noise model for this datum
  virtual vct3x3* DatumCovPtr(int datum) = 0; // return measurement noise model for this datum
  virtual vct3& DatumCovEig(int datum) = 0;         
  virtual vct3* DatumCovEigPtr(int datum) = 0;
  // may have to be manually called by user after defining the noise
  //  model of the datums
  //  (depending on the covariance tree type and constructor used)
  void ComputeNodeNoiseModels();
  void ComputeSubNodeNoiseModel(cisstCovTreeNode *node, bool useLocalVarsOverride);
#endif

};


#endif // _cisstCovTreeBase_h