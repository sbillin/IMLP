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

#include "DirPDTreeBase.h"

#include <stdio.h>
#include <limits>
#include <fstream>

#include <cisstVector.h>
#include <cisstCommon.h>
//#include <cisstNumerical/nmrLSSolver.h>

// quickly find an approximate initial match by dropping straight down the
//   tree to the node containing the sample point and picking a datum from there
int DirPDTreeBase::FastInitializeProximalDatum(
  const vct3 &v, const vct3 &n, vct3 &proxPoint, vct3 &proxNorm)
{
  // find proximal leaf node
  DirPDTreeNode *pNode;
  pNode = Top;

  while (!pNode->IsTerminalNode())
  {
    pNode = pNode->GetChildSplitNode(v);
  }

  int proxDatum = pNode->Datum(0);        // choose any datum from the leaf node
  proxPoint = DatumSortPoint(proxDatum);  // choose any point on the datum
  proxNorm = DatumNorm(proxDatum);

  return proxDatum;
}

double DirPDTreeBase::ComputeDatumMatchError( const vct3 &v, const vct3 &n, int datum )
{
  vct3 tmp1, tmp2;
  return pAlgorithm->FindClosestPointOnDatum(v, n, tmp1, tmp2, datum);
}

// Return the index for the datum in the tree that is closest to the given point
//  in terms of the complete error (orientation + distance) and set the closest point
int DirPDTreeBase::FindClosestDatum(
  const vct3 &v, const vct3 &n,
  vct3 &closestPoint, vct3 &closestPointNorm,
  int prevDatum,
  double &matchError,
  unsigned int &numNodesSearched,
  double currentMatchError)
{
  unsigned int numNodesVisited = 0;
  numNodesSearched = 0;

  // check if previous match was feasible
  // SDB: by specifying a good starting datum (such as previous closest datum)
  //      the search for new closest datum is more efficient because
  //      the bounds value will be a good initial guess => fewer datums are
  //      closely searched.
  if (prevDatum >= 0)
  { // previous match was feasible, make sure it is still feasible
    //  set error bound to a big number so that this call returns
    //  false only if the datum is infeasible rather than if the datum
    //  cannot be closer
    // TODO: refactor this
    if (!pAlgorithm->DatumMightBeCloser(v, n, prevDatum, 1e100))
    { // previous match is now infeasible
      prevDatum = -1;
      matchError = currentMatchError;
    }
    else
    {
      double tmpMatchError = pAlgorithm->FindClosestPointOnDatum(v, n, closestPoint, closestPointNorm, prevDatum);
      // update the new current matchError
      matchError = tmpMatchError < currentMatchError ? tmpMatchError : currentMatchError;
    }
  }
  else
  { // previous match was infeasible
    matchError = currentMatchError;
  }

  // since all datums must lie within the root node, we don't need to do a node bounds
  // check on the root => it is more efficient to explicitly search each child node of the
  // root rather than searching the root node itself.
  int datum;
  if (treeDepth > 0)
  {
    // 1st call to pLEq updates both distance bound and closest point
    //  before 2nd call to pMore. If pMore returns (-1), then pMore had
    //  nothing better than pLEq and the resulting datum should be
    //  the return value of pLEq (whether that is -1 or a closer datum index)
    int ClosestLEq = -1;
    int ClosestMore = -1;
    ClosestLEq = Top->pLEq->FindClosestDatum(v, n, closestPoint, closestPointNorm, matchError, numNodesVisited, numNodesSearched);
    ClosestMore = Top->pMore->FindClosestDatum(v, n, closestPoint, closestPointNorm, matchError, numNodesVisited, numNodesSearched);
    datum = (ClosestMore < 0) ? ClosestLEq : ClosestMore;
  }
  else
  { // if there is only one node, we must start from the root
    datum = Top->FindClosestDatum(v, n, closestPoint, closestPointNorm, matchError, numNodesVisited, numNodesSearched);
  }
  if (datum < 0)
  {
    datum = prevDatum;  // no feasible datum found closer than previous datum
  }

  //std::cout << "numNodesVisited: " << numNodesVisited << "\tnumNodesSearched: " << numNodesSearched << std::endl;
  return datum;
}

// Exhaustive linear search of all datums in the tree for validation of closest datum
int DirPDTreeBase::ValidateClosestDatum(
  const vct3 &v, const vct3 &n,
  vct3 &closestPoint, vct3 &closestPointNorm)
{
  double bestError = std::numeric_limits<double>::max();
  int    bestDatum = -1;
  double error;
  vct3 datumPoint;
  vct3 datumNorm;
  for (int datum = 0; datum < NData; datum++)
  {
    error = pAlgorithm->FindClosestPointOnDatum(v, n, datumPoint, datumNorm, datum);
    if (error < bestError)
    {
      bestError = error;
      bestDatum = datum;
      closestPoint = datumPoint;
      closestPointNorm = datumNorm;
    }
  }
  return bestDatum;
}

int DirPDTreeBase::FindTerminalNode(int datum, DirPDTreeNode **termNode)
{
  return Top->FindTerminalNode(datum, termNode);
}

void DirPDTreeBase::PrintTerminalNodes(std::ofstream &fs)
{
  Top->PrintTerminalNodes(fs);
}

