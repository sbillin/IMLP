// ****************************************************************************
//
//  $Id: mexInterface_Alg2D_DirPDTree_CP_PointCloud.cpp, v0.1 Exp $
//
//  Description:
//	        
//
//  Usage:
//	        
//
//  See Also:
//	        
//
//  Author(s):  Seth Billings
//
//  Created on:
//
//
//              Developed by the Engineering Research Center for
//          Computer-Integrated Surgical Systems & Technology (cisstST)
//
//               Copyright(c) 2001, The Johns Hopkins University
//                          All Rights Reserved.
//
//  Experimental Software - Not certified for clinical use.
//  For use only by cisstST sponsored research projects.  For use outside of
//  cisstST, please contact Dr. Russell Taylor, cisstST Director, at rht@cs.jhu.edu
//  
// ****************************************************************************

#include "mexInterface_Alg2D_DirPDTree_CP_PointCloud.h"

#include "cisstPointCloud2D.h"

// debug
//#define ValidatePDTreeSearch

#ifdef ValidatePDTreeSearch
std::ofstream validFS("./debugPDTreeSearchValidation.txt");
vct2 validPoint;
vct2 validNorm;
int  validDatum;
double validDist;
double validAng;
double validError;
double searchError;
unsigned int numValidDatums;
unsigned int numInvalidDatums;
double validPercent;
double doubleEps = 1e-16;
int validIter = 0;
#endif

// NOTES:
//
// mexGet: returns value of the specified property in the specified graphics
//         object on success. Returns NULL on failure. 
//         Do not modify the return argument from mexGet. Changing the data in 
//         a const (read-only) mxArray can produce undesired side effects.
//

// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //mexPrintf("I'm Running!\n");

  // Get the command string
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
  {
    MEX_ERROR("First input should be a command string less than 64 characters long.");
  }

  //--- Standard Functions ---//

  // Command: New
  if (!strcmp("new", cmd))
  {
    return CommandNew(nlhs, plhs, nrhs, prhs);
  }

  // All Other Comamnds:
  //  2nd input must be class instance handle
  if (nrhs < 2)
  {
    MEX_ERROR("ERROR: Second input should be a class instance handle.");
  }

  // Command: Delete
  if (!strcmp("delete", cmd))
  {
    return CommandDelete(nlhs, plhs, nrhs, prhs);
  }


  //--- ICP Methods ---//

  // Command: Initialize 
  if (!strcmp("Initialize", cmd))
  {
    return CommandInitialize(nlhs, plhs, nrhs, prhs);
  }

  // Command: Match
  if (!strcmp("ComputeMatches", cmd))
  {
    return CommandComputeMatches(nlhs, plhs, nrhs, prhs);
  }

  MEX_ERROR("ERROR: Command not recognized.");
}


void CommandNew(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ // create new class instance

  // Expected Input:
  //  cmd

  // Check parameters
  if (nlhs != 1 || nrhs != 1)
  {
    MEX_ERROR("New: requires 1 outputs, 1 inputs.");
  }
  // Return a handle to a new C++ instance
  plhs[0] = convertPtr2Mat < mexInterface_Alg2D_DirPDTree_CP_PointCloud >
    (new mexInterface_Alg2D_DirPDTree_CP_PointCloud);
}


void CommandDelete(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ // delete class instance

  // Expected Input:
  //  cmd
  //  class handle

  // Warn if other commands were ignored
  if (nlhs != 0 || nrhs != 2)
  {
    MEX_WARNING("Delete: requires 0 outputs, 2 inputs.");
  }

  // Destroy the C++ object
  destroyObject<mexInterface_Alg2D_DirPDTree_CP_PointCloud>(prhs[1]);
}

void CommandInitialize(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  points         ~ 2D sample point locations
  //

  std::stringstream ss;

  // Get the class instance from the second input
  mexInterface_Alg2D_DirPDTree_CP_PointCloud& obj =
    *convertMat2Ptr<mexInterface_Alg2D_DirPDTree_CP_PointCloud>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs != 3)
  {
    MEX_ERROR("Cmd Initialize: requires 0 outputs, 3 inputs.");
  }


  //--- Inputs ---//

  MEX_DEBUG("Parsing Inputs...\n");

  // Parse 2D Points
  MEX_DEBUG(" ...points\n");
  Parse2DVectorArray_Double(prhs[2], obj.points);

  MEX_DEBUG(" ...done\n");


  //--- Processing ---//

  // build 2D edge covariance tree
  int nThresh = 5;
  double diagThresh = 0.0;
  bool bUseOBB = false;
  obj.pDirTree = new DirPDTree2D_PointCloud(
    cisstPointCloud2D(obj.points),
    nThresh, diagThresh, bUseOBB);

  // create & initialize algorithm
  if (obj.pAlg)
  {
    delete obj.pAlg;
    obj.pAlg = NULL;
  }
  obj.pAlg = new alg2D_DirPDTree_CP_PointCloud(obj.pDirTree);
}

void CommandComputeMatches(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  samplePtsXfmd       ~ 2D sample points
  //  sampleNormsXfmd     ~ 2D sample orientations
  //  matchDatumsInit     ~ initial match datums [optional]
  // 
  // Output:
  //  MatchPts2D      ~ 2D match points (in image coords)     (2 x numSamps double)  
  //  MatchNorms2D    ~ 2D match norms  (in image coords)     (2 x numSamps double)
  //  MatchDatums     ~ match datums (numSamps int32)                                (optional)
  //

  // Check parameters
  if ((nlhs < 2 || nlhs > 3) || (nrhs < 4 || nrhs > 5))
  {
    MEX_ERROR("ComputeMatches: requires 2 to 3 outputs, 4 or 5 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_Alg2D_DirPDTree_CP_PointCloud& obj =
    *convertMat2Ptr<mexInterface_Alg2D_DirPDTree_CP_PointCloud>(prhs[1]);

  alg2D_DirPDTree_CP_PointCloud& alg = *obj.pAlg;


  //--- Inputs ---//

  MEX_DEBUG("Parsing Inputs...\n");

  // Parse Samples
  //vctDynamicVector<vct2>    samplePts;
  //vctDynamicVector<vct2>    sampleNorms;
  MEX_DEBUG(" ...samplePts\n");
  Parse2DVectorArray_Double(prhs[2], obj.samplePtsXfmd);
  MEX_DEBUG(" ...sampleNorms\n");
  Parse2DVectorArray_Double(prhs[3], obj.sampleNormsXfmd);

  // Parse Match Datum Initializers
  bool bDatumsInitialized = false;
  if (nrhs >= 5)
  {
    MEX_DEBUG(" ...matchDatumsInit\n");
    ParseVector_Int(prhs[4], obj.matchDatums);
    bDatumsInitialized = true;
  }

  MEX_DEBUG(" ...done\n");


  //--- Computation ---//

  obj.nSamples = obj.samplePtsXfmd.size();

  obj.matchPts.SetSize(obj.nSamples);
  obj.matchNorms.SetSize(obj.nSamples);
  obj.matchDatums.SetSize(obj.nSamples);
  obj.matchErrors.SetSize(obj.nSamples);

  if (!bDatumsInitialized)
  {
    // initialize match datums with accelerated approximate search
    for (unsigned int i = 0; i < obj.nSamples; i++)
    {
      obj.matchDatums[i] = obj.pDirTree->FastInitializeProximalDatum(
        obj.samplePtsXfmd[i], obj.sampleNormsXfmd[i],
        obj.matchPts[i], obj.matchNorms[i]);
    }
  }

  // compute matches
  obj.ComputeMatches();


  //--- Outputs ---//

  unsigned int nSamps = obj.nSamples;

  // outputs 1 & 2: 2D match points and normals
  plhs[0] = mxCreateDoubleMatrix(2, nSamps, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(2, nSamps, mxREAL);
  double *matchPts = mxGetPr(plhs[0]);
  double *matchNorms = mxGetPr(plhs[1]);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    (*matchPts) = obj.matchPts[i].Element(0);
    matchPts++;
    (*matchPts) = obj.matchPts[i].Element(1);
    matchPts++;

    (*matchNorms) = obj.matchNorms[i].Element(0);
    matchNorms++;
    (*matchNorms) = obj.matchNorms[i].Element(1);
    matchNorms++;
  }

  // output 3: match datums
  if (nlhs >= 3)
  {
    plhs[2] = mxCreateNumericMatrix(1, nSamps, mxINT32_CLASS, mxREAL);
    int *matchDatums = static_cast<int*> (mxGetData(plhs[2]));
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*matchDatums) = obj.matchDatums[i];
      matchDatums++;
    }
  }
}


void mexInterface_Alg2D_DirPDTree_CP_PointCloud::ComputeMatches()
{
  // Find the point on the model having lowest match error
  //  for each sample point

#ifdef ValidatePDTreeSearch  
  numInvalidDatums = 0;
  numValidDatums = 0;
#endif

  unsigned int nSamples = samplePtsXfmd.size();

  unsigned int nodesSearched = 0;
  minNodesSearched = std::numeric_limits<unsigned int>::max();
  maxNodesSearched = std::numeric_limits<unsigned int>::min();
  avgNodesSearched = 0;

  for (unsigned int s = 0; s < nSamples; s++)
  {
    // inform algorithm beginning new match
    //SamplePreMatch(s);

    // Find best match for this sample
    matchDatums.Element(s) = pDirTree->FindClosestDatum(
      samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
      matchPts.Element(s), matchNorms.Element(s),
      matchDatums.Element(s),
      matchErrors.Element(s),
      nodesSearched);

    avgNodesSearched += nodesSearched;
    minNodesSearched = (nodesSearched < minNodesSearched) ? nodesSearched : minNodesSearched;
    maxNodesSearched = (nodesSearched > maxNodesSearched) ? nodesSearched : maxNodesSearched;

#ifdef ValidatePDTreeSearch        
    validDatum = pDirTree->ValidateClosestDatum(
      samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
      validPoint, validNorm);
    if (validDatum != matchDatums.Element(s))
    {
      // It is possible to have different datums for same point if the match
      //  lies on a datum edge; if this is the case, then the search didn't
      //  actually fail
      // Note: cannot compare two double values for exact equality due to
      //       inexact representation of decimal values in binary arithmetic
      searchError = (validPoint - matchPts.Element(s)).NormSquare();
      if (searchError > doubleEps)
      {
        numInvalidDatums++;
        //matchDist = (matchPts.Element(s)-samplePtsXfmd.Element(s)).Norm();
        //matchAng  = acos( vctDotProduct(matchNorms.Element(s),sampleNormsXfmd.Element(s)) );
        validDist = (validPoint - samplePtsXfmd.Element(s)).Norm();
        validAng = acos(vctDotProduct(validNorm, sampleNormsXfmd.Element(s)));
        vct2 tmp1, tmp2;
        //searchError = algorithm->FindClosestPointOnDatum(
        //                    samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
        //                    tmp1, tmp2, matchDatums.Element(s));
        validError = pAlg->FindClosestPointOnDatum(
          samplePtsXfmd.Element(s), sampleNormsXfmd.Element(s),
          tmp1, tmp2, validDatum);
        validFS << "Correspondence Search Data:  (foundMatch / validMatch)" << std::endl
          << " MatchError = " << matchErrors.Element(s) << "/" << validError << std::endl
          //<< " dPos = " << dist_PostMatch.Element(s) << "/" << validDist << std::endl
          //<< " dAng = " << dTheta*180.0 / cmnPI << "/" << validAng*180.0 / cmnPI << std::endl
          << " XfmSamplePoint = [" << samplePtsXfmd.Element(s) << "]" << std::endl
          << " XfmSampleNorm  = [" << sampleNormsXfmd.Element(s) << "]" << std::endl
          << " MatchPoint =     [" << matchPts.Element(s) << "]" << std::endl
          << " MatchNorm =      [" << matchNorms.Element(s) << "]" << std::endl
          << " MatchDatum = " << matchDatums.Element(s) << std::endl
          << " ValidPoint =     [" << validPoint << "]" << std::endl
          << " ValidNorm =      [" << validNorm << "]" << std::endl
          << " ValidDatum = " << validDatum << std::endl
          << " SampleIndex = " << s << std::endl;

        //validFS << "Fact = [" << std::endl << Fact << "]" << std::endl;

        DirPDTree2DNode *termNode = 0;
        pDirTree->FindTerminalNode(validDatum, &termNode);
        if (!termNode)
        {
          std::cout << "ERROR: did not find terminal node for datum: " << validDatum << std::endl;
          assert(0);
        }
        validFS << " Valid Terminal Node:" << std::endl;
        validFS << "   MinCorner: " << termNode->Bounds.MinCorner << std::endl;
        validFS << "   MaxCorner: " << termNode->Bounds.MaxCorner << std::endl;
        validFS << "   NData: " << termNode->NData << std::endl;
        validFS << "Fnode_valid = [" << std::endl << termNode->F << "]" << std::endl;
      }
      else
      {
        numValidDatums++;
      }
    }
    else
    {
      numValidDatums++;
    }
#endif

    //SamplePostMatch(s);
  }

  avgNodesSearched /= nSamples;

#ifdef ValidatePDTreeSearch  
  validPercent = (double)numValidDatums / (double)nSamples * 100.0;
  validFS << "iter " << validIter << ":  NumMatches(valid/invalid): "
    << numValidDatums << "/" << numInvalidDatums << "  valid% = "
    << validPercent << std::endl;
  validIter++;
#endif

}
