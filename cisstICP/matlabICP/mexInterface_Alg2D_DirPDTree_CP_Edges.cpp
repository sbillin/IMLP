// ****************************************************************************
//
//  $Id: mexInterface_Alg2D_DirPDTree_CP_Edges.cpp, v0.1 Exp $
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

#include "mexInterface_Alg2D_DirPDTree_CP_Edges.h"

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
  plhs[0] = convertPtr2Mat < mexInterface_Alg2D_DirPDTree_CP_Edges >
    (new mexInterface_Alg2D_DirPDTree_CP_Edges);
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
  destroyObject<mexInterface_Alg2D_DirPDTree_CP_Edges>(prhs[1]);
}

void CommandInitialize(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  edgesV1         ~ first vertices of 2D edges      (in image coords)
  //  edgesV2         ~ second vertices of 2D edges     (in image coords)
  //  edgesNorm       ~ 2D edge normals                 (in image coords)
  //  edgesV1_3D      ~ first vertices of 3D edges      (in camera coords)
  //  edgesV2_3D      ~ second vertices of 3D edges     (in camera coords)
  //

  std::stringstream ss;

  // Get the class instance from the second input
  mexInterface_Alg2D_DirPDTree_CP_Edges& obj =
    *convertMat2Ptr<mexInterface_Alg2D_DirPDTree_CP_Edges>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs != 7)
  {
    MEX_ERROR("Cmd Initialize: requires 0 outputs, 7 inputs.");
  }


  //--- Inputs ---//

  MEX_DEBUG("Parsing Inputs...\n");

  // Parse 2D Edges
  //vctDynamicVector<vct2> edgesV1;
  //vctDynamicVector<vct2> edgesV2;
  //vctDynamicVector<vct2> edgesNorm;

  MEX_DEBUG(" ...edgesV1\n");
  Parse2DVectorArray_Double(prhs[2], obj.edgesV1);
  MEX_DEBUG(" ...edgesV2\n");
  Parse2DVectorArray_Double(prhs[3], obj.edgesV2);
  MEX_DEBUG(" ...edgesNorm\n");
  Parse2DVectorArray_Double(prhs[4], obj.edgesNorm);

  // Parse 3D Edges
  //vctDynamicVector<vct3> edgesV1_3D;
  //vctDynamicVector<vct3> edgesV2_3D;
  //vctDynamicVector<vct3> edgesNorm_3D;

  MEX_DEBUG(" ...edgesV1_3D\n");
  Parse3DVectorArray_Double(prhs[5], obj.edgesV1_3D);
  MEX_DEBUG(" ...edgesV2_3D\n");
  Parse3DVectorArray_Double(prhs[6], obj.edgesV2_3D);

  MEX_DEBUG(" ...done\n");


  //--- Processing ---//

  // build 2D edge covariance tree
  int nThresh = 5;
  double diagThresh = 0.0;
  bool bUseOBB = false;
  obj.pDirTree = new DirPDTree2D_Edges(
    obj.edgesV1, obj.edgesV2, obj.edgesNorm,
    nThresh, diagThresh, bUseOBB);

  // create & initialize algorithm
  if (obj.pAlg)
  {
    delete obj.pAlg;
    obj.pAlg = NULL;
  }
  obj.pAlg = new alg2D_DirPDTree_CP_Edges(obj.pDirTree);
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
  //  MatchPts2D      ~ 2D match points (in image corods)     (2 x numSamps double)  
  //  MatchNorms2D    ~ 2D match norms  (in image corods)     (2 x numSamps double)
  //  MatchPts3D      ~ 3D match points (in camera coords)    (3 x numSamps double)  (optional)
  //  MatchDatums     ~ match datums (numSamps int32)                                (optional)
  //

  // Check parameters
  if ((nlhs < 2 || nlhs > 4) || (nrhs < 4 || nrhs > 5))
  {
    MEX_ERROR("ComputeMatches: requires 2 to 4 outputs, 4 or 5 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_Alg2D_DirPDTree_CP_Edges& obj =
    *convertMat2Ptr<mexInterface_Alg2D_DirPDTree_CP_Edges>(prhs[1]);

  alg2D_DirPDTree_CP_Edges& alg = *obj.pAlg;


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
  obj.matchLambdas.SetSize(obj.nSamples);
  obj.MatchPts_3D.SetSize(obj.nSamples);

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

  // compute 3D match equivalents
  obj.ConvertMatchesTo3D();


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

  // output 3: 3D match points
  if (nlhs >= 3)
  {
    plhs[2] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
    //plhs[1] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
    double *matchPts3D = mxGetPr(plhs[2]);
    //double *matchNorms3D = mxGetPr(plhs[1]);
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*matchPts3D) = obj.MatchPts_3D[i].Element(0);
      matchPts3D++;
      (*matchPts3D) = obj.MatchPts_3D[i].Element(1);
      matchPts3D++;
      (*matchPts3D) = obj.MatchPts_3D[i].Element(2);
      matchPts3D++;

      //(*matchNorms3D) = obj.MatchNorms_3D[i].Element(0);
      //matchNorms3D++;
      //(*matchNorms3D) = obj.MatchNorms_3D[i].Element(1);
      //matchNorms3D++;
      //(*matchNorms3D) = obj.MatchNorms_3D[i].Element(2);
      //matchNorms3D++;
    }
  }

  // output 4: match datums
  if (nlhs >= 4)
  {
    plhs[3] = mxCreateNumericMatrix(1, nSamps, mxINT32_CLASS, mxREAL);
    int *matchDatums = static_cast<int*> (mxGetData(plhs[3]));
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*matchDatums) = obj.matchDatums[i];
      matchDatums++;
    }
  }
}


void mexInterface_Alg2D_DirPDTree_CP_Edges::ComputeMatches()
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

    // save lambda value of match for computing 3D match equivalent
    matchLambdas.Element(s) = pAlg->searchLambdas.Element(matchDatums.Element(s));
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
