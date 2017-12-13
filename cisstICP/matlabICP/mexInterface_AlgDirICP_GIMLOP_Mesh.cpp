// ****************************************************************************
//
//  $Id: mexInterface_AlgDirICP_GIMLOP_Mesh.cpp, v0.1 Exp $
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

#include "mexInterface_AlgDirICP_GIMLOP_Mesh.h"


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

  // debug:
  std::string cmdStr(cmd);
  cmdStr = "mexFunction: cmd = " + cmdStr + "\n";
  MEX_DEBUG(cmdStr.c_str());

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

  // Command: Initialize 
  if (!strcmp("Initialize", cmd))
  {
    return CommandInitialize(nlhs, plhs, nrhs, prhs);
  }

  // Command: Set Samples
  if (!strcmp("SetSamples", cmd))
  {
    return CommandSetSamples(nlhs, plhs, nrhs, prhs);
  }


  //--- ICP Methods ---//

  // Command: Initialize 
  if (!strcmp("ICP_InitializeParameters", cmd))
  {
    return CommandICP_InitializeParameters(nlhs, plhs, nrhs, prhs);
  }

  // Command: Match
  if (!strcmp("ICP_ComputeMatches", cmd))
  {
    return CommandICP_ComputeMatches(nlhs, plhs, nrhs, prhs);
  }

  // Command: Register Matches
  if (!strcmp("ICP_RegisterMatches", cmd))
  {
    return CommandICP_RegisterMatches(nlhs, plhs, nrhs, prhs);
  }

  // Command: Evaluate Error Function
  if (!strcmp("ICP_EvaluateErrorFunction", cmd))
  {
    return CommandICP_EvaluateErrorFunction(nlhs, plhs, nrhs, prhs);
  }

  // Command: Termination Test
  if (!strcmp("ICP_Terminate", cmd))
  {
    return CommandICP_Terminate(nlhs, plhs, nrhs, prhs);
  }


  //--- Debug Methods ---//


  // Command: Set Matches
  if (!strcmp("Debug_SetMatches", cmd))
  {
    return CommandDebug_SetMatches(nlhs, plhs, nrhs, prhs);
  }

  MEX_ERROR("ERROR: Command not recognized.");

  MEX_ERROR("ERROR: Command not recognized.");


  //// Command: Test 
  //if (!strcmp("Test", cmd))
  //{
  //  return CommandTest(  nlhs, plhs, nrhs, prhs );
  //}
  //// Command: Test2
  //if (!strcmp("Test2", cmd))
  //{
  //  return CommandTest2(  nlhs, plhs, nrhs, prhs );
  //}

  //// Command: Template 
  //if (!strcmp("template", cmd))
  //{
  //  return CommandInitialize(  nlhs, plhs, nrhs, prhs );
  //}
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
  plhs[0] = convertPtr2Mat < mexInterface_AlgDirICP_GIMLOP_Mesh >
    (new mexInterface_AlgDirICP_GIMLOP_Mesh);
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
  destroyObject<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);
}

// Initialize the algorithm class with registration data
void CommandInitialize(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  mesh
  //    V ~ 3D vertex positions     (3 x Nv double)
  //    T ~ triangle vertex indices (3 x Nt integer)
  //    N ~ 3D triangle normals     (3 x Nt double)
  //  sample points                 (3 x Ns double)
  //  sample normals                (3 x Ns double)
  //  M ~ covariance matrices of sample points    (3 x 3 x Ns double)
  //  k ~ concentration of sample orientations    (Ns x 1 double)
  //  E ~ eccentricity of sample orientations     (Ns x 1 double)
  //  L ~ major / minor axis of Kent distribution (2 x 3 x Ns double)
  //      (for sample orientations)
  //  dynamicParamEst   (integer)
  //
  std::stringstream ss;

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs < 12)
  {
    MEX_ERROR("Initialize: requires 0 outputs, 12 inputs.");
  }

  // Parse Mesh  
  vctDynamicVector<vct3> V;
  vctDynamicVector<vctInt3> T;
  vctDynamicVector<vct3> N;

  MEX_DEBUG("Parsing Inputs...\n");
  MEX_DEBUG(" ...V\n");
  Parse3DVectorArray_Double(prhs[2], V);
  MEX_DEBUG(" ...T\n");
  Parse3DVectorArray_Int(prhs[3], T);
  MEX_DEBUG(" ...N\n");
  Parse3DVectorArray_Double(prhs[4], N);

  // Parse sample data
  vctDynamicVector<vct3> samplePts;
  vctDynamicVector<vct3> sampleNorms;
  vctDynamicVector<vct3x3>  M;
  vctDynamicVector<double>  k;
  vctDynamicVector<double>  E;
  vctDynamicVector<vct3x2>  L;
  int dynamicParamEst;

  MEX_DEBUG(" ...SamplePts\n");
  if (mxGetN(prhs[5]) > 0)  // TODO: change these to if nrhs ...
  {
    Parse3DVectorArray_Double(prhs[5], samplePts);
  }
  MEX_DEBUG(" ...SampleNorms\n");
  if (mxGetN(prhs[6]) > 0)
  {
    Parse3DVectorArray_Double(prhs[6], sampleNorms);
  }
  MEX_DEBUG(" ...sample covariances (M)\n");
  if (mxGetN(prhs[7]) > 0)
  {
    ParseMatrixArray3x3_Double(prhs[7], M);
  }
  MEX_DEBUG(" ...sample orientation concentrations (k)\n");
  if (mxGetN(prhs[8]) > 0)
  {
    ParseVector_Double(prhs[8], k);
  }
  MEX_DEBUG(" ...sample orientation eccentricities (E)\n");
  if (mxGetN(prhs[9]) > 0)
  {
    ParseVector_Double(prhs[9], E);
  }
  MEX_DEBUG(" ...Kent major/minor axis (L)\n");
  if (mxGetN(prhs[10]) > 0)
  {
    ParseMatrixArray_KentL(prhs[10], L);
    //ParseMatrixArray3x2_Double(prhs[10], L);
    //ParseMatrixArray2x3_Double(prhs[10], L);
  }
  MEX_DEBUG(" ...dynamic parameter estimation (bDynamicParamEst)\n");
  if (mxGetN(prhs[11]) > 0)
  {
    ParseInteger(prhs[11], dynamicParamEst);
  }

  MEX_DEBUG(" ...done\n");

  //vctDynamicVector<double>  B(E.size());
  //for (unsigned int i = 0; i < E.size(); i++)
  //{
  //  B(i) = E(i) * k(i) / 2.0;
  //}
  algDirICP_GIMLOP::PARAM_EST_TYPE dynParamEstType;
  switch (dynamicParamEst)
  {
  case 0:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_FIXED;
    break;
  case 1:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_DYN_THRESH;
    break;
  case 2:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_FIXED_1stITER_POS;
    break;
  default:
    MEX_ERROR("ERROR: Dynamic Parameter Estimation Type Unrecognized!");
  }

  // build mesh
  MEX_DEBUG("Building mesh...\n");
  cisstMesh *pMesh = new cisstMesh();
  pMesh->LoadMesh(&V, &T, &N);
  if (pMesh->NumVertices() == 0)
  {
    MEX_ERROR("ERROR: Build mesh failed\n");
  }
  // build covariance tree
  int    nThresh = 20;       // Cov Tree Params
  double diagThresh = 20.0;  //  ''
  //int    nThresh = 5;       // Cov Tree Params
  //double diagThresh = 5.0;  //  ''
  MEX_DEBUG("Building mesh covariance tree...\n");
  DirPDTree_Mesh *pPDTree = new DirPDTree_Mesh(*pMesh, nThresh, diagThresh);
  ss.str("");
  ss << "Tree built:  NNodes = " << pPDTree->NumNodes() << "  NDatums = "
    << pPDTree->NumData() << "  TreeDepth = " << pPDTree->TreeDepth() << "\n";
  MEX_DEBUG(ss.str());

  // create & initialize algorithm object
  if (obj.pAlg)
  {
    delete obj.pAlg;
    obj.pAlg = NULL;
  }
  obj.pAlg = new algDirICP_GIMLOP_Mesh(
    pPDTree, samplePts, sampleNorms, k, E, L, M, dynParamEstType);
}

// Set the algorithm samples
void CommandSetSamples(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  sample points                 (3 x Ns double)
  //  sample normals                (3 x Ns double)
  //  M ~ covariance matrices of sample points    (3 x 3 x Ns double)
  //  k ~ concentration of sample orientations    (Ns x 1 double)
  //  E ~ eccentricity of sample orientations     (Ns x 1 double)
  //  L ~ major / minor axis of Kent distribution (2 x 3 x Ns double)
  //      (for sample orientations)
  //  dynamicParamEst   (int)
  //  

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs != 9)
  {
    MEX_ERROR("SetSamples: requires 0 outputs, 9 inputs.");
  }

  // Parse sample data
  vctDynamicVector<vct3> samplePts;
  vctDynamicVector<vct3> sampleNorms;
  vctDynamicVector<vct3x3>  M;
  vctDynamicVector<double>  k;
  vctDynamicVector<double>  E;  
  vctDynamicVector<vct3x2>  L;
  int dynamicParamEst;

  MEX_DEBUG(" ...SamplePts\n");
  Parse3DVectorArray_Double(prhs[2], samplePts);
  MEX_DEBUG(" ...SampleNorms2d\n");
  Parse3DVectorArray_Double(prhs[3], sampleNorms);
  MEX_DEBUG(" ...sample covariances (M)\n");
  ParseMatrixArray3x3_Double(prhs[4], M);
  MEX_DEBUG(" ...sample orientation concentrations (k)\n");
  ParseVector_Double(prhs[5], k);
  MEX_DEBUG(" ...sample orientation concentrations (E)\n");
  ParseVector_Double(prhs[6], E);
  MEX_DEBUG(" ...Kent major/minor axis\n");
  ParseMatrixArray_KentL(prhs[7], L);
  MEX_DEBUG(" ...dynamic parameter estimation (bDynamicParamEst)\n");
  ParseInteger(prhs[8], dynamicParamEst);

  MEX_DEBUG(" ...done\n");


  algDirICP_GIMLOP::PARAM_EST_TYPE dynParamEstType;
  switch (dynamicParamEst)
  {
  case 0:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_FIXED;
    break;
  case 1:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_DYN_THRESH;
    break;
  case 2:
    dynParamEstType = algDirICP_GIMLOP::PARAMS_FIXED_1stITER_POS;
    break;
  default:
    MEX_ERROR("ERROR: Dynamic Parameter Estimation Type Unrecognized!");
  }

  // set samples
  obj.pAlg->SetSamples(samplePts, sampleNorms, k, E, L, M, dynParamEstType);
}

// initialize algorithm parameters for the ICP loop
void CommandICP_InitializeParameters(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  Freg ~ [R,t]      (3 x 4 double)
  //

  // get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  // check parameters
  if (nlhs != 0 || nrhs != 3)
  {
    MEX_ERROR("Initialize: requires 0 outputs, 3 inputs.");
  }

  // parse transform
  vctFrm3 FGuess;
  ParseTransform(prhs[2], FGuess);

  // debug:
  std::stringstream ss;
  ss << "CommandICP_InitializeParameters (FGuess) = " << std::endl
    << FGuess << std::endl;
  MEX_DEBUG(ss.str());

  // initialize algorithm parameters
  obj.pAlg->ICP_InitializeParameters(FGuess);
}

void CommandICP_ComputeMatches(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  // 
  // Output:
  //  SampleMatchPts    ~ matches for 3D sample pts      (3 x numPts double)  [optional]
  //  SampleMatchNorms  ~ matches for 3D sample norms    (3 x numPts double)  [optional]
  //  

  // Check parameters
  if (nlhs > 2 || nrhs != 2)
  {
    MEX_ERROR("ComputeMatches: requires 0-2 outputs, 2 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  algDirICP_GIMLOP_Mesh& alg = *obj.pAlg;


  //--- Inputs ---//


  //--- Processing ---//

  // compute matches
  alg.ICP_ComputeMatches();

  // update post-match algorithm parameters
  alg.ICP_UpdateParameters_PostMatch();

  // filter matches
  alg.ICP_FilterMatches();


  //--- Outputs ---//

  unsigned int nSamps = alg.nSamples;

  if (nlhs >= 1)
  {
    // output 1: matchPts
    plhs[0] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
    double *matchPts = mxGetPr(plhs[0]);
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*matchPts) = alg.matchPts[i].Element(0);
      matchPts++;
      (*matchPts) = alg.matchPts[i].Element(1);
      matchPts++;
      (*matchPts) = alg.matchPts[i].Element(2);
      matchPts++;
    }
  }

  if (nlhs >= 2)
  {
    // output 2: matchNorms
    plhs[1] = mxCreateDoubleMatrix(3, nSamps, mxREAL);
    double *dataNorms = mxGetPr(plhs[1]);
    for (unsigned int i = 0; i < nSamps; i++)
    {
      (*dataNorms) = alg.matchNorms[i].Element(0);
      dataNorms++;
      (*dataNorms) = alg.matchNorms[i].Element(1);
      dataNorms++;
      (*dataNorms) = alg.matchNorms[i].Element(2);
      dataNorms++;
    }
  }
}


void CommandICP_RegisterMatches(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  // 
  // Output:
  //  Freg ~ [R,t]      (3 x 4 double)    (Y = Freg * X)
  //  sigma2    ~  updated average trace of covariances
  //  k         ~  updated average orientation concentration
  //  e         ~  updated average eccentricity

  std::stringstream ss;

  // Check parameters
  if (nlhs < 1 || nlhs > 4 || nrhs != 2)
  {
    MEX_ERROR("ComputeMatches: requires 1-4 outputs, 2 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  algDirICP_GIMLOP_Mesh& alg = *obj.pAlg;


  //--- Inputs ---//


  //--- Processing ---//

  // compute matches
  vctFrm3 Freg = alg.ICP_RegisterMatches();

  // update post-register algorithm parameters
  alg.ICP_UpdateParameters_PostRegister(Freg);

  //--- Outputs---//

  // output 1: transform
  plhs[0] = mxCreateDoubleMatrix(3, 4, mxREAL);
  double *Fout = mxGetPr(plhs[0]);
  // write out rotation
  vctRot3 R = Freg.Rotation();
  for (unsigned int i = 0; i < 3; i++)
  {
    (*Fout) = R.Column(i).Element(0);
    Fout++;
    (*Fout) = R.Column(i).Element(1);
    Fout++;
    (*Fout) = R.Column(i).Element(2);
    Fout++;
  }
  // write out translation
  vct3 t = Freg.Translation();
  (*Fout) = t.Element(0);
  Fout++;
  (*Fout) = t.Element(1);
  Fout++;
  (*Fout) = t.Element(2);
  Fout++;
    
  // output 2-3: noise model status
  if (nlhs > 1)
  {
    plhs[1] = mxCreateDoubleScalar(alg.meanSigma2);
    plhs[2] = mxCreateDoubleScalar(alg.meanK);
    plhs[3] = mxCreateDoubleScalar(alg.meanE);
  }

  // debug:  
  ss.str("");
  ss << "CommandICP_RegisterMatches (Freg) = " << std::endl
    << Freg << std::endl;
  MEX_DEBUG(ss.str());

}

void CommandICP_EvaluateErrorFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  // 
  // Output:
  //  fval  ~ value of error function (double)
  //  

  // Check parameters
  if (nlhs != 1 || nrhs != 2)
  {
    MEX_ERROR("ComputeMatches: requires 1 outputs, 2 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  algDirICP_GIMLOP_Mesh& alg = *obj.pAlg;


  //--- Inputs ---//


  //--- Processing ---//

  // evaluate error function
  double fval = alg.ICP_EvaluateErrorFunction();


  //--- Outputs---//

  // output 1: fval
  plhs[0] = mxCreateDoubleScalar(fval);

}


void CommandICP_Terminate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  // 
  // Output:
  //  bTerminate  ~ boolean specifying algorithm-specific termination (logical scalar)
  //  

  // Check parameters
  if (nlhs != 1 || nrhs != 2)
  {
    MEX_ERROR("ComputeMatches: requires 1 outputs, 2 inputs.");
  }

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  algDirICP_GIMLOP_Mesh& alg = *obj.pAlg;


  //--- Inputs ---//


  //--- Processing ---//

  // evaluate error function
  bool bTerminate = alg.ICP_Terminate( alg.Freg );


  //--- Outputs---//

  // output 1: fval
  plhs[0] = mxCreateLogicalScalar(bTerminate);

}

void CommandDebug_SetMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Expected Input:
  //  cmd
  //  class handle
  //  match points                 (3 x Ns double)
  //  match normals                (3 x Ns double)
  //  

  // Get the class instance from the second input
  mexInterface_AlgDirICP_GIMLOP_Mesh& obj =
    *convertMat2Ptr<mexInterface_AlgDirICP_GIMLOP_Mesh>(prhs[1]);

  // Check parameters
  if (nlhs != 0 || nrhs != 4)
  {
    MEX_ERROR("SetMatches: requires 0 outputs, 4 inputs.");
  }

  // Parse sample data
  vctDynamicVector<vct3> matchPts;
  vctDynamicVector<vct3> matchNorms;

  MEX_DEBUG(" ...MatchPts\n");
  if (mxGetN(prhs[2]) > 0)
  {
    Parse3DVectorArray_Double(prhs[2], matchPts);
  }
  MEX_DEBUG(" ...MatchNorms\n");
  if (mxGetN(prhs[3]) > 0)
  {
    Parse3DVectorArray_Double(prhs[3], matchNorms);
  }

  MEX_DEBUG(" ...done\n");

  // override matches
  obj.pAlg->matchPts = matchPts;
  obj.pAlg->matchNorms = matchNorms;
}

//void CommandTest(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
//{
//    // Expected Input:
//    //  cmd
//    //  class handle
//
//    // Get the class instance pointer from the second input
//    algDirICP_StdICP_Mesh* pAlg =
//        convertMat2Ptr<algDirICP_StdICP_Mesh>(prhs[1]);
//
//    // Check parameters
//    if (nlhs != 0 || nrhs != 2)
//    {
//        MEX_ERROR("Test: requires 0 outputs, 2 inputs.");
//    }
//
//    std::string str = pAlg->test();
//    MEX_PRINT(str.c_str());
//}
//
//void CommandTest2(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
//{
//    // Expected Input:
//    //  cmd
//    //  class handle
//    //  F = [R,t]
//
//    // Check parameters
//    if (nlhs != 0 || nrhs != 3)
//    {
//        MEX_ERROR("CommandTest2: requires 0 outputs, 3 inputs.");
//    }
//
//    // Get the class instance pointer from the second input
//    algDirICP_StdICP_Mesh* pAlg =
//        convertMat2Ptr<algDirICP_StdICP_Mesh>(prhs[1]);
//
//    //--- Inputs ---//
//    vctFrm3 F;
//    ParseTransform(prhs[2], F);
//}

//void CommandTemplate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
//{
//  // Expected Input:
//  //  cmd
//  //  class handle
//
//  // Get the class instance pointer from the second input
//  algDirICP_StdICP_Mesh* pAlg = 
//    convertMat2Ptr<algDirICP_StdICP_Mesh>( prhs[1] );
//}