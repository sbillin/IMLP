#ifndef _mexInterface_AlgPDTree_CP_Mesh_h
#define _mexInterface_AlgPDTree_CP_Mesh_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "algPDTree_CP_Mesh.h"


#ifdef DEBUG_MEX
#define MEX_DEBUG(x) MEX_PRINT((std::string("MEX_DEBUG: ") + x + "\n").c_str())
#else
#define MEX_DEBUG(x)
#endif;


// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// Command Interface
void CommandNew(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandDelete(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandInitialize(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandComputeMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

class mexInterface_AlgPDTree_CP_Mesh
{

  //--- Parameters ---//

public:

  PDTree_Mesh             *pTree;
  algPDTree_CP_Mesh *pAlg;

  // covariance tree inputs
  vctDynamicVector<vct3>      V;    // mesh vertices
  vctDynamicVector<vctInt3>   T;    // mesh triangles
  vctDynamicVector<vct3>      Tn;   // mesh triangle normals


  // per sample-set I/O
  unsigned int nSamples;

  vctDynamicVector<vct3> samplePtsXfmd;

  vctDynamicVector<vct3>  matchPts;
  vctDynamicVector<vct3>  matchNorms;
  vctDynamicVector<int>   matchDatums;
  vctDoubleVec            matchErrors;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
  mexInterface_AlgPDTree_CP_Mesh()
    : pAlg(NULL)
  {}

  ~mexInterface_AlgPDTree_CP_Mesh()
  {
    if (pAlg)
    {
      delete pAlg;
      pAlg = NULL;
    }
  }

  void ComputeMatches();
};
#endif