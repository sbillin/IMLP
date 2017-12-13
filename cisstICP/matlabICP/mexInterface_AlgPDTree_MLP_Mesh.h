#ifndef _mexInterface_AlgPDTree_MLP_Mesh_h
#define _mexInterface_AlgPDTree_MLP_Mesh_h

#include "mex.h"

#include "algPDTree_MLP_Mesh.h"


#ifdef DEBUG_MEX
  #define MEX_DEBUG(x) MEX_PRINT((std::string("MEX_DEBUG: ") + x + "\n").c_str())
#else
  #define MEX_DEBUG(x)
#endif


// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// Command Interface
void CommandNew(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandDelete(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandInitialize(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandComputeMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

class mexInterface_AlgPDTree_MLP_Mesh
{

  //--- Parameters ---//

public:

  PDTree_Mesh        *pTree;
  algPDTree_MLP_Mesh *pAlg;

  // covariance tree inputs
  vctDynamicVector<vct3>      V;    // mesh vertices
  vctDynamicVector<vctInt3>   T;    // mesh triangles
  vctDynamicVector<vct3>      Tn;   // mesh triangle normals

  vctDynamicVector<vct3x3>  TriangleCov;        // mesh triangle covariances
  vctDynamicVector<vct3>    TriangleCovEig;     // mesh triangle covariance eigenvalues (in decreasing order)


  // per sample-set I/O
  unsigned int nSamples;

  vctDynamicVector<vct3>    samplePtsXfmd;
  vctDynamicVector<vct3x3>  sampleCovXfmd;

  vctDynamicVector<vct3>  matchPts;
  vctDynamicVector<vct3>  matchNorms;
  vctDynamicVector<int>   matchDatums;
  vctDoubleVec            matchErrors;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
  mexInterface_AlgPDTree_MLP_Mesh()
    : pAlg(NULL)
  {}

  ~mexInterface_AlgPDTree_MLP_Mesh()
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
