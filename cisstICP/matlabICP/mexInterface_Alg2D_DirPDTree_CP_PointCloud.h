#ifndef _mexInterface_Alg2D_DirPDTree_CP_PointCloud_h
#define _mexInterface_Alg2D_DirPDTree_CP_PointCloud_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "alg2D_DirPDTree_CP_PointCloud.h"


// debug
//#define DEBUG_MEX

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

class mexInterface_Alg2D_DirPDTree_CP_PointCloud
{

  //--- Parameters ---//

public:

  DirPDTree2D_PointCloud *pDirTree;
  alg2D_DirPDTree_CP_PointCloud *pAlg;

  // model shape data
  // covariance tree inputs
  vctDynamicVector<vct2> points;
  vctDynamicVector<vct2> pointOrientations;


  // per sample-set I/O
  unsigned int nSamples;

  vctDynamicVector<vct2> samplePtsXfmd;
  vctDynamicVector<vct2> sampleNormsXfmd;

  vctDynamicVector<vct2> matchPts;
  vctDynamicVector<vct2> matchNorms;
  vctDynamicVector<int>  matchDatums;
  vctDoubleVec  matchErrors;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
	mexInterface_Alg2D_DirPDTree_CP_PointCloud() 
    : pAlg(NULL)
  {}

	~mexInterface_Alg2D_DirPDTree_CP_PointCloud()
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