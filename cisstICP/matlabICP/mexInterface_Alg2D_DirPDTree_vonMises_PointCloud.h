#ifndef _mexInterface_Alg2D_DirPDTree_vonMises_PointCloud_h
#define _mexInterface_Alg2D_DirPDTree_vonMises_PointCloud_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "alg2D_DirPDTree_vonMises_PointCloud.h"


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
void CommandSetNoiseModel(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);

class mexInterface_Alg2D_DirPDTree_vonMises_PointCloud
{

  //--- Parameters ---//

public:

  DirPDTree2D_PointCloud *pDirTree;
  alg2D_DirPDTree_vonMises_PointCloud *pAlg;

  // model shape data
  // covariance tree inputs
  vctDynamicVector<vct2> points;
  vctDynamicVector<vct2> pointOrientations;

  // sample data
  vctDynamicVector<vct2> samplePtsXfmd;
  vctDynamicVector<vct2> sampleNormsXfmd;
  // sample data noise models
  vctDynamicVector<double> sigma2;  // position variance
  vctDynamicVector<double> k;       // orientation concentration  


  // per sample-set I/O
  unsigned int nSamples;

  // output: matches
  vctDynamicVector<vct2> matchPts;
  vctDynamicVector<vct2> matchNorms;
  vctDynamicVector<int>  matchDatums;
  vctDoubleVec           matchErrors;
  vctDynamicVector<bool> matchPermitted;

  double minNodesSearched, maxNodesSearched, avgNodesSearched;



  //--- Methods ---//

public:

  // constructor
	mexInterface_Alg2D_DirPDTree_vonMises_PointCloud() 
    : pAlg(NULL)
  {}

	~mexInterface_Alg2D_DirPDTree_vonMises_PointCloud()
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