#ifndef _mexInterface_Alg2D_DirPDTree_vonMises_Edges_h
#define _mexInterface_Alg2D_DirPDTree_vonMises_Edges_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "alg2D_DirPDTree_vonMises_Edges.h"


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

class mexInterface_Alg2D_DirPDTree_vonMises_Edges
{

  //--- Parameters ---//

public:

  DirPDTree2D_Edges *pDirTree;
  alg2D_DirPDTree_vonMises_Edges *pAlg;

  // model shape data
  // covariance tree inputs
  vctDynamicVector<vct2>  edgesV1;
  vctDynamicVector<vct2>  edgesV2;
  vctDynamicVector<vct2>  edgesNorm;

  // sample data
  vctDynamicVector<vct2> samplePtsXfmd;
  vctDynamicVector<vct2> sampleNormsXfmd;

  // noise models of sample data
  vctDynamicVector<double> sigma2;  // position variance
  vctDynamicVector<double> k;       // orientation concentration  

  // match filtering
  //double match_ThetaMax;


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
	mexInterface_Alg2D_DirPDTree_vonMises_Edges() 
    : pAlg(NULL)//, match_ThetaMax(cmnPI)
  {}

	~mexInterface_Alg2D_DirPDTree_vonMises_Edges()
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