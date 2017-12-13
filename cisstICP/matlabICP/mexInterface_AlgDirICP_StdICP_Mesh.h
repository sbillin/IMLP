#ifndef _mexInterface_AlgDirICP_StdICP_Mesh_h
#define _mexInterface_AlgDirICP_StdICP_Mesh_h

#include "matlabClassHandle.h"
#include "mex.h"
#include "matlabParser.h"
#include "matlabExtras.h"

#include "algDirICP_StdICP_Mesh.h"


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
void CommandSetSamples(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_InitializeParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_ComputeMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_RegisterMatches(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_UpdateParameters_PostRegister(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_EvaluateErrorFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CommandICP_Terminate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void CommandTest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void CommandTest2(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

class mexInterface_AlgDirICP_StdICP_Mesh
{

  //--- Parameters ---//

public:

  algDirICP_StdICP_Mesh *pAlg;


  //--- Methods ---//

public:

  // constructor
	mexInterface_AlgDirICP_StdICP_Mesh() 
    : pAlg(NULL)
  {}

	~mexInterface_AlgDirICP_StdICP_Mesh()
  {
    if (pAlg)
    {
      delete pAlg;
      pAlg = NULL;
    }
  }

};


#endif