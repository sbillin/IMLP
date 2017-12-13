#ifndef _matlabExtras_h
#define _matlabExtras_h
#include "mex.h"

#define MEX_ERROR(x) { \
    /*mexErrMsgTxt((std::string("MEX_ERROR: ") + x + "\n").c_str()); \*/ \
    std::stringstream ss_mex_error; \
    ss_mex_error << "MEX_ERROR: " << x << "\n" \
        << " File: " << __FILE__ << " Line: " << __LINE__ << "\n"; \
    mexErrMsgTxt(ss_mex_error.str().c_str()); \
    mexEvalString("drawnow;"); \
    }

#define MEX_WARNING(x) { \
    mexWarnMsgTxt((std::string("MEX_WARNING: ") + x + "\n").c_str());   \
    mexEvalString("drawnow;"); \
    }

#define MEX_PRINT(x) { \
    mexPrintf(x);   \
    mexEvalString("drawnow;"); \
    }

#endif
