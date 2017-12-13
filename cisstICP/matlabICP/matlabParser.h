#ifndef _matlabParser_h
#define _matlabParser_h
#include "mex.h"
#include "cisstVector.h"
#include "matlabExtras.h"

// mxGetNumberOfElements( pArray )
// mxGetScalar
// mxGetPr
// mxCreateDoubleMatrix

void ParseTransform(const mxArray *pArray, vctFrm3 &F)
{
  // Expected input: [R,t]

  if (!mxIsDouble(pArray) ||
    mxIsComplex(pArray) ||
    mxGetM(pArray) != 3 ||
    mxGetN(pArray) != 4)
  {
    MEX_ERROR("ParseTransform: invalid input");
  }

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);

  // build transform
  vctRot3 R(
    pData[0], pData[3], pData[6],
    pData[1], pData[4], pData[7],
    pData[2], pData[5], pData[8]
    );
  vct3 t(pData[9], pData[10], pData[11]);
  F.Rotation() = R;
  F.Translation() = t;

  //mexPrintf("R:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
  //mexPrintf("t:\n %0.4f %0.4f %0.4f\n", pData[9], pData[10], pData[11]);
}

void ParseRotations(const mxArray *pArray, vctDynamicVector<vctRot3> &R)
{
  // Expected input: (3 x 3 x N double)

  mwSize nDimNum = mxGetNumberOfDimensions(pArray);
  const mwSize *pDims = mxGetDimensions(pArray);

  if (!mxIsDouble(pArray) || mxIsComplex(pArray) ||
    nDimNum != 3 || pDims[0] != 3 || pDims[1] != 3 || pDims[2] < 1)
  {
    MEX_ERROR("ParseRotations: invalid input");
  }
  unsigned int nRot = pDims[2];
  R.SetSize(nRot);

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);
  double *pRot;  // pointer to start of next rotation matrix

  // load rotation data
  for (unsigned int i = 0; i < nRot; i++)
  {
    pRot = pData + i * 9;
    R[i] = vctRot3(
      pRot[0], pRot[3], pRot[6],
      pRot[1], pRot[4], pRot[7],
      pRot[2], pRot[5], pRot[8]
      );
  }

  //char buf[300];
  //sprintf(buf, "R:\n %0.12f %0.12f %0.12f \n %0.12f %0.12f %0.12f \n %0.12f %0.12f %0.12f \n",
  //  pRot[0], pRot[3], pRot[6],
  //  pRot[1], pRot[4], pRot[7],
  //  pRot[2], pRot[5], pRot[8]);
  //MEX_PRINT(buf);

  //mexPrintf("R:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
}

void ParseDouble(const mxArray *pArray, double &X)
{
  unsigned int M = mxGetM(pArray);
  unsigned int N = mxGetN(pArray);

  if (!mxIsDouble(pArray) ||
    mxIsComplex(pArray))
  {
    MEX_ERROR("ParseDouble: invalid input type");
  }
  if (M != 1 || N != 1)
  {
    MEX_ERROR("ParseDouble: invalid input size");
  }

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);

  // transfer data
  X = *pData;
}

void ParseLogical(const mxArray *pArray, bool &X)
{
  if (!mxIsLogicalScalar(pArray))
  {
    MEX_ERROR("ParseLogical: invalid input type");
  }

  // create a pointer to logical data
  bool *pData = mxGetLogicals(pArray);

  // transfer data
  X = *pData;
}

void ParseInteger(const mxArray *pArray, int &X)
{
  unsigned int M = mxGetM(pArray);
  unsigned int N = mxGetN(pArray);

  if (!mxIsInt32(pArray) ||
    mxIsComplex(pArray))
  {
    MEX_ERROR("ParseInteger: invalid input type");
  }
  if (M != 1 || N != 1)
  {
    MEX_ERROR("ParseInteger: invalid input size");
  }

  // create a pointer to real matrix data
  //double *pData = mxGetPr(pArray);

  // create a pointer to integer data
  int *pData = static_cast<int*> (mxGetData(pArray));  

  // transfer data
  X = *pData;
}


void ParseVector_Double(const mxArray *pArray, vctDoubleVec &X)
{
  // Expected input: Vector ~ Nx1 or 1xN matrix

  unsigned int M = mxGetM(pArray);
  unsigned int N = mxGetN(pArray);
  unsigned int nSamps;

  if (!mxIsDouble(pArray) || mxIsComplex(pArray))
  {
    MEX_ERROR("ParseVector_Double: invalid input");
  }
  if (M >= N)
  {
    nSamps = M;
    if (N != 1 || M < 1)
    {
      MEX_ERROR("ParseVector_Double: invalid input");
    }
  }
  else
  {
    nSamps = N;
    if (M != 1 || N < 1)
    {
      MEX_ERROR("ParseVector_Double: invalid input");
    }
  }

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);

  // transfer data
  X.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    X.Element(i) = *pData;
    pData++;
  }
}


void ParseVector_Int(const mxArray *pArray, vctIntVec &X)
{
  // Expected input: Vector ~ Nx1 or 1xN matrix  (Int)

  unsigned int M = mxGetM(pArray);
  unsigned int N = mxGetN(pArray);
  unsigned int nSamps;

  if (!mxIsInt32(pArray))
  {
    MEX_ERROR("ParseVector_Int: invalid input");
  }
  if (M >= N)
  {
    nSamps = M;
    if (N != 1 || M < 1)
    {
      MEX_ERROR("ParseVector_Int: invalid input");
    }
  }
  else
  {
    nSamps = N;
    if (M != 1 || N < 1)
    {
      MEX_ERROR("ParseVector_Int: invalid input");
    }
  }

  // create a pointer to integer data
  int *pData = static_cast<int*> (mxGetData(pArray));

  // transfer data
  X.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    X.Element(i) = *pData;
    pData++;
  }
}


void Parse2DVectorArray_Double(const mxArray *pArray,
  vctDynamicVector<vct2> &X)
{
  // Expected input: Array ~ 2xN matrix

  unsigned int nSamps = mxGetN(pArray);

  if (!mxIsDouble(pArray) ||
    mxIsComplex(pArray) ||
    mxGetM(pArray) != 2 ||
    nSamps < 1)
  {
    std::stringstream ss;
    ss << "Parse2DVectorArray: invalid input:" << std::endl
      << "  !mxIsDouble(pArray): " << !mxIsDouble(pArray) << std::endl
      << "  mxIsComplex(pArray): " << mxIsComplex(pArray) << std::endl
      << "  mxGetM(pArray): " << mxGetM(pArray) << std::endl
      << "  nSamps: " << nSamps << std::endl;
    MEX_ERROR(ss.str());
    //MEX_ERROR("Parse2DVectorArray: invalid input");
  }

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);

  // transfer data
  X.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    X[i].Element(0) = *pData;
    pData++;
    X[i].Element(1) = *pData;
    pData++;
  }
}


void Parse3DVectorArray_Double(const mxArray *pArray,
  vctDynamicVector<vct3> &X)
{
  // Expected input: Array ~ 3xN matrix

  unsigned int nSamps = mxGetN(pArray);

  if (!mxIsDouble(pArray) ||
    mxIsComplex(pArray) ||
    mxGetM(pArray) != 3 ||
    nSamps < 1)
  {
    std::stringstream ss;
    ss << "mxIsDouble: " << mxIsDouble(pArray) << std::endl
      << " mxIsComplex: " << mxIsComplex(pArray) << std::endl
      << " mxGetM: " << mxGetM(pArray) << std::endl
      << " nSamps: " << nSamps << std::endl;
    MEX_PRINT(ss.str().c_str());
    MEX_ERROR("Parse3DVectorArray: invalid input");
  }

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);

  // transfer data
  X.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    X[i].Element(0) = *pData;
    pData++;
    X[i].Element(1) = *pData;
    pData++;
    X[i].Element(2) = *pData;
    pData++;
  }
}

void Parse3DVectorArray_Int(const mxArray *pArray,
  vctDynamicVector<vctInt3> &X)
{
  // Expected input: Array ~ 3xN matrix

  unsigned int nSamps = mxGetN(pArray);

  if (!mxIsInt32(pArray) ||
    mxIsComplex(pArray) ||
    mxGetM(pArray) != 3 ||
    nSamps < 1)
  {
    MEX_ERROR("Parse3DVectorArray: invalid input");
  }

  // create a pointer to real matrix data
  int *pData = static_cast<int*> (mxGetData(pArray));

  // transfer data
  X.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    X[i].Element(0) = *pData;
    pData++;
    X[i].Element(1) = *pData;
    pData++;
    X[i].Element(2) = *pData;
    pData++;
  }
}

void ParseMatrixArray3x3_Double(const mxArray *pArray, vctDynamicVector<vct3x3> &X)
{
  // Expected input: (3 x 3 x N double)

  mwSize nDimNum = mxGetNumberOfDimensions(pArray);
  const mwSize *pDims = mxGetDimensions(pArray);

  if (!mxIsDouble(pArray) || mxIsComplex(pArray) ||
    nDimNum != 3 || pDims[0] != 3 || pDims[1] != 3 || pDims[2] < 1)
  {
    MEX_ERROR("ParseMatrixArray3x3_Double: invalid input");
  }
  unsigned int nMat = pDims[2];
  X.SetSize(nMat);

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);
  double *pMat;  // pointer to start of next matrix

  // load data
  for (unsigned int i = 0; i < nMat; i++)
  {
    pMat = pData + i * 9;
    X[i] = vct3x3(
      pMat[0], pMat[3], pMat[6],
      pMat[1], pMat[4], pMat[7],
      pMat[2], pMat[5], pMat[8]
      );
  }

  //mexPrintf("M:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
}


void ParseMatrixArray2x3_Double(const mxArray *pArray, vctDynamicVector<vct2x3> &X)
{
  // Expected input: (2 x 3 x N double)

  mwSize nDimNum = mxGetNumberOfDimensions(pArray);
  const mwSize *pDims = mxGetDimensions(pArray);

  if (!mxIsDouble(pArray) || mxIsComplex(pArray) ||
    nDimNum != 3 || pDims[0] != 2 || pDims[1] != 3 || pDims[2] < 1)
  {
    MEX_ERROR("ParseMatrixArray2x3_Double: invalid input");
  }
  unsigned int nMat = pDims[2];
  X.SetSize(nMat);

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);
  double *pMat;  // pointer to start of next matrix

  // load data
  for (unsigned int i = 0; i < nMat; i++)
  {
    pMat = pData + i * 6;
    X[i] = vct2x3(
      pMat[0], pMat[2], pMat[4],
      pMat[1], pMat[3], pMat[5]
      );
  }

  //mexPrintf("M:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
}


void ParseMatrixArray3x2_Double(const mxArray *pArray, vctDynamicVector<vct3x2> &X)
{
  // Expected input: (3 x 2 x N double)

  mwSize nDimNum = mxGetNumberOfDimensions(pArray);
  const mwSize *pDims = mxGetDimensions(pArray);

  if (!mxIsDouble(pArray) || mxIsComplex(pArray) ||
    nDimNum != 3 || pDims[0] != 3 || pDims[1] != 2 || pDims[2] < 1)
  {
    MEX_ERROR("ParseMatrixArray3x2_Double: invalid input");
  }
  unsigned int nMat = pDims[2];
  X.SetSize(nMat);

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);
  double *pMat;  // pointer to start of next matrix

  // load data
  for (unsigned int i = 0; i < nMat; i++)
  {
    pMat = pData + i * 6;
    X[i] = vct3x2(
      pMat[0], pMat[3],
      pMat[1], pMat[4],
      pMat[2], pMat[5]
      );
  }

  //mexPrintf("M:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
}

void ParseMatrixArray_KentL(const mxArray *pArray, vctDynamicVector<vct3x2> &L)
{
  // Expected input: (2 x 3 x N double)
  //  need to convert this to (3 x 2 x N)

  mwSize nDimNum = mxGetNumberOfDimensions(pArray);
  const mwSize *pDims = mxGetDimensions(pArray);

  if (!mxIsDouble(pArray) || mxIsComplex(pArray) ||
    nDimNum != 3 || pDims[0] != 2 || pDims[1] != 3 || pDims[2] < 1)
    //nDimNum < 2 || pDims[0] != 2 || pDims[1] != 3 || pDims[2] < 1)    
  {
    std::stringstream ss;
    ss << "ParseMatrixArray_KentL: invalid input" << std::endl
      << " mxIsDouble: " << mxIsDouble(pArray) << std::endl
      << " mxIsComplex: " << mxIsComplex(pArray) << std::endl
      << " nDimNum: " << nDimNum << std::endl
      << " pDims: " << pDims[0] << " " << pDims[1] << " " << pDims[2] << std::endl;
    MEX_ERROR(ss.str().c_str());
  }
  unsigned int nMat = pDims[2];
  L.SetSize(nMat);

  // create a pointer to real matrix data
  double *pData = mxGetPr(pArray);
  double *pMat;  // pointer to start of next matrix

  // load data
  for (unsigned int i = 0; i < nMat; i++)
  {
    pMat = pData + i * 6;
    // store 2x3 input as 3x2 in L
    L[i] = vct3x2(
      pMat[0], pMat[1],
      pMat[2], pMat[3],
      pMat[4], pMat[5]
      );
  }

  //mexPrintf("M:\n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n %0.4f %0.4f %0.4f \n",
  //  pData[0], pData[3], pData[6],
  //  pData[1], pData[4], pData[7],
  //  pData[2], pData[5], pData[8]);
}


#endif // _matlabParser_h
