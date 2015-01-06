// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************
#include "utilities.h"

#include <algorithm>
#include "cisstNumerical.h"
#include "Wm5NoniterativeEigen3x3.h"  // WildMagic5

//#define ENABLE_UTILITIES_DEBUG

void ComputeCovEigenDecomposition_NonIter(const vct3x3 &M, vct3 &eigenValues, vct3x3 &eigenVectors)
{
  // Calls the non-iterative eigen solver of the WildMagic5 library
  // 
  // eigen values in descending order
  // eigen vectors listed by column   (has determinant = 1)
  //

  bool rowMajor = 1;
  Wm5::Matrix3<double> A(M.Pointer(0,0), rowMajor);

  Wm5::NoniterativeEigen3x3<double> solver(A);

  //  change order from ascending to descending
  eigenValues.Assign(
    solver.GetEigenvalue(2),
    solver.GetEigenvalue(1),
    solver.GetEigenvalue(0));
  eigenVectors.Column(0).Assign(solver.GetEigenvector(2));
  eigenVectors.Column(1).Assign(solver.GetEigenvector(1));
  eigenVectors.Column(2).Assign(solver.GetEigenvector(0));

  // check that eigen vector matrix is a proper rotation matrix
  //  NOTE: this is probably not important, but a nice property to have
  double det = 
    solver.GetEigenvector(2).Dot(solver.GetEigenvector(1).Cross(
    solver.GetEigenvector(0)));

  if (det < 0.0)
  {
    eigenVectors.Column(2).Multiply(-1.0);
  }

  //for (i = 0; i < 3; ++i)
  //{
  //  mNoniterativeEigenvalues[i] = noniterativeSolver.GetEigenvalue(i);
  //  mNoniterativeEigenvectors[i] = noniterativeSolver.GetEigenvector(i);
  //  result = A*mNoniterativeEigenvectors[i] -
  //    mNoniterativeEigenvalues[i] * mNoniterativeEigenvectors[i];
  //  length = result.Length();
  //  if (length > noniterativeError)
  //  {
  //    noniterativeError = length;
  //  }
  //}
  //noniterativeDeterminant = Mathf::FAbs(
  //  mNoniterativeEigenvectors[0].Dot(mNoniterativeEigenvectors[1].Cross(
  //  mNoniterativeEigenvectors[2])));
}


void ComputeCovEigenDecomposition_SVD(const vct3x3 &M, vct3 &eigenValues, vct3x3 &eigenVectors)
{
  //
  // eigen values in descending order
  // eigen vectors listed by column
  //

  // Compute SVD of M
  //   M = U*diag(S)*V'   where U = V
  //
  //  NOTE: matrices must be column major
  //
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Mcopy;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  static nmrSVDFixedSizeData<3, 3, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try
  {
    Mcopy.Assign(M);  // must use "assign" rather than equals to properly transfer between different vector orderings
    nmrSVD(Mcopy, U, eigenValues, Vt, workspace);
}
  catch (...)
  {
    std::cout << std::endl << "========> ERROR: ComputeCovEigenDecomposition_SVD() failed!" << std::endl << std::endl;
    assert(0);
  }

  // copy eigen vectors to output using "assign" to preserve ordering
  //  of output matrix
  eigenVectors.Assign(U);

#ifdef ENABLE_UTILITIES_DEBUG
  if (eigenValues(2) < 0.0 || eigenValues(2) > eigenValues(1) || eigenValues(1) > eigenValues(0))
  {
    std::cout << std::endl << "========> ERROR: ComputeCovEigenDecomposition_SVD() eigen values misordered or less than zero!" << std::endl << std::endl;
    assert(0);
  }
#endif
}

void ComputeCovEigenDecomposition_SEP(const vct3x3 &M, vct3 &eigenValues, vct3x3 &eigenVectors)
{
  //
  // Symmetric Eigenproblems (SEP) 
  //
  // Computes all the eigen values and eigen vectors of a symmetric matrix
  //   A = V D V^T. The eigen values are sorted in ascending order. This
  //   function uses LAPACK dsyevr.
  //
  //   eigen values in descending order
  //   eigen vectors listed by column
  //
  // This method is much less efficient than the SVD method
  //   SVD: 6.42449 (sec)
  //   SEP: 9.30579 (sec)
  //

  static vctDynamicMatrix<double> Mcopy(3, 3, VCT_COL_MAJOR);
  static vctDynamicMatrix<double> eigVct(3, 3, VCT_COL_MAJOR);
  static vctDynamicVector<double> eigVal(3);
  static nmrSymmetricEigenProblem::Data workspace = nmrSymmetricEigenProblem::Data(Mcopy, eigVal, eigVct);

  Mcopy.Assign(M);
  if (nmrSymmetricEigenProblem::EFAILURE == nmrSymmetricEigenProblem(Mcopy, eigVal, eigVct, workspace))
  {
    std::cout << std::endl << "========> ERROR: ComputeCovEigenDecomposition_SEP() failed!" << std::endl << std::endl;
    assert(0);
  }

  // assign to fixed size containers
  //  change order from ascending to descending 
  eigenValues.Assign(eigVal[2],eigVal[1],eigVal[0]);
  eigenVectors.Column(0).Assign(eigVct.Column(2));
  eigenVectors.Column(1).Assign(eigVct.Column(1));
  eigenVectors.Column(2).Assign(eigVct.Column(0));
}

void ComputeCovInverse_NonIter(const vct3x3 &M, vct3x3 &Minv)
{
  vct3    eigenValues;
  vct3x3  eigenVectors;
  
  // Compute Minv
  //   M = U*diag(S)*V'   where U = V
  //   Minv = V*diag(1/S)*V'
  
  ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> V_Sinv;
  static vct3 Sinv;
  Sinv[0] = 1.0 / eigenValues[0];
  Sinv[1] = 1.0 / eigenValues[1];
  Sinv[2] = 1.0 / eigenValues[2];
  V_Sinv.Column(0) = eigenVectors.Column(0)*Sinv[0];
  V_Sinv.Column(1) = eigenVectors.Column(1)*Sinv[1];
  V_Sinv.Column(2) = eigenVectors.Column(2)*Sinv[2];
  Minv.Assign(V_Sinv * eigenVectors.TransposeRef());
}

void ComputeCovInverse_SVD(const vct3x3 &M, vct3x3 &Minv)
{
  // Compute SVD of M
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Mcopy;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  static vct3 S;
  static nmrSVDFixedSizeData<3, 3, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try
  {
    Mcopy.Assign(M);
    nmrSVD(Mcopy, U, S, Vt, workspace);
  }
  catch (...)
  {
    assert(0);
  }

  // Compute Minv
  //   M = U*diag(S)*V'   where U = V
  //   Minv = V*diag(1/S)*U' = U*diag(1/S)*V'
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Sinv_Ut;
  static vct3 Sinv;
  Sinv[0] = 1.0 / S[0];
  Sinv[1] = 1.0 / S[1];
  Sinv[2] = 1.0 / S[2];
  Sinv_Ut.Row(0) = Sinv[0] * Vt.Row(0);
  Sinv_Ut.Row(1) = Sinv[1] * Vt.Row(1);
  Sinv_Ut.Row(2) = Sinv[2] * Vt.Row(2);
  Minv.Assign(U*Sinv_Ut);
}

// this is only slightly slower than the non-iterative
//  method based on eigen decomposition, but much faster
//  than the SVD method
void ComputeCovInverse_Nmr(const vct3x3 &M, vct3x3 &Minv)
{
  Minv = M;
  nmrInverse(Minv); // computes inverse in-place
}

void ComputeCovEigenValues_SVD(const vct3x3 &M, vct3 &eigenValues)
{
  // Compute SVD of M
  //   M = U*diag(S)*V'   where U = V
  //
  //  NOTE: matrices must be column major
  //        eigen values are in descending order
  //
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Mcopy;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  static nmrSVDFixedSizeData<3, 3, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try
  {
    Mcopy.Assign(M);
    nmrSVD(Mcopy, U, eigenValues, Vt, workspace);
  }
  catch (...)
  {
    std::cout << std::endl << "========> ERROR: ComputeCovEigenValues_SVD() failed!" << std::endl << std::endl;
    assert(0);
  }

}

void ComputeCovEigenValues_Trig(const vct3x3 &M, vct3 &eigenValues)
{
  //
  // Code derived from an algorithm posted on Wikipedia for computing 
  //  eigenvalues of a 3x3 real symmetric matrix based on the article:
  //
  //  Oliver Smith, "Eigenvalues of a symmetric 3x3 matrix", 
  //   J. Comm. of ACS Vol. 4, Issue 4, pg. 168, 1961
  //
  //  NOTE: eigen values are in descending order
  //

  static vctDeterminant<3> detCalc;

  double p, p1, p2;
  double q, q1, q2, q3;
  double r, phi;
  double eig1, eig2, eig3;
  vct3x3 B;

  p1 = M.Element(0, 1) * M.Element(0, 1)
    + M.Element(0, 2) * M.Element(0, 2)
    + M.Element(1, 2) * M.Element(1, 2);

  if (p1 <= std::numeric_limits<double>::epsilon() * 10.0)
    //if (p1 == 0.0)    
  {
    // M is diagonal
    eig1 = M.Element(0, 0);
    eig2 = M.Element(1, 1);
    eig3 = M.Element(2, 2);

    // sort eigenvalues in descending order
    if (eig1 > eig2)
    {
      if (eig3 > eig1) std::swap(eig1, eig3);
    }
    else
    {
      if (eig2 > eig3) std::swap(eig1, eig2);
      else std::swap(eig1, eig3);
    }
    // now eig1 is largest; order the remaining two
    if (eig3 > eig2) std::swap(eig2, eig3);
  }
  else
  {
    q = M.Trace() / 3.0;
    q1 = M.Element(0, 0) - q;
    q2 = M.Element(1, 1) - q;
    q3 = M.Element(2, 2) - q;
    p2 = q1*q1 + q2*q2 + q3*q3 + 2.0 * p1;
    p = sqrt(p2 / 6.0);
    B = (1.0 / p) * (M - q * vct3x3::Eye());
    r = detCalc.Compute(B) / 2.0;

    // In exact arithmetic for a symmetric matrix - 1 <= r <= 1
    //  but computation error can leave it slightly outside this range.
    if (r <= -1.0)
    {
      phi = cmnPI / 3.0;
    }
    else if (r >= 1)
    {
      phi = 0.0;
    }
    else
    {
      phi = acos(r) / 3.0;
    }

    // the eigenvalues satisfy eig1 >= eig2 >= eig3
    eig1 = q + 2.0 * p * cos(phi);
    eig3 = q + 2.0 * p * cos(phi + (2.0 * cmnPI / 3.0));
    eig2 = 3.0 * q - eig1 - eig3;   // since trace(A) = eig1 + eig2 + eig3
  }

  eigenValues.Assign(eig1, eig2, eig3);
}