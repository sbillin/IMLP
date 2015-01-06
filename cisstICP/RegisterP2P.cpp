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


#include "RegisterP2P.h"
#include "utilities.h"


vct3 vctCentroid(const vctDynamicVector<vct3>& A)
{
  vct3 mean(0.0);
  unsigned int N = A.size();
  for (unsigned int i = 0; i < N; i++)
  {
    mean.Add(A[i]);
  }
  mean.Divide(N);
  return mean;
};

// Compute weighted vector mean:  Amean = Sum(Wi*Ai)/Sum(Wi)
vct3 vctWeightedMean(const vctDynamicVector<vct3>& A, const vctDoubleVec &W)
{
  vct3 mean(0.0);
  unsigned int N = A.size();
  assert(N == W.size());
  for (unsigned int i = 0; i < N; i++)
  {
    mean.AddProductOf(W[i], A[i]);
  }
  mean.Divide(W.SumOfElements());
  return mean;
};

// Compute the least squares rigid body transform to minimize P2P distances:
//  minimize:  Sum( ||Yi - T*Xi||^2 )
void RegisterP2P_LSQ(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctFrm3 &F)
{
  unsigned int numPts = X.size();
  assert(numPts == Y.size());

  // weighted centroids of each point set
  vct3 Xmean = vctCentroid(X);
  vct3 Ymean = vctCentroid(Y);

  // recenter point sets about the centroids
  vctDynamicVector<vct3> Xp(numPts);
  vctDynamicVector<vct3> Yp(numPts);
  for (unsigned int i = 0; i < numPts; i++)
  {
    Xp[i] = X[i] - Xmean;
    Yp[i] = Y[i] - Ymean;
  }

  // Rotation
  RotateP2P_LSQ_SVD(Xp, Yp, F.Rotation());
  //RotateP2P_LSQ_Quaternion(Xp, Yp, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Non-linear total least squares Point-to-Point rigid body registration
//  based on Mahalanobis distance assuming a multivariate Gaussian
//  distribution on both the match position errors and sample position errors
//  and assuming independence between the match and sample errors.
//  Solved via a linearized iterative scheme based on the Jacobian and
//  using Lagrange multipliers to enforce the model constraint.
//
//  Minimizes:  sum( dXi'*inv(Mxi)*dXi + dYi'*inv(Myi)*dYi )
//                where Yi = R*Xi + t  (model constraint)
//                      dX = Xcalc - Xobs
//                      dY = Ycalc - Yobs
//
//   x   - observed sample values
//   y   - observed model values
//   Mxi - covariance matrix for dX  (assumed same for all dXi)
//   Myi - covariance matrix for dY  (assumed same for all dYi)
//
void RegisterP2P_TLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vct3x3 &Mxi, const vct3x3 &Myi,
  vctFrm3 &Fact)
{

  // termination conditions
  double dTheta_term = 0.001*cmnPI / 180.0;   // delta rotation (radians)
  double dt_term = 0.001;       // delta translation (mm)
  int    maxIter = 10;        // max iterations

  // Workspace for SVD computations
  vctMat A(6, 6, VCT_COL_MAJOR);
  vctMat U(6, 6, VCT_COL_MAJOR);
  vctVec S(6);
  vctMat Vt(6, 6, VCT_COL_MAJOR);
  vctDynamicVector<double> workspaceSVD6x6(nmrSVDDynamicData::WorkspaceSize(A));

  //// Workspace for Inverse computations
  //nmrInverseFixedSizeData<6, VCT_COL_MAJOR> workspaceInv6;

  unsigned int nSamps = x.size();

  vctMat J(3 * nSamps, 6, 0.0);
  vctMat Jt_Minv(6, 3 * nSamps);
  vctMat Jt_Minv_J(6, 6);
  vctVec neg_Jt_Minv_F0(6);
  vctVec F0(3 * nSamps);
  //vctMat Fxi( 3,3 );
  vct3x3 Fxi;
  vct3x3 Mi;
  vct3x3 Minv;
  vctRot3 dR;
  vctRodRot3 dAlpha;
  vct3 dt;
  unsigned int idx;

  // set last three columns of J to negative block identity
  for (unsigned int i = 0; i<3 * nSamps; i += 3)
  {
    J.Element(i, 3) = -1.0;
    J.Element(i + 1, 4) = -1.0;
    J.Element(i + 2, 5) = -1.0;
  }

  // Step 1: Initialize transform params
  vct3x3 R = vct3x3::Eye();
  vct3 t(0.0);
  vctDynamicMatrixRef<double> Rref(R);

  int numIter = 0;
  double dTheta = 1.0;
  double dt_norm = 1.0;

  while (((dTheta > dTheta_term) || (dt_norm > dt_term)) && (numIter++ < maxIter))
  {
    // Step 2: F0
    vct3 res;
    vctDynamicVector<vct3> x_rot(nSamps);
    for (unsigned int i = 0; i < nSamps; i++)
    {
      x_rot.Element(i) = R*x.Element(i);
      res = y.Element(i) - (x_rot.Element(i) + t);
      // stack residuals into a single vector [f1x f1y f1z f2x...]'  
      idx = 3 * i;
      F0.Element(idx++) = res.X();
      F0.Element(idx++) = res.Y();
      F0.Element(idx) = res.Z();
    }

    // Step 3: J & Fx
    //   Note: Fx is block diagonal with all sub-blocks being equal
    //         to each other (equal to Fxi)
    vctDynamicMatrixRef<double> Jref3x3;
    for (unsigned int i = 0; i < nSamps; i++)
    {
      // populate first 3 cols of J
      Jref3x3.SetRef(J, 3 * i, 0, 3, 3);
      skew(x_rot.Element(i), Jref3x3);  // Jref3x3 = skew(x_rot.Element(i))
    }
    Fxi = -R;

    // Step 4: solve dP by least squares
    //  Note: Mi & Minv are block diagonal with same block
    //        repeated; therefore just calculate this sub-block
    Mi = Fxi*Mxi*Fxi.TransposeRef() + Myi;
    Minv = Mi;
    nmrInverse(Minv); // computes inverse in-place 
    //nmrInverse(Minv, workspaceInv6);   // Anton (why doesn't this work?)

    vctDynamicMatrixRef<double> JrowsRef;
    vctDynamicMatrixRef<double> Jt_MinvRef;
    vctDynamicMatrixRef<double> MinvRef(Minv);
    for (unsigned int i = 0; i < nSamps; i++)
    { // Minv is block diagonal => only multiply the sub-blocks
      //  Jt_Minv:  6x3N
      //  J:        3Nx6
      //  M:        3Nx3N  (block diagonal)
      idx = 3 * i;
      JrowsRef.SetRef(J, idx, 0, 3, 6);         // next 3 rows of J
      Jt_MinvRef.SetRef(Jt_Minv, 0, idx, 6, 3);  // next 3 cols of Jt_Minv
      Jt_MinvRef = JrowsRef.TransposeRef()*MinvRef;
    }
    Jt_Minv_J = Jt_Minv * J;
    neg_Jt_Minv_F0 = Jt_Minv * (-F0);

    // solve Ax = b
    //  Note: most of the run-time is for this SVD call
    try
    {
      A.Assign(Jt_Minv_J);
      // this changes the storage order of A (must use assign instead)
      //A = Jt_Minv_J;
      nmrSVD(A, U, S, Vt, workspaceSVD6x6);
    }
    catch (...)
    {
      std::cerr << "ERROR: Compute SVD failed" << std::endl;
      assert(0);
    }

    // dP = V * Sinv * U' * neg_Jt_Minv_F0
    vctVec Ut_F0(6);
    vctVec Sinv_Ut_F0(6);
    vctVec dP(6);
    Ut_F0 = U.TransposeRef()*neg_Jt_Minv_F0;
    for (unsigned int i = 0; i < 6; i++)
    {
      Sinv_Ut_F0.Element(i) = Ut_F0.Element(i) / S.Element(i);
    }
    dP = Vt.TransposeRef()*Sinv_Ut_F0;
    dAlpha.X() = dP.Element(0);
    dAlpha.Y() = dP.Element(1);
    dAlpha.Z() = dP.Element(2);
    dt.X() = dP.Element(3);
    dt.Y() = dP.Element(4);
    dt.Z() = dP.Element(5);

    // Step 5: [R,t]
    dR = vctRot3(dAlpha);
    R = dR*R;
    t = t + dt;

    //std::cout << "dt: " << dt[0] << " " << dt[1] << " " << dt[2] << "  dA: " << dAlpha[0] << " " << dAlpha[1]
    //<< " " << dAlpha[2] << "  " << dAlpha.Norm()*180/cmnPI << std::endl;

    dt_norm = dt.Norm();
    dTheta = dAlpha.Norm();
  }

  Fact.Rotation() = R;
  Fact.Translation() = t;
}


// Non-linear total least squares Point-to-Point rigid body registration
//  based on Mahalanobis distance assuming a multivariate Gaussian
//  distribution on both the match position errors and sample position errors
//  and assuming independence between the match and sample errors.
//  Solved via a linearized iterative scheme based on the Jacobian and
//  using Lagrange multipliers to enforce the model constraint.
//
//  Minimizes:  sum( dXi'*inv(Mxi)*dXi + dYi'*inv(Myi)*dYi )
//                where Yi = R*Xi + t  (model constraint)
//                      dX = Xcalc - Xobs
//                      dY = Ycalc - Yobs
//
//   x   - observed sample values
//   y   - observed model values
//   Mxi - covariance matrix for dX  (assumed different for all dXi)
//   Myi - covariance matrix for dY  (assumed different for all dYi)
//
void RegisterP2P_TLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3x3> &Mxi, const vctDynamicVector<vct3x3> &Myi,
  vctFrm3 &Fact)
{

  // termination conditions
  double dTheta_term = 0.001*cmnPI / 180.0;   // delta rotation (radians)
  double dt_term = 0.001;       // delta translation (mm)
  int    maxIter = 10;        // max iterations

  // Workspace for SVD computations
  vctMat A(6, 6, VCT_COL_MAJOR);
  vctMat U(6, 6, VCT_COL_MAJOR);
  vctVec S(6);
  vctMat Vt(6, 6, VCT_COL_MAJOR);
  vctDynamicVector<double> workspaceSVD6x6(nmrSVDDynamicData::WorkspaceSize(A));

  //// Workspace for Inverse computations
  //nmrInverseFixedSizeData<6, VCT_COL_MAJOR> workspaceInv6;

  unsigned int nSamps = x.size();

  vctMat J(3 * nSamps, 6, 0.0);
  vctMat Jt_Minv(6, 3 * nSamps);
  vctMat Jt_Minv_J(6, 6);
  vctVec neg_Jt_Minv_F0(6);
  vctVec F0(3 * nSamps);
  //vctMat Fxi( 3,3 );
  vct3x3 Fxi;
  vctDynamicVector<vct3x3> Mi(nSamps);
  vctDynamicVector<vct3x3> Minv(nSamps);
  vctRot3 dR;
  vctRodRot3 dAlpha;
  vct3 dt;
  unsigned int idx;

  // J (cols 4 - 6)
  // set last three cols of J to negative block identity
  for (unsigned int i = 0; i<3 * nSamps; i += 3)
  {
    J.Element(i, 3) = -1.0;
    J.Element(i + 1, 4) = -1.0;
    J.Element(i + 2, 5) = -1.0;
  }

  // Step 1: Initialize transform params
  vct3x3 R = vct3x3::Eye();
  vct3 t(0.0);
  vctDynamicMatrixRef<double> Rref(R);

  int numIter = 0;
  double dTheta = 1.0;
  double dt_norm = 1.0;

  while (((dTheta > dTheta_term) || (dt_norm > dt_term)) && (numIter++ < maxIter))
  {
    // Step 2: F0
    vct3 res;
    vctDynamicVector<vct3> x_rot(nSamps);
    for (unsigned int i = 0; i < nSamps; i++)
    {
      x_rot.Element(i) = R*x.Element(i);
      res = y.Element(i) - (x_rot.Element(i) + t);
      // stack residuals into a single vector [f1x f1y f1z f2x...]'  
      idx = 3 * i;
      F0.Element(idx++) = res.X();
      F0.Element(idx++) = res.Y();
      F0.Element(idx) = res.Z();
    }

    // Step 3: J (cols 4 - 6) & Fx
    //   Note: Fx is block diagonal with all sub-blocks being equal
    //         to each other (equal to Fxi)
    vctDynamicMatrixRef<double> Jref3x3;
    for (unsigned int i = 0; i < nSamps; i++)
    {
      // populate first 3 cols of J
      Jref3x3.SetRef(J, 3 * i, 0, 3, 3);
      skew(x_rot.Element(i), Jref3x3);  // Jref3x3 = skew(x_rot.Element(i))
    }
    Fxi = -R;

    // Step 4: solve dP by least squares
    //  Note: Mi & Minv are block diagonal with same block
    //        repeated; therefore just calculate this sub-block
    for (unsigned int i = 0; i < nSamps; i++)
    {
      Mi[i] = Fxi*Mxi[i] * Fxi.TransposeRef() + Myi[i];
      ComputeCovInverse_NonIter(Mi[i], Minv[i]);

      //Mi[i] = Fxi*Mxi[i] * Fxi.TransposeRef() + Myi[i];
      //Minv[i] = Mi[i];
      //nmrInverse(Minv[i]); // computes inverse in-place
      ////nmrInverse(Minv, workspaceInv6);   // Anton (why doesn't this work?)
    }
    
    vctDynamicMatrixRef<double> JrowsRef;
    vctDynamicMatrixRef<double> Jt_MinvRef;
    vctDynamicMatrixRef<double> MinvRef;

    for (unsigned int i = 0; i < nSamps; i++)
    { // Minv is block diagonal => only multiply the sub-blocks
      //  Jt_Minv:  6x3N
      //  J:        3Nx6
      //  M:        3Nx3N  (block diagonal)
      MinvRef.SetRef(Minv[i]);
      idx = 3 * i;
      JrowsRef.SetRef(J, idx, 0, 3, 6);          // next 3 rows of J
      Jt_MinvRef.SetRef(Jt_Minv, 0, idx, 6, 3);  // next 3 cols of Jt_Minv
      Jt_MinvRef = JrowsRef.TransposeRef()*MinvRef;
    }
    Jt_Minv_J = Jt_Minv * J;
    neg_Jt_Minv_F0 = Jt_Minv * (-F0);

    // solve Ax = b
    //  Note: most of the run-time is for this SVD call
    try
    {
      A.Assign(Jt_Minv_J);
      // this changes the storage order of A (must use assign instead)
      //A = Jt_Minv_J;
      nmrSVD(A, U, S, Vt, workspaceSVD6x6);
    }
    catch (...)
    {
      std::cerr << "ERROR: Compute SVD failed" << std::endl;
      assert(0);
    }

    // dP = V * Sinv * U' * neg_Jt_Minv_F0
    vctVec Ut_F0(6);
    vctVec Sinv_Ut_F0(6);
    vctVec dP(6);
    Ut_F0 = U.TransposeRef()*neg_Jt_Minv_F0;
    for (unsigned int i = 0; i < 6; i++)
    {
      Sinv_Ut_F0.Element(i) = Ut_F0.Element(i) / S.Element(i);
    }
    dP = Vt.TransposeRef()*Sinv_Ut_F0;
    dAlpha.X() = dP.Element(0);
    dAlpha.Y() = dP.Element(1);
    dAlpha.Z() = dP.Element(2);
    dt.X() = dP.Element(3);
    dt.Y() = dP.Element(4);
    dt.Z() = dP.Element(5);

    // Step 5: [R,t]
    dR = vctRot3(dAlpha);
    R = dR*R;
    t = t + dt;

    //std::cout << "dt: " << dt[0] << " " << dt[1] << " " << dt[2] << "  dA: " << dAlpha[0] << " " << dAlpha[1]
    //<< " " << dAlpha[2] << "  " << dAlpha.Norm()*180/cmnPI << std::endl;

    dt_norm = dt.Norm();
    dTheta = dAlpha.Norm();
  }

  //vctRot3 R1( R );
  //vctAxAnRot3 Rrod( R1 );
  //std::cout << " TLS-Reg:  " << numIter << " iterations" << std::endl
  //  << "   t: " << t << std::endl
  //  << "   axis: " << Rrod.Axis() << " angle: " << Rrod.Angle()*180.0/cmnPI << std::endl;

  Fact.Rotation() = R;
  Fact.Translation() = t;
}


// --- Rotation --- //


// quaternion (Horn's) method of solving for rotation
vctRot3 SolveRotation_HornsMethod(vct3x3 &H)
{
  // Build Quaternion matrix N
  int i, j;
  vct4x4 N;

  N(0, 0) = H(0, 0) + H(1, 1) + H(2, 2);
  N(1, 0) = N(0, 1) = H(1, 2) - H(2, 1);
  N(2, 0) = N(0, 2) = H(2, 0) - H(0, 2);
  N(3, 0) = N(0, 3) = H(0, 1) - H(1, 0);
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      N(i + 1, j + 1) = H(i, j) + H(j, i);
    }
    N(i + 1, i + 1) -= N(0, 0);
  }

  // Compute SVD of N
  //   N = U*diag(S)*V'   where U = V
  static vctFixedSizeMatrix<double, 4, 4, VCT_COL_MAJOR> Ncopy;
  static vctFixedSizeMatrix<double, 4, 4, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 4, 4, VCT_COL_MAJOR> Vt;
  static vct4 S;
  static nmrSVDFixedSizeData<4, 4, VCT_COL_MAJOR>::VectorTypeWorkspace workspace;
  try
  {
    Ncopy.Assign(N);  // must use "assign" rather than equals to properly transfer between different vector orderings
    nmrSVD(Ncopy, U, S, Vt, workspace);
  }
  catch (...)
  {
    assert(0);
  }

  // optimal quaternion rotation is the eigenvector associated with the
  //  max eigenvalue of N:  N = W*diag(E)*W'
  //
  //  NOTE: singular values are always positive whereas eigenvalues can be both positive and negative
  //
  // compute eigenvalues from the SVD decomposition
  vct4 Eigen;
  for (int k = 0; k < 4; k++)
  {
    if (vctDotProduct(U.Column(k), Vt.Row(k)) > 0.0)
      Eigen[k] = S[k];
    else
      Eigen[k] = -S[k];
  }
  // find eigenvector associated with the largest eigenvalue
  vctQuaternionRotation3<double> rq;
  double maxEv = Eigen[0];
  int	maxEk = 0;
  for (int k = 1; k<4; k++)
  {
    if (Eigen[k]>maxEv) { maxEv = Eigen[k]; maxEk = k; }
  }
  rq.R() = U(0, maxEk);   // Do it this way to avoid problems from the choice of component order
  rq.X() = U(1, maxEk);
  rq.Y() = U(2, maxEk);
  rq.Z() = U(3, maxEk);

  return vctRot3(rq);
}

// SVD (Arun's) method of solving for rotation
vctRot3 SolveRotation_ArunsMethod(vct3x3 &H)
{
  // Compute SVD of H
  //  H = USV'
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Hcopy;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  static vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  vct3 S;
  try
  {
    Hcopy.Assign(H);  // must use "assign()" to transfer values properly between different vector orderings
    nmrSVD(Hcopy, U, S, Vt);
  }
  catch (...)
  {
    std::cerr << std::endl
      << "=====> ERROR: Compute SVD of H matrix failed in ComputeRotation_SVD()" << std::endl << std::endl;
    assert(0);
  }

  // Calculate R
  vctRot3 R;
  vct3x3 Ut;
  Ut.Assign(U.TransposeRef());
  vct3x3 V(3, 3);
  V.Assign(Vt.TransposeRef());
  R.ProductOf(V, Ut);

  // If determinant of R is -1, that means we have a reflection
  //   rather than rotation matrix. This can be fixed if at least
  //   one of the singular values ~= 0
  vct3x3 Dt = R;
  double Rdet = vctDeterminant<3>::Compute(Dt);
  if (Rdet < 0.0)
  {
    // negate the smallest singular value to flip reflection into a rotation matrix
    V.Column(2) = -V.Column(2);
    R.ProductOf(V, Ut);
  }

  return R;
}

// Compute the least-squares rotation matrix that minimizes the sum 
//  of point-to-point square distances.
//
//  Minimize:  Sum_i( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes an SVD solution for R following the method of Arun.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Arun et al, "Least-Squares Fitting of Two 3-D Point Sets", PAMI, 1987
//
void RotateP2P_LSQ_SVD(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R)
{
  int numPts = X.size();
  int i, j;
  vct3x3 tmp;
  vct3x3 H(0.0);

  for (i = 0; i < numPts; i++)
  {
    tmp.OuterProductOf(X[i], Y[i]);
    H += tmp;
  };

  R = SolveRotation_ArunsMethod(H);
}

// Compute the least-squares rotation matrix that minimizes the sum 
//  of point-to-point square distances.
//
//  Minimize:  Sum_i( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes an quaternion solution for R following the method of Horn.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Horn, "Closed-form solution of absolute orientation using unit quaternions", JOSA, 1987
//
void RotateP2P_LSQ_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R)
{
  int numPts = X.size();
  int i, j;
  vct3x3 tmp;
  vct3x3 H(0.0);

  for (i = 0; i < numPts; i++)
  {
    tmp.OuterProductOf(X[i], Y[i]);
    H += tmp;
  };

  H.Divide(numPts);  // normalization important?

  R = SolveRotation_HornsMethod(H);
}
