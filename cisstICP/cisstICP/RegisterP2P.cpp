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
  RotateP2P_LSQ_SVD( Xp, Yp, F.Rotation() );  
  //RotateP2P_LSQ_Quaternion(Xp, Yp, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Compute the least-squares rigid body transform that minimizes the sum
//  of weighted point-to-point distances.
//
//  Minimize:  Sum( Wi*||Yi - T*Xi||^2 )
//
//  This function follows the methods of Maurer:
//
//  Calvin R Maurer, Jr, et.al. "Registration of 3D Images Using Weighted Geometrical Features", IEEE TMI, 1996
//
void RegisterP2P_WLSQ(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctFrm3 &F)
{
  int rv;
  unsigned int N = X.size();
  assert(N == Y.size() && N == W.size());

  // weighted centroids of each point set
  vct3 Xmean = vctWeightedMean(X, W);
  vct3 Ymean = vctWeightedMean(Y, W);

  // recenter point sets about the centroids
  vctDynamicVector<vct3> Xp(N);
  vctDynamicVector<vct3> Yp(N);
  for (unsigned int i = 0; i < N; i++)
  {
    Xp[i] = X[i] - Xmean;
    Yp[i] = Y[i] - Ymean;
  }

  // Rotation
  //RotateP2P_WLSQ_SVD( Xp,Yp,W,F.Rotation() );
  RotateP2P_WLSQ_Quaternion(Xp, Yp, W, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Compute the least-squared rigid body transform to minimize weighted
//  P2P and normal distances assuming a vonMises-Fisher distribution
//  on the orientation noise and an independent, isotropic Gaussian 
//  distribution on the position noise.
//
//   Minimize:  -k*Sum_i( Wi*dot(Nyi,R*Nxi)) + B*Sum_i( Wi*||Yi - T*Xi||^2 )
//
//   Equivalently Minimizes:  
//        
//         k*Sum_i(1 - Wi*dot(Nyi,R*Nxi)) + B*Sum_i( Wi*||Yi - T*Xi||^2 )
//
void RegisterP2P_Normals_vMFG(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctFrm3 &F)
{
  int rv;
  unsigned int N = X.size();
  assert(N == Y.size() && N == Nx.size() && N == Ny.size() && N == W.size());

  // weighted centroids of each point set
  vct3 Xmean = vctWeightedMean(X, W);
  vct3 Ymean = vctWeightedMean(Y, W);

  // recenter point sets about the centroids
  vctDynamicVector<vct3> Xp(N);
  vctDynamicVector<vct3> Yp(N);
  for (unsigned int i = 0; i < N; i++)
  {
    Xp[i] = X[i] - Xmean;
    Yp[i] = Y[i] - Ymean;
  }

  // Rotation
  RotateP2P_Normals_vMFG_Quaternion(Xp, Yp, Nx, Ny, W, B, k, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Same as other function by this name, but without weights
//  (i.e. all weights set to 1)
void RegisterP2P_Normals_vMFG(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  double B, double k,
  vctFrm3 &F)
{
  int rv;
  unsigned int N = X.size();
  assert(N == Y.size() && N == Nx.size() && N == Ny.size());

  // weighted centroids of each point set
  vct3 Xmean = vctCentroid(X);
  vct3 Ymean = vctCentroid(Y);

  // recenter point sets about the centroids
  vctDynamicVector<vct3> Xp(N);
  vctDynamicVector<vct3> Yp(N);
  for (unsigned int i = 0; i < N; i++)
  {
    Xp[i] = X[i] - Xmean;
    Yp[i] = Y[i] - Ymean;
  }

  // Rotation
  RotateP2P_Normals_vMFG_Quaternion(Xp, Yp, Nx, Ny, B, k, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Compute the least-squared rigid body transform to minimize weighted
//  P2P and normal distances assuming an independent, isotropic 
//  Gaussian distribution on both the orientation and position noise.
//
//   Minimize:   normWeight*Sum_i( Wi*||Nyi - R*Nxi||^2 ) + posWeight*Sum_i( Wi*||Yi - (R*Xi+t)||^2 )
//
void RegisterP2P_Normals_VarEst(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double posWeight, double normWeight,
  vctFrm3 &F)
{
  int rv;
  unsigned int N = X.size();
  assert(N == Y.size() && N == Nx.size() && N == Ny.size() && N == W.size());

  // weighted centroids of each point set
  vct3 Xmean = vctWeightedMean(X, W);
  vct3 Ymean = vctWeightedMean(Y, W);

  // recenter point sets about the centroids
  vctDynamicVector<vct3> Xp(N);
  vctDynamicVector<vct3> Yp(N);
  for (unsigned int i = 0; i < N; i++)
  {
    Xp[i] = X[i] - Xmean;
    Yp[i] = Y[i] - Ymean;
  }

  // Rotation
  RotateP2P_Normals_VarEst_Quaternion(Xp, Yp, Nx, Ny, W, posWeight, normWeight, F.Rotation());

  // Translation
  //  optimal translation always aligns the centroid locations:  Ymean = R*Xmean + t
  //  =>  t = Ymean - R*Xmean
  F.Translation() = (Ymean - F.Rotation()*Xmean);
}

// Non-linear least squares Point-to-Point rigid body registration
//  based on Mahalanobis distance assuming a multivariate Gaussian
//  distribution on the match position errors and zero noise on
//  the sample position errors.
//
//    minimize:  Sum( (Yi - T*Xi)'*inv(My)*(Yi - T*Xi) )
//       where My is the covariance matrix for errors in Yi;
//       Xi is assumed to have no error
//
//  Note: an alternative to this function is to call RegisterP2P_TLS
//        and setting argument Mx = 0
//  Note: this function takes the inverted covariance matrix inv(My)
//        as input, whereas RegisterP2P_TLS takes the non-inverted 
//        covariance matrices as input
//
void RegisterP2P_LSQ_CovEst(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  vct3x3 Myinv, vctFrm3 &Fact)
{
  // Note: this follows the algorithm for total least squares with simplifications
  //       made due to Mxi = 0.

  // termination conditions
  double dTheta_term = 0.001;   // delta rotation (degrees)
  double dt_term = 0.001;       // delta translation (mm)
  int    maxIter = 30;        // max iterations

  // Workspace for SVD computations
  vctMat A(6, 6, VCT_COL_MAJOR);
  vctMat U(6, 6, VCT_COL_MAJOR);
  vctVec S(6);
  vctMat Vt(6, 6, VCT_COL_MAJOR);
  vctDynamicVector<double> workspaceSVD6x6(nmrSVDDynamicData::WorkspaceSize(A));

  unsigned int nSamps = x.size();

  vctMat J(3 * nSamps, 6, 0.0);
  vctMat Jt_Minv(6, 3 * nSamps);
  vctMat Jt_Minv_J(6, 6);
  vctVec neg_Jt_Minv_F0(6);
  vctVec F0(3 * nSamps);
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
  double dTheta = 100.0;
  double dt_norm = 100.0;

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

    // Step 3: J
    vctDynamicMatrixRef<double> Jref3x3;
    for (unsigned int i = 0; i < nSamps; i++)
    {
      // populate next rows in first 3 cols of J
      Jref3x3.SetRef(J, 3 * i, 0, 3, 3);
      skew(x_rot.Element(i), Jref3x3);  // Jref3x3 = skew(x_rot.Element(i))
    }

    // Step 4: solve dP by least squares
    vctDynamicMatrixRef<double> JrowsRef;
    vctDynamicMatrixRef<double> Jt_MinvRef;
    vctDynamicMatrixRef<double> MinvRef(Myinv);
    for (unsigned int i = 0; i < nSamps; i++)
    { // Myinv is block diagonal => only multiply the sub-blocks
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
      std::cout << "ERROR: Compute SVD failed" << std::endl;
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
      std::cout << "ERROR: Compute SVD failed" << std::endl;
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
      std::cout << "ERROR: Compute SVD failed" << std::endl;
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


// Non-linear Weighted Total Least Squares Point-to-Point rigid body registration
//  based on Mahalanobis distance assuming a multivariate Gaussian
//  distribution on both the match position errors and sample position errors
//  and assuming independence between the match and sample errors.
//  Solved via a linearized iterative scheme based on the Jacobian and
//  using Lagrange multipliers to enforce the model constraint.
//
//  Minimizes:  sum_i( Wi*dXi'*inv(Mxi)*dXi + Wi*dYi'*inv(Myi)*dYi )
//                where Yi = R*Xi + t  (model constraint)
//                      dX = Xcalc - Xobs
//                      dY = Ycalc - Yobs
//                      Wi = match weight
//
//   x   - observed sample values
//   y   - observed model values
//   Mxi - covariance matrix for dX  (assumed different for all dXi)
//   Myi - covariance matrix for dY  (assumed different for all dYi)
//   Wi  - vector of weights for each point pair
//
void RegisterP2P_WTLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3x3> &Mxi, const vctDynamicVector<vct3x3> &Myi,
  const vctDynamicVector<double> &Wi,
  vctFrm3 &Fact)
{
  // termination conditions
  double dTheta_term = 0.001;   // delta rotation (degrees)
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
    for (unsigned int i = 0; i < nSamps; i++)
    {
      Mi[i] = (Fxi*Mxi[i] * Fxi.TransposeRef() + Myi[i]) / Wi[i];
      Minv[i] = Mi[i];
      nmrInverse(Minv[i]); // computes inverse in-place
    }
    //nmrInverse(Minv, workspaceInv6);   // Anton (why doesn't this work?)

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
      std::cout << "ERROR: Compute SVD failed" << std::endl;
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

// Total Least Squares registration of point pair positions & orientations
//  (assumes error in both sample and model values)
//  solved via an iterative linear scheme based on Jacobian
void RegisterP2P_Normals_TLS()
{
  std::cout << "Method not implemented" << std::endl;
  exit(EXIT_FAILURE);
}


// Nonlinear registration of point pair positions
//  solved via an iterative linear scheme using Jacobian
void RegisterP2P_IterLSQ(const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  vctFrm3 &Fact,
  P2PIterativeCallback funCmptAb,
  vctDynamicVector<double> *pWorkspace)
{
  static const double dAlpha_Converged = 0.1*cmnPI / 180.0;
  static const double dt_Converged = 0.1;
  static const unsigned int maxIter = 10;

  unsigned int nSamps = x.size();
  vctFrm3 dF;
  vctRodRot3 dRrod;

  // storage for moving sample
  vctDynamicVector<vct3> x_mvg(nSamps);
  // initialize sample positions
  for (unsigned int i = 0; i < nSamps; i++)
  {
    x_mvg.Element(i) = Fact * x.Element(i);
  }

  // least squares variables
  //  A*g=b
  //  A = U*S*Vt
  vctMat A(3 * nSamps, 6, VCT_COL_MAJOR);
  vctVec b(3 * nSamps);
  vctVec gamma(6);
  vctMat U(3 * nSamps, 6, VCT_COL_MAJOR);
  vctVec S(6);
  vctMat Vt(6, 6, VCT_COL_MAJOR);
  vctVec temp(6);
  vctDynamicVectorRef<double> alpha(gamma, 0, 3); // first three parameters
  vctDynamicVectorRef<double> t(gamma, 3, 3);     // last three parameters

  //// SVD workspace
  //vctDynamicVector<double> workspaceTemp;
  //vctDynamicVectorRef<double> workspaceRef;
  //if (pWorkspace == 0)
  //{ // no workspace provided => create one
  //  workspaceTemp.SetSize(nmrSVDEconomyDynamicData::WorkspaceSize(A));
  //  workspaceRef.SetRef(workspaceTemp);
  //}
  //else
  //{ // use workspace provided
  //  workspaceRef.SetRef(*pWorkspace);
  //}

  // Run iterated linear least squares
  unsigned int iter = 0;
  while (iter++ < maxIter)
  {
    // compute A & b for current sample positions
    funCmptAb(x_mvg, y, A, b);

    // SVD decomposition of A
    try
    {
      nmrSVDEconomy(A, U, S, Vt);// workspaceRef);
    }
    catch (...)
    {
      std::cout << "ERROR: Compute SVD Economy failed" << std::endl;
    }

    // compute least squares parameter estimate
    //  Ax=b
    //  A = U*S*Vt
    //  x = V*inv(S)*Ut*b
    //
    temp.ProductOf(U.TransposeRef(), b);
    temp.ElementwiseDivide(S);
    gamma.ProductOf(Vt.TransposeRef(), temp);

    // update coordinate transfrom
    dRrod.Assign(alpha);
    dF.Rotation() = vctRot3(dRrod);
    dF.Translation() = t;
    Fact = dF * Fact;

    // update sample positions
    for (unsigned int i = 0; i < nSamps; i++)
    {
      x_mvg.Element(i) = Fact * x.Element(i);
    }

    // check for convergence
    double dAlpha = alpha.Norm();
    double dt = t.Norm();
    if (dAlpha < dAlpha_Converged && dt < dt_Converged)
    {
      //std::cout << "...P2P converged in " << iter << " iterations" << std::endl;
      break;
    }
  }

  if (iter == maxIter)
  {
    std::cout << "WARNING: maxIter (" << maxIter << ") reached for iterative P2P registration" << std::endl;
  }
}

// Nonlinear registration of point pair positions & orientations
//  solved via an iterative linear scheme using Jacobian
void RegisterP2P_Normals_IterLSQ(const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3> &Nx, const vctDynamicVector<vct3> &Ny,
  vctFrm3 &Fact,
  Callback_P2PNormals_IterLSQ funCmptAb,
  vctDynamicVector<double> *pWorkspace)
{
  static const double dAlpha_Converged = 0.05*cmnPI / 180.0;
  static const double dt_Converged = 0.05;
  static const unsigned int maxIter = 50;

  unsigned int nSamps = x.size();
  vctFrm3 dF;
  vctRodRot3 dRrod;

  // storage for moving sample
  vctDynamicVector<vct3> x_mvg(nSamps);
  vctDynamicVector<vct3> Nx_mvg(nSamps);
  // initialize sample positions
  for (unsigned int i = 0; i < nSamps; i++)
  {
    x_mvg.Element(i) = Fact * x.Element(i);
    Nx_mvg.Element(i) = Fact.Rotation() * Nx.Element(i);
  }

  // least squares variables
  //  A*g=b
  //  A = U*S*Vt
  vctMat A(6 * nSamps, 6, VCT_COL_MAJOR);
  vctVec b(6 * nSamps);
  vctVec g(6);
  vctMat U(6 * nSamps, 6, VCT_COL_MAJOR);
  vctVec S(6);
  vctMat Vt(6, 6, VCT_COL_MAJOR);
  vctVec temp(6);
  vctDynamicVectorRef<double> alpha(g, 0, 3); // first three parameters
  vctDynamicVectorRef<double> t(g, 3, 3);     // last three parameters

  //// SVD workspace
  //vctDynamicVector<double> workspaceTemp;
  //vctDynamicVectorRef<double> workspaceRef;
  //if (pWorkspace == 0)
  //{ // no workspace provided => create one
  //  workspaceTemp.SetSize(nmrSVDEconomyDynamicData::WorkspaceSize(A));
  //  workspaceRef.SetRef(workspaceTemp);
  //}
  //else
  //{ // use workspace provided
  //  workspaceRef.SetRef(*pWorkspace);
  //}

  // Run iterated linear least squares
  unsigned int iter = 0;
  while (iter++ < maxIter)
  {
    // compute A & b for current sample positions
    funCmptAb(x_mvg, y, Nx_mvg, Ny, A, b);

    // SVD decomposition of A
    try
    {
      nmrSVDEconomy(A, U, S, Vt);// workspaceRef);
    }
    catch (...)
    {
      std::cout << "ERROR: Compute SVD Economy failed" << std::endl;
    }

    // compute least squares parameter estimate
    //  Ax=b
    //  A = U*S*Vt
    //  x = V*inv(S)*Ut*b
    //
    temp.ProductOf(U.TransposeRef(), b);
    temp.ElementwiseDivide(S);
    g.ProductOf(Vt.TransposeRef(), temp);

    // update coordinate transfrom
    dRrod.Assign(alpha);
    dF.Rotation() = vctRot3(dRrod);
    dF.Translation() = t;
    Fact = dF * Fact;

    // update sample positions
    for (unsigned int i = 0; i < nSamps; i++)
    {
      x_mvg.Element(i) = Fact * x.Element(i);
      Nx_mvg.Element(i) = Fact.Rotation() * Nx.Element(i);
    }

    // check for convergence
    double dAlpha = alpha.Norm();
    double dt = t.Norm();
    if (dAlpha < dAlpha_Converged && dt < dt_Converged)
    {
      //std::cout << "...P2P converged in " << iter << " iterations" << std::endl;
      break;
    }
  }

  if (iter == maxIter)
  {
    std::cout << "WARNING: maxIter (" << maxIter << ") reached for iterative P2P registration" << std::endl;
  }
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

  //vct4x4 V;
  //vct4 Eigen;
  //int rc = nmrJacobi(N, Eigen, V);
  //double maxEv = Eigen[0];
  //int	maxEk = 0;
  //for (int k = 1; k<4; k++)
  //{
  //  if (Eigen[k]>maxEv) { maxEv = Eigen[k]; maxEk = k; }
  //};
  //vctQuaternionRotation3<double> rq;
  //rq.R() = V(0, maxEk);   // Do it this way to avoid problems from the choice of component order
  //rq.X() = V(1, maxEk);
  //rq.Y() = V(2, maxEk);
  //rq.Z() = V(3, maxEk);
  //std::cout << "rq (nmrJacobi): " << rq << std::endl;  
  //std::cout << " V: " << std::endl << V << std::endl;
  //std::cout << " Ev: " << std::endl << Eigen << std::endl;
  //std::cout << " maxEv: " << maxEv << std::endl;
  //std::cout << " maxEk: " << maxEk << std::endl;
  //return  vctRot3(rq);
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
    std::cout << std::endl
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
    //// check for smallest singular value ~= 0
    ////  nmrSVD returns singular values in descending order
    //if (S.Element(2) > 0.001 && S.Element(2) / S.Element(0) > 0.001)
    //{
    //  std::cout << std::endl
    //    << "=====> WARNING: Arun's method failed; determinant of R is -1 and there exist no singular value ~= 0" << std::endl
    //    << "                Largetst singular value: " << S.Element(0) << " Smallest singular value : " << S.Element(2) << std::endl << std::endl;
    //}

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

  //H.Divide(numPts);  // normalization important?

  R = SolveRotation_HornsMethod(H);
}

// Compute the least-squares rotation matrix that minimizes the weighted sum 
//  of point-to-point distances.
//
//  Minimize:  Sum_i( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes a quaternion solution for R from the sum of 3x3 covariance matrices
//   following the method of Horn.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Horn, "Closed-form solution of absolute orientation using unit quaternions", JOSA, 1987
//
void RotateP2P_WLSQ_Quaternion(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctRot3 &R)
{
  unsigned int nSamps = X.size();
  assert(nSamps == Y.size() && nSamps == W.size());

  // Build the covariance matrix
  //  H = iSum( Wi*(Xi*Bi') )
  vct3x3 H(0.0);
  vct3x3 Htemp;
  for (unsigned int i = 0; i < nSamps; i++)
  {
    Htemp.OuterProductOf(W[i] * X[i], Y[i]);
    H.Add(Htemp);
  }

  H.Divide(nSamps);   // is this desirable? (changes the soln very minutely)

  R = SolveRotation_HornsMethod(H);
};

// Compute the least-squares rotation matrix that minimizes the weighted sum 
//  of point-to-point distances.
//
//  Minimize:  Sum_i( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes an SVD solution for R from the sum of 3x3 covariance matrices
//   following the methods of Arun & Maurer.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Arun, et.al., "Least-squares fitting of two 3-D point sets", PAMI, 1987
//  Calvin R Maurer, Jr, et.al. "Registration of 3D Images Using Weighted Geometrical Features", IEEE TMI, 1996
//
void RotateP2P_WLSQ_SVD(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctRot3 &R)
{
  unsigned int N = X.size();
  assert(N == Y.size() && N == W.size());

  // Build the covariance matrix
  //  H = iSum( Wi*(Xi*Bi') )
  vct3x3 H(0.0);
  vct3x3 Htemp;
  for (unsigned int i = 0; i < N; i++)
  {
    Htemp.OuterProductOf(W[i] * X[i], Y[i]);
    H.Add(Htemp);
  }

  R = SolveRotation_ArunsMethod(H);
};


// Compute the least-squares rotation matrix that minimizes the weighted sum 
//  of point-to-point and normal distances.
//
//  Minimize:  -k*Sum( Wi*dot(Nyi,R*Nxi) ) + B*Sum( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  k*Sum( Wi*dot(Nyi,R*Nxi) ) + 2*B*Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes R using an SVD solution, extending the methods of Arun & Maurer.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Arun, et.al., "Least-squares fitting of two 3-D point sets", PAMI, 1987
//  Calvin R Maurer, Jr, et.al. "Registration of 3D Images Using Weighted Geometrical Features", IEEE TMI, 1996
//
void RotateP2P_Normals_vMFG_SVD(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctRot3 &R)
{
  unsigned int N = X.size();
  assert(N == Y.size() && N == Nx.size() && N == Ny.size() && N == W.size());

  // Build the covariance matrices
  vct3x3 H(0.0), H1(0.0), H2(0.0);
  vct3x3 Htemp;
  //  points:  H1 = Sum( Wi*(Ai*Bi') )
  for (unsigned int i = 0; i < N; i++)
  {
    Htemp.OuterProductOf(W[i] * X[i], Y[i]);
    H1.Add(Htemp);
  }
  //  normals:  H2 = Sum( Wi*(Nxi*Nyi') )
  for (unsigned int i = 0; i < N; i++)
  {
    Htemp.OuterProductOf(W[i] * Nx[i], Ny[i]);
    H2.Add(Htemp);
  }
  H = H1.Multiply(2.0*B) + H2.Multiply(k);

  R = SolveRotation_ArunsMethod(H);
}


// Compute the least-squares rotation matrix that minimizes the weighted sum 
//  of point-to-point and normal distances.
//
//  Minimize:  -k*Sum( Wi*dot(Nyi,R*Nxi) ) + B*Sum( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  k*Sum( Wi*dot(Nyi,R*Nxi) ) + 2*B*Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes R using a quaternion solution, extending Horn's method.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Horn, "Closed-form solution of absolute orientation using unit quaternions", JOSA, 1987
//
void RotateP2P_Normals_vMFG_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctRot3 &R)
{
  unsigned int numPts = X.size();
  assert(numPts == Y.size() && numPts == Nx.size() && numPts == Ny.size() && numPts == W.size());

  // Build the covariance matrices
  //  Note: H is M in Horn's paper
  vct3x3 H(0.0), H1(0.0), H2(0.0);
  vct3x3 Htemp;
  //  points:  H1 = Sum( Wi*(Ai*Bi') )
  for (unsigned int i = 0; i < numPts; i++)
  {
    Htemp.OuterProductOf(W[i] * X[i], Y[i]);
    H1.Add(Htemp);
  }
  //  normals:  H2 = Sum( Wi*(Nxi*Nyi') )
  for (unsigned int i = 0; i < numPts; i++)
  {
    Htemp.OuterProductOf(W[i] * Nx[i], Ny[i]);
    H2.Add(Htemp);
  }
  H = H1.Multiply(2.0*B) + H2.Multiply(k);

  R = SolveRotation_HornsMethod(H);
}

// same as above but without weights
void RotateP2P_Normals_vMFG_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  double B, double k,
  vctRot3 &R)
{
  RotateP2P_Normals_vMFG_Quaternion(X, Y, Nx, Ny, vctDoubleVec(X.size(), 1.0), B, k, R);
}


// Compute the least-squares rotation matrix that minimizes the weighted sum 
//  of point-to-point and normal distances.
//
//  Minimize:  normWeight*Sum( Wi*||Nyi - R*Nxi||^2 ) + posWeight*Sum( Wi*||Yi - R*Xi||^2 )
//
//      => Maximize:  normWeight*Sum( Wi*dot(Nyi,R*Nxi) ) + posWeight*Sum( Wi*dot(Yi,R*Xi) )
//
//  Computes R using a quaternion solution, extending Horn's method.
//
//  NOTE: This method assumes that X and Y have already been
//        centered about their respective centroids. => the values
//        of Xi & Yi lead directly to their covariance matrices without
//        having to subtract off the mean values.
//
//  Horn, "Closed-form solution of absolute orientation using unit quaternions", JOSA, 1987
//
void RotateP2P_Normals_VarEst_Quaternion(const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double posWeight, double normWeight,
  vctRot3 &R)
{
  unsigned int numPts = X.size();
  assert(numPts == Y.size() && numPts == Nx.size() && numPts == Ny.size() && numPts == W.size());

  // Build the covariance matrices
  //  Note: H is M in Horn's paper
  vct3x3 H(0.0), H1(0.0), H2(0.0);
  vct3x3 Htemp;
  //  points:  H1 = Sum( Wi*(Ai*Bi') )
  for (unsigned int i = 0; i < numPts; i++)
  {
    Htemp.OuterProductOf(W[i] * X[i], Y[i]);
    H1.Add(Htemp);
  }
  //  normals:  H2 = Sum( Wi*(Nxi*Nyi') )
  for (unsigned int i = 0; i < numPts; i++)
  {
    Htemp.OuterProductOf(W[i] * Nx[i], Ny[i]);
    H2.Add(Htemp);
  }
  H = H1.Multiply(posWeight) + H2.Multiply(normWeight);

  SolveRotation_ArunsMethod(H);
}
