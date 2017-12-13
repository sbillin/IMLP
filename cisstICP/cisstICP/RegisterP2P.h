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
#ifndef _RegisterP2P_H
#define _RegisterP2P_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include <cisstNumerical.h>
#include <cisstVector.h>


// --- P2P Rigid Body Registration --- //

// Standard linear least squares Point-to-Point rigid body registration
//  based on Euclidean distance.
//    minimize:  Sum( ||Yi - T*Xi||^2 )
void RegisterP2P_LSQ(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctFrm3 &F);

// Linear weighted least squares Point-to-Point rigid body registration
//  based on weighted Euclidean distance.
//    minimize:  Sum( Wi*||Yi - T*Xi||^2 )
void RegisterP2P_WLSQ(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctFrm3 &F);

// Linear weighted least squares Point-to-Point rigid body registration
//   for oriented points assuming a vonMises-Fisher & independent, isotropic 
//   Gaussian noise model for orientations and positions respectively.
// minimize:              -k*Sum_i( Wi*dot(Nyi,R*Nxi)) + B*Sum_i( Wi*||Yi - T*Xi||^2 )
// equivalently minimize:  k*Sum_i(1 - Wi*dot(Nyi,R*Nxi)) + B*Sum_i( Wi*||Yi - T*Xi||^2 )
void RegisterP2P_Normals_vMFG(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctFrm3 &F);
void RegisterP2P_Normals_vMFG(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  double B, double k,
  vctFrm3 &F);

// Linear weighted least squares Point-to-Point rigid body registration
//   for oriented points using independent, isotropic Gaussian noise model.
// minimize:   normWeight*Sum_i( Wi*||Nyi - R*Nxi||^2 ) + posWeight*Sum_i( Wi*||Yi - (R*Xi+t)||^2 )
void RegisterP2P_Normals_VarEst(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double posWeight, double normWeight,
  vctFrm3 &F);

// Non-linear least squares Point-to-Point rigid body registration
//  based on Mahalanobis distance.
//    minimize:  Sum( (Yi - T*Xi)'*inv(M)*(Yi - T*Xi) )
//       where M is the covariance matrix for errors in Yi
//       (Xi is assumed to have no error)
void RegisterP2P_LSQ_CovEst(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  vct3x3 Myinv, vctFrm3 &Fact);


// Total Least Squares Point-to-Point rigid body registration
//  (assumes error in both sample and model values)
//  solved via an iterative linear scheme based on Jacobian and Lagrange multipliers
//
//  Minimizes:  sum( dXi'*Mxi*dXi + dYi'*Myi*dYi )
//                where Yi = R*Xi + t  (model constraint)
//                      dX = Xcalc - Xobs, dY = Ycalc - Yobs
//
//   x   - observed sample values
//   y   - observed model values
//   Mxi - covariance matrix for dX  (assumed same for all dXi)
//   Myi - covariance matrix for dY  (assumed same for all dYi)
//
void RegisterP2P_TLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vct3x3 &Mxi, const vct3x3 &Myi,
  vctFrm3 &Fact);
//  ...
//  
//   Mxi - covariance matrix for dX  (not same for all dXi)
//   Myi - covariance matrix for dY  (not same for all dYi)
//   Wi  - vector of weighs for each point pair (optional)
//
//  Note: this version requires 50% more run-time than above version where noise
//        model is the same for every point
//
void RegisterP2P_TLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3x3> &Mxi, const vctDynamicVector<vct3x3> &Myi,
  vctFrm3 &Fact);
void RegisterP2P_WTLS(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3x3> &Mxi, const vctDynamicVector<vct3x3> &Myi,
  const vctDynamicVector<double> &Wi,
  vctFrm3 &Fact);
void RegisterP2P_Normals_TLS();


// Non-linear least squares Point-to-Point registration by linearizing the cost function
//  and computing the linear least squares estimate of transform parmaters x
//  from Ax=b at each iteration until convergence.
// Note: this function may be used to solve a non-linear (i.e., Mahalanobis distance) least
//       squares problem. Alternatively, a Mahalanobis distance may be transformed to Euclidean
//       space and solved using standard linear least squares.
//
// Type for callback function which computes A & b at each iteration
//  where Ax=b and x is the transform parameters [alpha, t]
//  and alpha is a Rodrigues vector representation of the rotation
typedef  void(*P2PIterativeCallback)(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  vctMat &A, vctVec &B);
void RegisterP2P_IterLSQ(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  vctFrm3 &Fact,
  P2PIterativeCallback funCmptAb,
  vctDynamicVector<double> *pWorkspace = 0);

// Non-linear least squares Point-to-Point registration for oriented points solved by linearizing the 
//  cost function and computing the linear least squares estimate of transform parmaters x
//  from Ax=b at each iteration until convergence.
// Note: this function may be used to solve a non-linear (i.e., Mahalanobis distance) least
//       squares problem. Alternatively, a Mahalanobis distance may be transformed to Euclidean
//       space and solved using standard linear least squares.
//
// Type for callback function which computes A & b at each iteration
//  where Ax=b and x is the transform parameters [alpha, t]
//  and alpha is a Rodrigues vector representation of the rotation
typedef  void(*Callback_P2PNormals_IterLSQ)(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3> &Nx, const vctDynamicVector<vct3> &Ny,
  vctMat &A, vctVec &b);
void RegisterP2P_Normals_IterLSQ(
  const vctDynamicVector<vct3> &x, const vctDynamicVector<vct3> &y,
  const vctDynamicVector<vct3> &Nx, const vctDynamicVector<vct3> &Ny,
  vctFrm3 &Fact,
  Callback_P2PNormals_IterLSQ funCmptAb,
  vctDynamicVector<double> *pWorkspace = 0);

// --- P2P Rotations --- //

// internal routines
vctRot3 SolveRotation_HornsMethod(vct3x3 &H);
vctRot3 SolveRotation_ArunsMethod(vct3x3 &H);

// uses Horn's method
void RotateP2P_LSQ_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R);

// uses Arun's method
void RotateP2P_LSQ_SVD(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R);

// Compute weighted least squares point-to-point rotation using quaternion (Arun's) method
void RotateP2P_WLSQ_SVD(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctRot3 &R);

// Compute weighted least squares point-to-point rotation using quaternion (Horn's) method
void RotateP2P_WLSQ_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDoubleVec &W,
  vctRot3 &R);

// Compute weighted least squares point-to-point rotation for oriented points
//  using SVD method (extension to Arun's method) 
//  assuming a vonMises-Fisher / independent, isotropic Gaussian noise model.
void RotateP2P_Normals_vMFG_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctRot3 &R);
void RotateP2P_Normals_vMFG_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  double B, double k,
  vctRot3 &R);

// Compute weighted least squares point-to-point rotation for oriented points
//  using SVD method (extension to Arun's method)
//  assuming a vonMises-Fisher / independent, isotropic Gaussian noise model
void RotateP2P_Normals_vMFG_SVD(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double B, double k,
  vctRot3 &R);

// Compute weighted least squares point-to-point rotation for oriented points
//  using SVD method (extension to Arun's method)
//  assuming a independent, isotropic Gaussian noise model.
// minimize:  normWeight*Sum( Wi*||Nyi - R*Nxi||^2 ) + posWeight*Sum( Wi*||Yi - R*Xi||^2 )
void RotateP2P_Normals_VarEst_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  const vctDynamicVector<vct3> &Nx,
  const vctDynamicVector<vct3> &Ny,
  const vctDoubleVec &W,
  double posWeight, double normWeight,
  vctRot3 &R);

#endif // _RegisterP2P_H
