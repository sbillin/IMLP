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


// Compute the centroid for a set of vetors
vct3 vctCentroid(const vctDynamicVector<vct3>& A);
// Compute the weighted centroid for a set of vectors
vct3 vctWeightedMean(const vctDynamicVector<vct3>& A, const vctDoubleVec &W);


// --- P2P Rigid Body Registration --- //

// Standard linear least squares Point-to-Point rigid body registration
//  based on Euclidean distance.
//    minimize:  Sum( ||Yi - T*Xi||^2 )
void RegisterP2P_LSQ(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctFrm3 &F);

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


// --- P2P Rotations --- //

// internal routines
vctRot3 SolveRotation_HornsMethod(vct3x3 &H);
vctRot3 SolveRotation_ArunsMethod(vct3x3 &H);

// Compute least squares point-to-point rotation using quaternion (Horn's) method
void RotateP2P_LSQ_Quaternion(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R);

// Compute least squares point-to-point rotation using SVD (Arun's) method
void RotateP2P_LSQ_SVD(
  const vctDynamicVector<vct3> &X,
  const vctDynamicVector<vct3> &Y,
  vctRot3 &R);

#endif // _RegisterP2P_H
