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
#ifndef _Utility_H
#define _Utility_H

#include "cisstVector.h"
#include "cisstOSAbstraction.h"

#include "cisstMesh.h"
#include "cisstICP.h"

void CreateDir(const std::string &dir);

void transform_write(vctFrm3 &t, std::string &filename);

vctDynamicVector<vctRot3> rotations_read(std::string &filepath);
vctDynamicVector<vct3x3> cov_read(std::string &filepath);
void cov_write(vctDynamicVector<vct3x3> &cov, std::string &filename);

double ExtractGaussianRVFromStream(std::ifstream &randnStream);

// This function uses the Box-Muller method to generate
//  zero mean, unit variance Gaussian random variables
//  from a source of uniform random variables
//double GenerateGaussianRV( cmnRandomSequence &rvGen );

void ComputeGroundTruthStats(
  const vctDynamicVector<vct3> &samples,
  const vctDynamicVector<vct3> &sampleNorms,
  const vctFrm3 &Fgt, const vctFrm3 &Freg,
  double &posErr_mean, double &posErr_SD,
  double &normErr_mean, double &normErr_SD);

void GenerateRandomTransform(
  unsigned int randSeed, unsigned int &randSeqPos,
  double minOffsetPos, double maxOffsetPos,
  double minOffsetAng, double maxOffsetAng,
  vctFrm3 &F);

void GenerateRandomRotation(
  unsigned int randSeed, unsigned int &randSeqPos,
  double minOffsetAng, double maxOffsetAng,
  vctRot3 &R);

void GenerateRandomL(
  unsigned int randSeed, unsigned int &randSeqPos,
  const vctDynamicVector<vct3> &uDir, vctDynamicVector<vct3x2> &L);

void CreateMesh(
  cisstMesh &mesh,
  const std::string &meshLoadPath,
  std::string *SavePath_Mesh = 0);

void GenerateNoisyMesh(
  cisstMesh &mesh,
  cisstMesh &noisyMesh,
  std::ifstream &randnStream,
  //unsigned int randSeed, unsigned int &randSeqPos
  double noiseStdDev,
  std::string *SavePath_NoisyMesh);


// Defines the measurement noise for each triangle in a mesh
//  (not the measurement noise of samples, but of the mesh itself)
void SetMeshTriangleCovariances(
  cisstMesh &mesh,
  double stdDevPerpPlane, double stdDevInPlane);

//  cisstRandomSeq - source for uniform distributed random variables
void GenerateSamples(
  cisstMesh &mesh,
  unsigned int randSeed, unsigned int &randSeqPos,
  unsigned int nSamps,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<unsigned int>  &sampleDatums,
  std::string *SavePath_Samples = 0);

// Generate noisy samples having the specified Gaussian distributions
void GenerateNoisySamples_Gaussian(
  std::ifstream &randnStream,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3>   &noisySamples,
  std::string *SavePath_NoisySamples = 0,
  std::string *SavePath_Cov = 0);

void GenerateOutlierSamples_SurfaceOffset(
  unsigned int randSeed, unsigned int &randSeqPos,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &outlierSamples,
  double minPosOffsetOutlier, double maxPosOffsetOutlier,
  std::string *SavePath_OutlierSamples = 0);

void ComputeCovariances_Random(
  unsigned int randSeed, unsigned int &randSeqPos,
  const vct3 &StdDev,
  vctDynamicVector<vct3x3> &cov,
  unsigned int nCov,
  std::string *SavePath_Cov = 0);

void ComputeCovariances_SurfaceModel(
  double StdDevInPlane, double StdDevPerpPlane,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3x3> &modelCov,
  //vctDynamicVector<vct3x3> &modelInvCov,
  std::string *SavePath_Cov = 0);

void GenerateSampleErrors_Covariance(
  std::ifstream &randnStream,
  const vctDynamicVector<vct3x3> &sampleCov,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &noisySamples,
  std::string *SavePath_NoisySamples);

//// Generate sample noise having a random position covariance with the
////  specified eigenvalues and random orientation with the specified
////  circular standard deviation and eccentricity.
//void GenerateSampleRandomNoise(
//  unsigned int randSeed, unsigned int &randSeqPos,
//  std::ifstream &randnStream,
//  vct3 &covEigenvalues,
//  vctDynamicVector<vct3>   &samples,
//  vctDynamicVector<vct3>   &noisySamples,
//  vctDynamicVector<vct3x3> &sampleCov,
//  vctDynamicVector<vct3x3> &sampleInvCov,
//  double percentOutliers,
//  double minPosOffsetOutlier, double maxPosOffsetOutlier,
//  std::string *SavePath_NoisySamples = 0,
//  std::string *SavePath_Cov = 0);


// Generate sample noise having the specified standard deviation of
//  noise in directions parallel and perpendicular to the triangle
//  plane from which the sample was drawn.
void GenerateSampleErrors_SurfaceNoise(
  unsigned int randSeed, unsigned int &randSeqPos,
  std::ifstream &randnStream,
  double StdDevInPlane, double StdDevPerpPlane,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &noisySamples,
  vctDynamicVector<vct3x3> &sampleCov,
  //vctDynamicVector<vct3x3> &sampleInvCov,
  double percentOutliers,
  double minPosOffsetOutlier, double maxPosOffsetOutlier,
  std::string *SavePath_NoisySamples = 0,
  std::string *SavePath_Cov = 0);

// Generate oriented sample noise having a random position covariance with the
//  specified eigenvalues and random orientation with the specified
//  circular standard deviation and eccentricity.
void GenerateSampleRandomNoise(
  unsigned int randSeed, unsigned int &randSeqPos,
  std::ifstream &randnStream,
  vct3 &covEigenvalues,
  double circStdDev, double circEccentricity,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &noisySamples,
  vctDynamicVector<vct3>   &noisySampleNorms,
  vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3x3> &sampleInvCov,
  vctDynamicVector<vct3x2> &noiseL,
  double percentOutliers,
  double minPosOffsetOutlier, double maxPosOffsetOutlier,
  double minAngOffsetOutlier, double maxAngOffsetOutlier,
  std::string *SavePath_NoisySamples = 0,
  std::string *SavePath_Cov = 0,
  std::string *savePath_L = 0);

// Generate oriented sample noise having the specified standard deviation of
//  noise in directions parallel and perpendicular to the triangle
//  plane from which the sample was drawn.
void GenerateSampleSurfaceNoise(
  unsigned int randSeed, unsigned int &randSeqPos,
  std::ifstream &randnStream,
  double StdDevInPlane, double StdDevPerpPlane,
  double circStdDev, double circEccentricity,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &noisySamples,
  vctDynamicVector<vct3>   &noisySampleNorms,
  vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3x3> &sampleInvCov,
  vctDynamicVector<vct3x2> &noiseL,
  double percentOutliers,
  double minPosOffsetOutlier, double maxPosOffsetOutlier,
  double minAngOffsetOutlier, double maxAngOffsetOutlier,
  std::string *SavePath_NoisySamples = 0,
  std::string *SavePath_Cov = 0,
  std::string *SavePath_L = 0);

void Draw3DGaussianSample(
  std::ifstream &randnStream,
  const vct3x3 &M, vct3 &x);

void DrawFisherSample(
  std::ifstream &randnStream,
  double k, const vct3 &mean, vct3 &n);

void DrawGIMLOPSample(
  std::ifstream &randnStream,
  double k, double B, const vct3 &mean,
  const vctFixedSizeMatrix<double, 3, 2> &L, vct3 &n);

void Callback_SaveIterationsToFile_Utility(cisstICP::CallbackArg &arg, void *userData);
void Callback_TrackRegPath_Utility(cisstICP::CallbackArg &arg, void *userData);

void WriteToFile_Cov(
  const vctDynamicVector<vct3x3> &cov,
  std::string &filePath);

void WriteToFile_L(
  const vctDynamicVector<vctFixedSizeMatrix<double, 3, 2> > &L,
  std::string &filePath);

vctRot3 XProdRotation(const vct3 &a, const vct3 &b);


double  ComputeAvgNeighborDistance(std::string meshFile);
double  ComputeAvgNeighborDistance(cisstMesh mesh);
bool    TrianglesAreNeighbors(vctFixedSizeVector<int, 3> vi, vctFixedSizeVector<int, 3> vj);
bool    vct3HasValue(int x, vctFixedSizeVector<int, 3> vj);

#endif
