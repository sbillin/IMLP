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

#include <fstream>

//#include <boost/filesystem.hpp>

#include "cisstNumerical.h"

#include "utility.h"
#include "cisstPointCloud.h"
#include "utilities.h"


//void CreateDir(const std::string &dir)
//{
//  boost::filesystem::create_directories(dir.c_str());
//}

//#include <windows.h>
//void CreateDir(const std::string &dir)
//{
//  int rv = CreateDirectory(dir.c_str(), NULL);
//  if (rv || ERROR_ALREADY_EXISTS == GetLastError())
//  {
//    if (rv)
//    { // created directory
//      std::cout << "Creating directory: \"" << dir << "\"" << std::endl;
//    }
//    // else directory already exists
//  }
//  else
//  { // failed to create directory
//    std::cout << "ERROR! create directory failed: " << dir << std::endl;
//    assert(0);
//  }
//}

void transform_write(vctFrm3 &F, std::string &filename)
{
  // output format:
  //  r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz  
  std::ofstream out(filename.c_str());
  if (out) {
    out << F.Rotation().Row(0) << " " << F.Rotation().Row(1) << " "
      << F.Rotation().Row(2) << " " << F.Translation();
    //out << t.Rotation().Row(0) << "  " << t.Translation()[0] << std::endl 
    //  << t.Rotation().Row(1) << "  " << t.Translation()[1] << std::endl
    //  << t.Rotation().Row(2) << "  " << t.Translation()[2] << std::endl;
  }
  else {
    out << "ERROR: cannot open transform file for writing" << std::endl;
    assert(0);
  }
  out.close();
}


vctDynamicVector<vctRot3> rotations_read(std::string &filepath)
{
  unsigned int itemsRead;
  std::string line;
  float f1, f2, f3, f4, f5, f6, f7, f8, f9;

  std::cout << "Reading vctRot3 matrices from file: " << filepath << std::endl;
  std::ifstream fs(filepath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filepath << std::endl;
    assert(0);
  }

  // read matrices
  unsigned int numRot = 0;
  std::vector<vctRot3> rotVct;
  vctRot3 rot;
  while (fs.good())
  {
    std::getline(fs, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f %f %f %f %f %f %f",
      &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9);
    if (itemsRead != 9)
    {
      break;
      //std::cerr << "ERROR: expeced a covariance entry at line: " << line << std::endl;
      //assert(0);
    }
    rot.Assign(f1, f2, f3, f4, f5, f6, f7, f8, f9);
    rotVct.push_back(rot);
    numRot++;
  }
  //if (fs.bad() || fs.fail())
  //{
  //  std::cerr << "ERROR: read pts from file failed; last line read: " << line << std::endl;
  //  assert(0);
  //}
  fs.close();

  vctDynamicVector<vctRot3> rotArray(numRot);
  for (unsigned int i = 0; i < numRot; i++)
  {
    rotArray[i] = rotVct[i];
  }

  return rotArray;
}

vctDynamicVector<vct3x3> cov_read(std::string &filepath)
{
  unsigned int itemsRead;
  std::string line;
  float f1, f2, f3, f4, f5, f6, f7, f8, f9;

  std::cout << "Reading 3x3 matrices from file: " << filepath << std::endl;
  std::ifstream fs(filepath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filepath << std::endl;
    assert(0);
  }

  // read matrices
  unsigned int numCov = 0;
  std::vector<vct3x3> covVct;
  vct3x3 cov;
  while (fs.good())
  {
    std::getline(fs, line);
    itemsRead = std::sscanf(line.c_str(), "%f %f %f %f %f %f %f %f %f",
      &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9);
    if (itemsRead != 9)
    {
      std::cerr << "ERROR: expeced a covariance entry at line: " << line << std::endl;
      assert(0);
    }
    cov.Assign(f1, f2, f3, f4, f5, f6, f7, f8, f9);
    covVct.push_back(cov);
    numCov++;
  }
  if (fs.bad() || fs.fail())
  {
    std::cerr << "ERROR: read pts from file failed; last line read: " << line << std::endl;
    assert(0);
  }
  fs.close();

  vctDynamicVector<vct3x3> covArray(numCov);
  for (unsigned int i = 0; i < numCov; i++)
  {
    covArray[i] = covVct[i];
  }

  return covArray;
}

void cov_write(vctDynamicVector<vct3x3> &cov, std::string &filename)
{
  std::ofstream out(filename.c_str());
  if (out) {
    unsigned int numCov = cov.size();
    for (unsigned int i = 0; i < numCov; i++)
    {
      out << cov[i].Row(0) << " " << cov[i].Row(1) << " " << cov[i].Row(2) << std::endl;
    }
  }
  else {
    out << "ERROR: cannot open covariance file for writing" << std::endl;
    assert(0);
  }
  out.close();
}

// Write cov to file
void WriteToFile_Cov(const vctDynamicVector<vct3x3> &cov,
  std::string &filePath)
{
  // Text file format:
  //
  //  c00 c01 c02 c10 c11 c12 c20 c21 c22
  //   ...
  //  c00 c01 c02 c10 c11 c12 c20 c21 c22
  //
  //  where vx's are indices into the pts array

  std::cout << "Saving cov to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filePath << std::endl;
    assert(0);
  }
  //fs << "COV " << cov.size() << "\n";
  unsigned int nCov = cov.size();
  for (unsigned int i = 0; i < nCov; i++)
  {
    fs << cov.at(i)(0, 0) << " " << cov.at(i)(0, 1) << " " << cov.at(i)(0, 2) << " "
      << cov.at(i)(1, 0) << " " << cov.at(i)(1, 1) << " " << cov.at(i)(1, 2) << " "
      << cov.at(i)(2, 0) << " " << cov.at(i)(2, 1) << " " << cov.at(i)(2, 2) << "\n";
  }
  fs.close();
}

// Write L to file
void WriteToFile_L(const vctDynamicVector<vctFixedSizeMatrix<double, 3, 2> > &L,
  std::string &filePath)
{
  // Text file format:
  //
  //  L1x L1y L1z L2x L2y L2z
  //   ...
  //  L1x L1y L1z L2x L2y L2z
  //
  //  where vx's are indices into the pts array

  std::cout << "Saving L to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cerr << "ERROR: failed to open file: " << filePath << std::endl;
    assert(0);
  }
  //fs << "L " << L.size() << "\n";
  unsigned int numL = L.size();
  for (unsigned int i = 0; i < numL; i++)
  {
    fs << L.at(i)(0, 0) << " " << L.at(i)(1, 0) << " " << L.at(i)(2, 0) << " "
      << L.at(i)(0, 1) << " " << L.at(i)(1, 1) << " " << L.at(i)(2, 1) << "\n";
  }
  fs.close();
}


double ExtractGaussianRVFromStream(std::ifstream &randnStream)
{
  double p;
  randnStream >> p;    // N(0,1) distributed RV
  if (!randnStream.good())
  {
    std::cout << "ERROR: could not read from Gaussian random number file!   Flags: "
      << randnStream.eofbit << " " << randnStream.failbit << " " << randnStream.badbit << std::endl;
    assert(0);
  }
  return p;
}


//// This function uses the Box-Muller method to generate
////  zero mean, unit variance Gaussian random variables
////  from a source of uniform random variables
//// Ref: http://www.springer.com/cda/content/document/cda_downloaddocument/9781447129929-c2.pdf?SGWID=0-0-45-1314437-p174313273
//double GenerateGaussianRV( cmnRandomSequence &rvGen )
//{
//  // 1: generate two uniform variables (U1 & U2) on (0,1]
//  //     Note: since log(U1) = inf, we must generate variables
//  //           on (0,1] rather than [0,1]
//  double U1 = rvGen.ExtractRandomDouble(1e-6, 1.0);
//  double U2 = rvGen.ExtractRandomDouble(1e-6, 1.0);
//
//  // 2: compute theta & p
//  double theta = 2*cmnPI*U2;
//  double p = sqrt(-2.0*log(U1));
//
//  // 3: compute gaussian RV
//  //     Note: Z1 = p*cos(theta) and Z2 = p*sin(theta)
//  //           are both Gaussian variables; to reduce
//  //           interdependency between generated variables,
//  //           we only use one of these.
//  double Z1 = p*cos(theta);
//  return Z1;
//}


void ComputeGroundTruthStats(const vctDynamicVector<vct3> &samples,
  const vctDynamicVector<vct3> &sampleNorms,
  const vctFrm3 &Fgt, const vctFrm3 &Freg,
  double &posErr_mean, double &posErr_SD,
  double &normErr_mean, double &normErr_SD)
{
  unsigned int nSamps = samples.size();

  // Note: samples represents validation samples, which are at the
  //       ground truth position.
  //       Fgt is the xfm that takes non-registered samples to the ground truth position
  //        (i.e. it is what Freg should be)
  //       Freg is the computed registration that takes non-registered samples to the
  //        estimated registered position.
  vctFrm3 dF(Freg*Fgt.Inverse());

  // compute ground truth distance statistics
  vctDynamicVector<double> posErr(samples.size());
  vctDynamicVector<double> normErr(sampleNorms.size());

  // position errors
  posErr_mean = 0.0;
  vct3 t1;
  double t2;
  for (unsigned int i = 0; i < nSamps; i++)
  { // compute mean distances
    t1 = samples[i] - dF*samples[i];
    posErr[i] = t1.Norm();
    posErr_mean += posErr[i];
  }
  posErr_mean /= nSamps;
  double posErr_Var = 0.0;
  for (unsigned int i = 0; i < nSamps; i++)
  { // compute standard deviation of distance error
    t2 = (posErr[i] - posErr_mean);
    posErr_Var += t2*t2;
  }
  posErr_Var /= nSamps;
  posErr_SD = sqrt(posErr_Var);

  // orientation errors
  normErr_mean = 0.0;
  double dotProd;
  for (unsigned int i = 0; i < nSamps; i++)
  { // compute mean angular offset
    dotProd = vctDotProduct(sampleNorms[i], dF.Rotation()*sampleNorms[i]);
    normErr[i] = acos(dotProd) * (180.0 / cmnPI);
    normErr_mean += normErr[i];
  }
  normErr_mean /= nSamps;

  double normErr_Var = 0.0;
  for (unsigned int i = 0; i < nSamps; i++)
  { // compute standard deviation of angular offset
    t2 = (normErr[i] - normErr_mean);
    normErr_Var += t2*t2;
  }
  normErr_Var /= nSamps;
  normErr_SD = sqrt(normErr_Var);
}


void GenerateRandomTransform(unsigned int randSeed, unsigned int &randSeqPos,
  double minOffsetPos, double maxOffsetPos,
  double minOffsetAng, double maxOffsetAng,
  vctFrm3 &F)
{
  //initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  // Generate random rotation
  vctRot3 R;
  GenerateRandomRotation(randSeed, randSeqPos, minOffsetAng, maxOffsetAng, R);

  // Generate random translation
  // get translation direction
  vct3 dir(0.0);
  while (dir.Norm() < 1e-6)
  {
    double pa[3];
    cisstRandomSeq.ExtractRandomDoubleArray(-1.0, 1.0, pa, 3);
    dir.Assign(pa[0], pa[1], pa[2]);
  }
  dir.NormalizedSelf();
  // get translation magnitude
  double mag = cisstRandomSeq.ExtractRandomDouble(minOffsetPos, maxOffsetPos);
  vct3 p = mag*dir;

  randSeqPos = cisstRandomSeq.GetSequencePosition();

  F.Assign(vctRot3(R), p);

  //// Generate random translation
  //double pa[3];
  //double sign;
  //cisstRandomSeq.ExtractRandomDoubleArray(minOffsetPos,maxOffsetPos, pa,3);
  //sign = cisstRandomSeq.ExtractRandomDouble(0.0,1.0);
  //if (sign < 0.5) sign = -1.0;
  //else sign = 1.0;
  //vct3 pi(sign*pa[0],sign*pa[1],sign*pa[2]);
}

void GenerateRandomRotation(unsigned int randSeed, unsigned int &randSeqPos,
  double minOffsetAng, double maxOffsetAng,
  vctRot3 &R)
{
  //initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  // Generate random rotation on the unit sphere with uniform rotation angle
  //  Note: using a fixed polar reference, such as the z-axis,
  //        results in rotation vectors that are more concentrated about
  //        the poles (i.e. about rotation angles near zero/180 degrees) than
  //        about the equator (i.e. rotation angles near 90 degrees), due to 
  //        increased surface area at the equator.
  //  Note: one way to fix this would be a to choose a random polar reference
  //
  //  use z-axis as starting point for random rotation axis
  //  set xy axis along random axis direction on the x-y plane
  //  generate uniform rotation angle 0-180 deg
  //    (covers both upper and lower hemisphere since only positive rotations
  //     will be generated for actual angle of rotation)
  //  apply rotation about the xy axis to the z-axis to get the random rotation axis
  //  generate unsigned uniform rotation angle around the random axis
  vct3 z(0.0, 0.0, 1.0);
  double xyDir = cisstRandomSeq.ExtractRandomDouble(0.0, 359.999)*cmnPI / 180.0;
  vct3 xyAx(cos(xyDir), sin(xyDir), 0.0);
  double xyAn = cisstRandomSeq.ExtractRandomDouble(0.0, 180.0)*cmnPI / 180.0;
  vctAxAnRot3 Rxy(xyAx, xyAn);
  vct3 rndAx = vctRot3(Rxy)*z;
  double rndAn = cisstRandomSeq.ExtractRandomDouble(minOffsetAng, maxOffsetAng)*cmnPI / 180.0;
  vctAxAnRot3 Ri(rndAx, rndAn);

  randSeqPos = cisstRandomSeq.GetSequencePosition();
  R.Assign(vctRot3(Ri));
}

void GenerateRandomL(unsigned int randSeed, unsigned int &randSeqPos,
  const vctDynamicVector<vct3> &uDir, vctDynamicVector<vct3x2> &L)
{
  //initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  unsigned int nL = uDir.size();
  L.SetSize(nL);

  for (unsigned int i = 0; i < nL; i++)
  {
    // Generate random major / minor axis for GIMLOP distribution where the
    //  major / minor axis are perpendicular to each other and perpendicular
    //  to the central direction
    // generate any vector not parallel to central direction
    vct3 v(0.0);
    double rv[3];
    while (v.Norm() < 1e-10)
    {
      cisstRandomSeq.ExtractRandomDoubleArray(0.1, 1.0, rv, 3);
      v.Assign(rv[0], rv[1], rv[2]);
      v.NormalizedSelf();
      // compute 1st vector perpendicular to central direction
      v = vctCrossProduct(v, uDir(i));
    }    
    L(i).Column(0) = v.Normalized();
    // comptue 2nd vector perpendicular to central direction and the 1st vector
    L(i).Column(1) = vctCrossProduct(uDir(i),L(i).Column(0)).Normalized();
  }
}

void Draw3DGaussianSample(std::ifstream &randnStream,
  const vct3x3 &M, vct3 &x)
{
  // Compute eigen decomposition of the covariance
  //  this method is stable under small perturbations of M across
  //  different machines
  vct3 eigenValues;
  vct3x3 eigenVectors;

  ComputeCovEigenDecomposition_SVD(M, eigenValues, eigenVectors);
  //ComputeCovEigenDecomposition_NonIter(M, eigenValues, eigenVectors);

  vct3x3 Ninv;
  Ninv.Column(0) = eigenVectors.Column(0)*sqrt(eigenValues[0]);
  Ninv.Column(1) = eigenVectors.Column(1)*sqrt(eigenValues[1]);
  Ninv.Column(2) = eigenVectors.Column(2)*sqrt(eigenValues[2]);

  // Generate an isotropic Gaussian noise and then transform it 
  //  to the noise model defined by covariance M
  //  f(x) ~ N(0,M) = x'inv(M)x = (Nx)'I(Nx)
  //  p = Nx ~ N(0,I)
  //  x = inv(N)*p
  //  M = VSV'  =>  Minv = V*S^-1*V'  =>  N = S^(-1/2)*V'  =>  Ninv = V*S^(1/2)
  double p1, p2, p3;
  p1 = ExtractGaussianRVFromStream(randnStream);
  p2 = ExtractGaussianRVFromStream(randnStream);
  p3 = ExtractGaussianRVFromStream(randnStream);
  vct3 p(p1, p2, p3);
  x = Ninv*p;

  // This is not stable for small differences in M
  //  i.e. M can have slightly different values when running code on
  //  different machines (like e-10 differences) and this creates large
  //  differences in the generated noise vector when computing Ninv
  //// compute covariance decomposition:
  ////    M = U*S*V' = U*S*U'
  ////    inv(N) = U*sqrt(S)
  ////    M = inv(N)*inv(N'),  inv(M) = N'*N
  //vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> A;
  //vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> U;
  //vctFixedSizeMatrix<double, 3, 3, VCT_COL_MAJOR> Vt;
  //vct3 S;
  //try
  //{
  //  A.Assign(M);
  //  nmrSVD(A, U, S, Vt);
  //}
  //catch (...)
  //{
  //  assert(0);
  //}
  //vct3x3 Ninv;
  //Ninv.Column(0) = U.Column(0)*sqrt(S[0]);
  //Ninv.Column(1) = U.Column(1)*sqrt(S[1]);
  //Ninv.Column(2) = U.Column(2)*sqrt(S[2]);
}

// Fisher Parameters
//  k     ~ concentration  (Note: here we use an alternative approximation k ~= 1.0/circSD^2)
//  mean  ~ central direction
void DrawFisherSample(std::ifstream &randnStream, double k, const vct3 &mean, vct3 &n)
{
  // Use the approximation that the tangential to the mean direction
  //  is approximately Gaussian distributed as: 
  //     sqrt(1/k)*Xperp ~ N(0,I2)  where  1/k ~= sigma2 (circular std dev squared)

  // Compute noisy version of z-axis as: [N1,N2,t]'
  double t, N1, N2;
  double circSD = 1.0 / sqrt(k);
  N1 = ExtractGaussianRVFromStream(randnStream) * circSD * (cmnPI / 180.0);
  N2 = ExtractGaussianRVFromStream(randnStream) * circSD * (cmnPI / 180.0);
  t = sqrt(1.0 - N1*N1 - N2*N2);
  vct3 zn(N1, N2, t);
  zn.NormalizedSelf();  // apply unit normalization to noisy z just to be safe
  // find rotation that rotates z to the mean direction
  vct3 z(0.0, 0.0, 1.0);
  vct3 xProd = vctCrossProduct(z, mean);
  if (xProd.Norm() < 1.0e-6)    // protect from divide-by-zero
  { // sample norm is parallel to z
    if (vctDotProduct(z, mean) > 0.0)
      n = zn;
    else
      n = -zn;
  }
  else
  {
    // the norm of the cross product is the same for angles of x deg & x+180 deg
    //  between two vectors => use dot product to determine the angle
    //   NOTE: the angle corresponding to the cross product axis is always > 0;
    //         acos of the dot product gives the correct form
    //   NOTE: the problem with using norm of cross product isn't that we aren't
    //         going the right direction, but rather that we don't rotate far enough
    //         if A & B are seperated by more than 90 degrees.  I.e. if angular
    //         seperation between A & B is 100 degrees, then asin(norm(AxB)) gives
    //         the same angle as if A & B are seperated by 80 degrees => we don't
    //         know if the actual angle is X or 90+X using the norm of cross product.
    double dProd = vctDotProduct(z, mean);
    vctAxAnRot3 Rs(xProd.Normalized(), acos(dProd));
    // apply rotation to noisy Z to obtain the noisy sample
    n = vctRot3(Rs)*zn;
  }

  //// Apply a Gaussian-distributed offset angle about a random rotation axis oriented
  ////  perpendicular to the orientation direction
  ////  Note: this method is more likely to generate angular offsets close to zero
  ////        than a proper Fisher or binomial Gaussian approximation. (To do this correctly,
  ////        the offset angle should be chi-square distributed.)
  //// Compute noisey version of +z axis:
  //// compute a random rotation axis on the xy plane to add noise to base (+Z) direction
  //double xyDir = cisstRandomSeq.ExtractRandomDouble(0.0, 359.999)*cmnPI/180.0;
  //vct3 xyAx(cos(xyDir),sin(xyDir),0.0);
  //double xyAn;      
  //xyAn = ExtractGaussianRVFromStream( randnStream ) * circStdDev * (cmnPI/180.0);
  //vctAxAnRot3 Rn(xyAx,xyAn);
  //vct3 zn = vctRot3(Rn)*z;     // noisy z             // *** why must convert to vctRot3 type?
  //// find rotation that rotates z to the sample norm
  //xProd = vctCrossProduct(z,sampleNorms.at(k));
  //if (xProd.Norm() < 1e-6)
  //{ // sample norm is parallel to z
  //  // protect from divide-by-zero
  //  if (vctDotProduct(z,sampleNorms.at(k)) > 0)
  //    noisySampleNorms.at(k) = zn;
  //  else
  //    noisySampleNorms.at(k) = -zn;
  //}
  //else
  //{
  //  // the norm of the cross product is the same for angles of x deg & x+180 deg
  //  //  between two vectors => use dot product to determine the angle
  //  //   NOTE: the angle corresponding to the cross product axis is always > 0;
  //  //         acos of the dot product gives the correct form
  //  //   NOTE: the problem with using norm of cross product isn't that we aren't
  //  //         going the right direction, but rather that we don't rotate far enough
  //  //         if A & B are seperated by more than 90 degrees.  I.e. if angular
  //  //         seperation between A & B is 100 degrees, then asin(norm(AxB)) gives
  //  //         the same angle as if A & B are seperated by 80 degrees => we don't
  //  //         know if the actual angle is X or 90+X using the norm of cross product.
  //  double dProd = vctDotProduct(z,sampleNorms.at(k));
  //  vctAxAnRot3 Rs(xProd.Normalized(), acos(dProd));
  //  // apply rotation to noisy Z to obtain the noisy sample
  //  noisySampleNorms.at(k) = vctRot3(Rs)*zn;
  //}
}

// GIMLOP Parameters
//  k  ~ concentration       (Note the approximation 1/k ~= circSD^2)
//  B  ~ ovalness parameter  (Note eccentricity e = 2*Beta/k) 
//  mean ~ central direction
//  L  ~ [l1,l2]
//        l1 ~ major axis
//        l2 ~ minor axis
void DrawGIMLOPSample(std::ifstream &randnStream,
  double k, double B, const vct3 &mean,
  const vctFixedSizeMatrix<double, 3, 2> &L, vct3 &n)
{
  // Use the approximation that the tangential to the mean direction
  //  is approximately Gaussian distributed as: 
  //     Z ~ N(0,[k*I3-2*A]^-1)  where  A = Beta*[l1*l1' - l2*l2'] 
  //                                    Z is the vector component tangential to the central direction
  //
  //  Since Z only has 2DOF, we can also compute Z using a 2D distribution:
  //   define Zp = [l1'Z, l2'Z]'
  //   then Zp is approximately 2D Gaussian distributed as: Zp ~ N(0,diag(k-2B,k+2B)^-1)
  //
  //  Alternatively, the 3D form may be expressed as:
  //   Z ~ N(0,M^-1)  where  M = (k-2B)*l1*l1' + (k+2B)*l2*l2' + k*u*u'
  //   which can be generated using 2 Gaussian distributed values
  //   for lengths along l1 & l2 since Z is known to lie in the plane of l1 & l2
  //
  //   Variance along l1 = 1/(k-2B)
  //   Variance along l2 = 1/(k+2B)
  //
  // Noisy vector X ~= u + Z  for concentrated points
  //   more accurately, X = Zwrap where Zwrap is the "wrapped" position of Z on the unit
  //   sphere starting from the central direction u

  if (k <= 2 * B)
  {
    std::cout << "ERROR: can only simulate for k > 2*B" << std::endl;
  }
  double N1, N2, SD1, SD2;
  N1 = ExtractGaussianRVFromStream(randnStream);
  N2 = ExtractGaussianRVFromStream(randnStream);
  SD1 = 1.0 / sqrt(k - 2.0*B);
  SD2 = 1.0 / sqrt(k + 2.0*B);
  // compute tangential vector Z
  vct3 Z = N1*SD1*L.Column(0) + N2*SD2*L.Column(1);
  // wrap Z onto unit sphere
  //  do this by rotating u until it aligns with Zwrap
  vct3 Raxis = vctCrossProduct(mean, Z);
  Raxis.NormalizedSelf();
  double Rangle = Z.Norm();
  if (Rangle > cmnPI)
  {
    std::cout << "WARNING: generated GIMLOP sample more than 180 degrees"
      << " offset from the central direction" << std::endl;
  }
  vctRot3 R(vctAxAnRot3(Raxis, Rangle));
  n = R*mean;   // Zwrap
}


void CreateMesh(cisstMesh &mesh,
  const std::string &meshLoadPath,
  std::string *SavePath_Mesh)
{
  // load mesh
  mesh.LoadPLY(meshLoadPath);
  if (mesh.NumVertices() == 0)
  {
    printf("ERROR: Read mesh resulted in 0 triangles\n");
    assert(0);
  }

  // save mesh
  if (SavePath_Mesh)
  {
    // TODO: update
    //if (mesh.SaveMeshFile((*SavePath_Mesh).append(".mesh")) < 0)
    {
      std::cout << "ERROR: Save mesh failed" << std::endl;
      assert(0);
    }
  }
}


void GenerateSamples(cisstMesh &mesh,
  unsigned int randSeed, unsigned int &randSeqPos,
  unsigned int nSamps,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<unsigned int>  &sampleDatums,
  std::string *SavePath_Samples)
{
  //initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  samples.SetSize(nSamps);
  sampleNorms.SetSize(nSamps);
  sampleDatums.SetSize(nSamps);

  // generating samples
  int Fx, Nf;
  double mu0, mu1, mu2, mu;
  double lam0, lam1, lam2;
  // initialize samples
  Nf = mesh.NumTriangles();
  for (unsigned int s = 0; s < nSamps; s++) {
    // sampling the object surface and generating second pointcloud from it
    Fx = cisstRandomSeq.ExtractRandomInt(0, Nf - 1);
    mu0 = cisstRandomSeq.ExtractRandomDouble(0.0, 1.0);
    mu1 = cisstRandomSeq.ExtractRandomDouble(0.0, 1.0);
    mu2 = cisstRandomSeq.ExtractRandomDouble(0.0, 1.0);
    mu = mu0 + mu1 + mu2;
    lam0 = mu0 / mu;
    lam1 = mu1 / mu;
    lam2 = mu2 / mu;
    vct3 v0, v1, v2;
    mesh.FaceCoords(Fx, v0, v1, v2);
    samples.at(s) = lam0*v0 + lam1*v1 + lam2*v2;
    sampleNorms.at(s) = mesh.faceNormals(Fx);
    sampleDatums.at(s) = Fx;
  }
  randSeqPos = cisstRandomSeq.GetSequencePosition();

  // save samples
  if (SavePath_Samples)
  {
    if (cisstPointCloud::WritePointCloudToFile(*SavePath_Samples, samples, sampleNorms) < 0)
    {
      std::cout << "ERROR: Samples save failed" << std::endl;
      assert(0);
    }
  }
}

void GenerateNoisySamples_Gaussian(
  std::ifstream &randnStream,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3x3> &sampleCov,
  vctDynamicVector<vct3>   &noisySamples,
  std::string *SavePath_NoisySamples,
  std::string *SavePath_Cov)
{
  unsigned int numSamps = samples.size();
  noisySamples.SetSize(numSamps);
  for (unsigned int i = 0; i < numSamps; i++)
  {
    //=== Generate position noise ===//

    // generate Guassian noise
    vct3 noise;
    Draw3DGaussianSample(randnStream, sampleCov(i), noise);

    // apply Gaussian noise to sample
    noisySamples(i) = samples(i) + noise;
  }

  // save noisy samples
  if (SavePath_NoisySamples)
  {
    std::string noisySampsSavePath = *SavePath_NoisySamples;
    cisstPointCloud::WritePointCloudToFile(noisySampsSavePath, noisySamples);
  }
  // save sample covariances
  if (SavePath_Cov)
  {
    WriteToFile_Cov(sampleCov, *SavePath_Cov);
  }
}

void GenerateOutlierSamples_SurfaceOffset(
  unsigned int randSeed, unsigned int &randSeqPos,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &outlierSamples,
  double minPosOffsetOutlier, double maxPosOffsetOutlier,
  std::string *SavePath_OutlierSamples)
{
  // add outliers to samples such that outliers are offset outward from
  //  the shape (offset along the point normal direction)

  // initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  unsigned int numOutliers = samples.size();
  outlierSamples.SetSize(numOutliers);
  for (unsigned int k = 0; k < numOutliers; k++)
  {
    //=== Generate position outlier ===//

    // generate a random outlier
    double offset = cisstRandomSeq.ExtractRandomDouble(minPosOffsetOutlier, maxPosOffsetOutlier);
    outlierSamples(k) = samples(k) + offset*sampleNorms(k);
  }

  randSeqPos = cisstRandomSeq.GetSequencePosition();

  // save noisy samples
  if (SavePath_OutlierSamples)
  {
    std::string savePath = *SavePath_OutlierSamples;
    cisstPointCloud::WritePointCloudToFile(savePath, outlierSamples);
  }
}

// Generate noisy samples having the specified standard deviation of
//  noise in directions parallel and perpendicular to the surface
void ComputeCovariances_Random( 
  unsigned int randSeed, unsigned int &randSeqPos,
  const vct3 &StdDev,
  vctDynamicVector<vct3x3> &cov,
  unsigned int nCov,
  std::string *SavePath_Cov)
{
  // initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  vctRot3 R;
  vct3x3 M0;

  // define diagonal matrix of eigenvalues
  M0.SetAll(0.0);
  M0.Element(0, 0) = StdDev[0] * StdDev[0];
  M0.Element(1, 1) = StdDev[1] * StdDev[1];
  M0.Element(2, 2) = StdDev[2] * StdDev[2];

  // Compute covariances
  cov.SetSize(nCov);
  for (unsigned int i = 0; i < nCov; i++)
  {
    // generate random rotation axis
    vct3 axis;    
    cisstRandomSeq.ExtractRandomDoubleArray(-1.0, 1.0, axis.Pointer(0), 3);    
    if (axis.Norm() < 1e-12)  // protect from division by zero
    {
      axis = (1.0, 0.0, 0.0);
    }
    axis.NormalizedSelf();
    // generate random rotation angle
    double angle = cisstRandomSeq.ExtractRandomDouble(-cmnPI, cmnPI);
    // form rotation matrix
    vctRot3 R = vctRot3(vctAxAnRot3(axis, angle));

    // form covariance 
    cov(i) = R*M0*R.Transpose();
  }

  // save model covariances
  if (SavePath_Cov)
  {
    WriteToFile_Cov(cov, *SavePath_Cov);
  }

  randSeqPos = cisstRandomSeq.GetSequencePosition();
}

// Generate noisy samples having the specified standard deviation of
//  noise in directions parallel and perpendicular to the surface
void ComputeCovariances_SurfaceModel(
  double StdDevInPlane, double StdDevPerpPlane,
  vctDynamicVector<vct3>   &samples,
  vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3x3> &modelCov,
  //vctDynamicVector<vct3x3> &modelInvCov,
  std::string *SavePath_Cov)
{
  unsigned int nSamps = samples.size();

  vctRot3 R;
  vct3x3 M, M0;
  //vct3x3 invM0, Ninv;
  vct3 x(1.0, 0.0, 0.0);
  vct3 y(0.0, 1.0, 0.0);
  vct3 z(0.0, 0.0, 1.0);

  // Covariance Model
  //  define model such that noise perpendicular to the surface
  //  is different than noise in-plane with the surface
  //  assume a z-axis surface orientation
  M0.SetAll(0.0);
  M0.Element(0, 0) = StdDevInPlane * StdDevInPlane;
  M0.Element(1, 1) = StdDevInPlane * StdDevInPlane;
  M0.Element(2, 2) = StdDevPerpPlane * StdDevPerpPlane;
  //invM0.SetAll(0.0);
  //invM0.Element(0, 0) = 1.0 / M0.Element(0, 0);
  //invM0.Element(1, 1) = 1.0 / M0.Element(1, 1);
  //invM0.Element(2, 2) = 1.0 / M0.Element(2, 2);

  // Compute model covariances
  modelCov.SetSize(nSamps);
  //modelInvCov.SetSize(nSamps);
  for (unsigned int i = 0; i < nSamps; i++)
  {
    //  find a rotation to rotate the sample normal to the z-axis
    R = XProdRotation(sampleNorms(i), z);
    // compute model covariance M of this sample
    //   Note: rotate to align normal with z-axis, apply covariance model, rotate back
    M = R.Transpose()*M0*R;

    modelCov(i) = M;
    //modelInvCov(i) = R.Transpose()*invM0*R;
  }

  // save model covariances
  if (SavePath_Cov)
  {
    WriteToFile_Cov(modelCov, *SavePath_Cov);
  }
}

// Generate noisy samples with noise distributed according the
//  given covariance
void GenerateSampleErrors_Covariance(
  std::ifstream &randnStream,
  const vctDynamicVector<vct3x3> &sampleCov,
  const vctDynamicVector<vct3>   &samples,
  const vctDynamicVector<vct3>   &sampleNorms,
  vctDynamicVector<vct3>   &noisySamples,
  std::string *SavePath_NoisySamples)
{
  unsigned int nSamps = samples.size();  

  // Compute noise
  noisySamples.SetSize(nSamps);
  GenerateNoisySamples_Gaussian(randnStream, samples, sampleCov, noisySamples);

  // save noisy samples
  if (SavePath_NoisySamples)
  {
    std::string savePath = *SavePath_NoisySamples;
    cisstPointCloud::WritePointCloudToFile(savePath, noisySamples);
  }
}

// Generate noisy samples having the specified standard deviation of
//  noise in directions parallel and perpendicular to the surface
//  with the specified percentage of samples being outliers
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
  std::string *SavePath_NoisySamples,
  std::string *SavePath_Cov)
{
  unsigned int nSamps = samples.size();
  unsigned int nOutliers = (unsigned int)(percentOutliers*(double)nSamps);
  unsigned int nNoisy = nSamps - nOutliers;

  // Compute Covariance Model
  //  define model such that noise perpendicular to the surface
  //  is different than noise in-plane with the surface
  ComputeCovariances_SurfaceModel(
    StdDevInPlane, StdDevPerpPlane,
    samples, sampleNorms,
    sampleCov,
    SavePath_Cov);

  // Compute noise
  vctDynamicVector<vct3> samplesNoisy(nNoisy);
  GenerateNoisySamples_Gaussian(randnStream, samples, sampleCov, samplesNoisy);

  // Compute outliers
  vctDynamicVectorRef<vct3> samplesRef;
  vctDynamicVectorRef<vct3> sampleNormsRef;
  samplesRef.SetRef(samples, nNoisy, nOutliers);
  sampleNormsRef.SetRef(sampleNorms, nNoisy, nOutliers);
  vctDynamicVector<vct3> samplesOutliers(nOutliers);
  GenerateOutlierSamples_SurfaceOffset(
    randSeed, randSeqPos,
    samplesRef, sampleNormsRef,
    samplesOutliers,
    minPosOffsetOutlier, maxPosOffsetOutlier);

  // compose output samples
  noisySamples.SetSize(nSamps);
  for (unsigned int i = 0; i < nNoisy; i++)
  {
    noisySamples(i) = samplesNoisy(i);
  }
  for (unsigned int i = nNoisy; i < nSamps; i++)
  {
    noisySamples(i) = samplesOutliers(i - nNoisy);
  }

  // save noisy samples
  if (SavePath_NoisySamples)
  {
    std::string savePath = *SavePath_NoisySamples;
    cisstPointCloud::WritePointCloudToFile(savePath,noisySamples);
  }
}


// Generate sample noise having a random position covariance with the
//  specified eigenvalues and random orientation with the specified
//  circular standard deviation and eccentricity.
void GenerateSampleRandomNoise(unsigned int randSeed, unsigned int &randSeqPos,
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
  std::string *SavePath_NoisySamples,
  std::string *SavePath_Cov,
  std::string *savePath_L)
{
  // initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  unsigned int nSamps = samples.size();
  unsigned int nOutliers = (unsigned int)(percentOutliers*(double)nSamps);
  unsigned int nNoisy = nSamps - nOutliers;

  noisySamples.SetSize(nSamps);
  noisySampleNorms.SetSize(nSamps);
  sampleCov.SetSize(nSamps);
  sampleInvCov.SetSize(nSamps);
  noiseL.SetSize(nSamps);

  vct3x3 M, M0, invM0, Ninv;
  vct3 n;
  vctRot3 R;
  vct3 x(1.0, 0.0, 0.0);
  vct3 y(0.0, 1.0, 0.0);
  vct3 z(0.0, 0.0, 1.0);

  // set eigenvalues of noise covariance
  //  assume a z-axis normal orientation
  M0.SetAll(0.0);
  M0.Element(0, 0) = covEigenvalues[0];
  M0.Element(1, 1) = covEigenvalues[1];
  M0.Element(2, 2) = covEigenvalues[2];
  invM0.SetAll(0.0);
  invM0.Element(0, 0) = 1.0 / M0.Element(0, 0);
  invM0.Element(1, 1) = 1.0 / M0.Element(1, 1);
  invM0.Element(2, 2) = 1.0 / M0.Element(2, 2);

  // add random noise to samples such that noise perpendicular to the
  //  sample triangle is different than noise in-plane with triangle
  for (unsigned int i = 0; i < nNoisy; i++)
  {
    //=== Generate position noise ===//

    // Generate random noise covariance for this sample
    vctRot3 Rcov;
    GenerateRandomRotation(randSeed, randSeqPos, 0.0, 180.0, Rcov);
    M = Rcov.Transpose()*M0*Rcov;

    // generate Guassian noise
    vct3 p;
    Draw3DGaussianSample(randnStream, M, p);

    // apply Gaussian noise to the sample
    noisySamples(i) = samples(i) + p;
    sampleCov(i) = M;
    sampleInvCov(i) = Rcov.Transpose()*invM0*Rcov;

    //=== Generate orientation noise ===//

    // Fisher Noise Model
    // Uses the approximation that 1/k ~= circSD^2
    //double k = 1.0/(circStdDev*circStdDev);
    //DrawFisherSample( randnStream, k, sampleNorms(i), noisySampleNorms(i) );

    // Generate a major/minor axis for this sample
    vctFixedSizeMatrix<double, 3, 2> L;
    // generate a random major/minor axis in the x-y plane
    double xyAngle = cisstRandomSeq.ExtractRandomDouble(0.0, cmnPI);
    vctRot3 Rz(vctAxAnRot3(z, xyAngle));
    vct3 majorXY = Rz*x;
    vct3 minorXY = Rz*y;
    // rotate major/minor axis to be perpendicular to the sample orientation
    //  apply rotation that takes z-axis to the sample norm
    R = XProdRotation(sampleNorms(i), z);  // rotation takes sample norm to z-axis
    L.Column(0) = R.Transpose()*majorXY;
    L.Column(1) = R.Transpose()*minorXY;

    // GIMLOP Noise Model
    // Use the approximation that k ~= 1/circSD^2
    double k = 1.0 / (circStdDev*circStdDev);
    double B = circEccentricity*k / 2.0;
    DrawGIMLOPSample(randnStream, k, B, sampleNorms(i), L, noisySampleNorms(i));

    // TODO: change to 3x3 L allowing for a different central direction
    //       than the observed sample
    // Reorient "reported" L to be perpendicular to the noisySample rather than
    //  the non-noisy sample (do this because the noisy sample will be used as
    //  the central direction in the GIMLOP noise model to simulate an observed
    //  data point with unknown ground truth orientation)
    // Compute rotation from sample norm to noisy sample norm
    vctRot3 Rl = XProdRotation(sampleNorms(i), noisySampleNorms(i));
    noiseL(i) = Rl*L;
  }

  // TODO: how to define noise model for outliers (i.e. what to assign
  //       as covariance, k, e, and L values?)
  // TODO: need seperate k, e, and L for each sample (not just seperate L) if
  //       using outliers?
  // add outliers to samples such that outliers are offset outward from
  //  the shape (offset along the point normal direction)
  vctRot3 Routlier;
  for (unsigned int k = nNoisy; k < nSamps; k++)
  {
    //=== Generate position outlier ===//

    double offset = cisstRandomSeq.ExtractRandomDouble(minPosOffsetOutlier, maxPosOffsetOutlier);
    noisySamples.at(k) = samples.at(k) + offset*sampleNorms.at(k);
    sampleCov(k) = vct3x3::Eye();
    sampleInvCov(k) = vct3x3::Eye();


    //=== Generate orientation outlier ===//

    // Generate random rotation on the unit sphere with uniform rotation angle
    //  within the specified range for the outliers
    GenerateRandomRotation(randSeed, randSeqPos, minAngOffsetOutlier, maxAngOffsetOutlier, Routlier);
    noisySampleNorms.at(k) = Routlier*sampleNorms.at(k);
    // define L
    // find any two axis perpendicular to the noisy sample and to each other
    vct3 xProd = vctCrossProduct(z, noisySampleNorms.at(k));
    vct3 tmp1, tmp2;
    if (xProd.Norm() <= 1.0e-8)  // protect from divide by zero
    { // noisy samples points down z axis
      tmp1 = x;
      tmp2 = y;
    }
    else
    {
      tmp1 = xProd.Normalized();
      tmp2 = vctCrossProduct(tmp1, noisySampleNorms.at(k));
      tmp2.NormalizedSelf();
    }
    noiseL(k).Column(0) = tmp1;
    noiseL(k).Column(1) = tmp2;
  }
  randSeqPos = cisstRandomSeq.GetSequencePosition();

  // save noisy samples
  if (SavePath_NoisySamples)
  {
    std::string noisySampsSavePath = *SavePath_NoisySamples;
    if (cisstPointCloud::WritePointCloudToFile(noisySampsSavePath, noisySamples, noisySampleNorms) < 0)
    {
      std::cout << "ERROR: Samples save failed" << std::endl;
      assert(0);
    }
  }
  // save cov
  if (SavePath_Cov)
  {
    WriteToFile_Cov(sampleCov, *SavePath_Cov);
  }
  // save L
  if (savePath_L)
  {
    WriteToFile_L(noiseL, *savePath_L);
  }
}


// Generate sample noise having the specified standard deviation of
//  noise in directions parallel and perpendicular to the triangle
//  plane from which the sample was drawn.
void GenerateSampleSurfaceNoise(unsigned int randSeed, unsigned int &randSeqPos,
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
  std::string *SavePath_NoisySamples,
  std::string *SavePath_Cov,
  std::string *savePath_L)
{
  // initialize random numbers
  cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  cisstRandomSeq.SetSeed(randSeed);
  cisstRandomSeq.SetSequencePosition(randSeqPos);

  unsigned int nSamps = samples.size();
  unsigned int nOutliers = (unsigned int)(percentOutliers*(double)nSamps);
  unsigned int nNoisy = nSamps - nOutliers;

  noisySamples.SetSize(nSamps);
  noisySampleNorms.SetSize(nSamps);
  sampleCov.SetSize(nSamps);
  sampleInvCov.SetSize(nSamps);
  noiseL.SetSize(nSamps);

  vct3x3 M, M0, invM0, Ninv;
  vct3 n;
  vctRot3 R;
  vct3 x(1.0, 0.0, 0.0);
  vct3 y(0.0, 1.0, 0.0);
  vct3 z(0.0, 0.0, 1.0);

  // set eigenvalues of noise covariance
  //  assume a z-axis normal orientation
  M0.SetAll(0.0);
  M0.Element(0, 0) = StdDevInPlane * StdDevInPlane;
  M0.Element(1, 1) = StdDevInPlane * StdDevInPlane;
  M0.Element(2, 2) = StdDevPerpPlane * StdDevPerpPlane;
  invM0.SetAll(0.0);
  invM0.Element(0, 0) = 1.0 / M0.Element(0, 0);
  invM0.Element(1, 1) = 1.0 / M0.Element(1, 1);
  invM0.Element(2, 2) = 1.0 / M0.Element(2, 2);

  // add random noise to samples such that noise perpendicular to the
  //  sample triangle is different than noise in-plane with triangle
  for (unsigned int i = 0; i < nNoisy; i++)
  {
    //=== Generate position noise ===//

    // Define the noise covariance for this sample
    //  find rotation to rotate the sample normal to the z-axis
    R = XProdRotation(sampleNorms(i), z);
    // compute noise covariance M of this sample
    //   Note: rotate to align normal with z-axis, apply noise covariance, rotate back
    M = R.Transpose()*M0*R;

    // generate Guassian noise
    vct3 p;
    Draw3DGaussianSample(randnStream, M, p);

    // apply Gaussian noise to the sample
    noisySamples(i) = samples(i) + p;
    sampleCov(i) = M;
    sampleInvCov(i) = R.Transpose()*invM0*R;

    if (circStdDev > 0.0)
    {
      //=== Generate orientation noise ===//

      // Fisher Noise Model
      // Uses the approximation that 1/k ~= circSD^2
      //double k = 1.0/(circStdDev*circStdDev);
      //DrawFisherSample( randnStream, k, sampleNorms(i), noisySampleNorms(i) );

      // Generate a major/minor axis for this sample
      vctFixedSizeMatrix<double, 3, 2> L;
      // generate a random major/minor axis in the x-y plane
      double xyAngle = cisstRandomSeq.ExtractRandomDouble(0.0, cmnPI);
      vctRot3 Rz(vctAxAnRot3(z, xyAngle));
      vct3 majorXY = Rz*x;
      vct3 minorXY = Rz*y;
      // rotate major/minor axis to be perpendicular to the sample orientation
      //  apply rotation that takes z-axis to the sample norm
      L.Column(0) = R.Transpose()*majorXY;
      L.Column(1) = R.Transpose()*minorXY;

      // GIMLOP Noise Model
      // Use the approximation that k ~= 1/circSD^2
      double k = 1.0 / (circStdDev*circStdDev);
      double B = circEccentricity*k / 2.0;
      DrawGIMLOPSample(randnStream, k, B, sampleNorms(i), L, noisySampleNorms(i));

      // Reorient "reported" L to be perpendicular to the noisySample rather than
      //  the non-noisy sample (do this because the noisy sample will be used as
      //  the central direction in the GIMLOP noise model to simulate an observed
      //  data point with unknown ground truth orientation)
      // Compute rotation from sample norm to noisy sample norm
      vctRot3 Rl = XProdRotation(sampleNorms(i), noisySampleNorms(i));
      noiseL(i) = Rl*L;
    }
    else
    {
      noisySampleNorms(i) = sampleNorms(i);
      noiseL(i) = vct3x2(0.0);
    }
  }

  // add outliers to samples such that outliers are offset outward from
  //  the shape (offset along the point normal direction)
  for (unsigned int k = nNoisy; k < nSamps; k++)
  {
    //=== Generate position outlier ===//

    // generate a random outlier
    double offset = cisstRandomSeq.ExtractRandomDouble(minPosOffsetOutlier, maxPosOffsetOutlier);
    noisySamples.at(k) = samples.at(k) + offset*sampleNorms.at(k);

    // Define the noise covariance for this sample
    //  find rotation to rotate the sample normal to the z-axis
    R = XProdRotation(sampleNorms.at(k), z);
    // compute noise covariance M of this sample
    //   Note: rotate to align normal with z-axis, apply noise covariance, rotate back
    M = R.Transpose()*M0*R;
    sampleCov.at(k) = M;
    sampleInvCov.at(k) = R.Transpose()*invM0*R;

    if (maxAngOffsetOutlier > 0.0)
    {
      //=== Generate orientation outlier ===//

      // Generate random rotation on the unit sphere with uniform rotation angle
      //  use z-axis as starting point for random rotation axis
      //  set rotation axis along random direction on the x-y plane
      //  generate uniform rotation angle 0-180 deg
      //  apply rotation about the rotation axis to the z-axis to get the random rotation axis
      //  generate a uniformly distributed rotation angle for the random rotation axis
      vct3 z(0.0, 0.0, 1.0);
      double xyDir = cisstRandomSeq.ExtractRandomDouble(0.0, 359.999)*cmnPI / 180.0;
      vct3 xyAx(cos(xyDir), sin(xyDir), 0.0);
      double xyAn = cisstRandomSeq.ExtractRandomDouble(0.0, 180.0)*cmnPI / 180.0;
      vctAxAnRot3 Rxy(xyAx, xyAn);
      vct3 rndAx = vctRot3(Rxy)*z;
      double rndAn = cisstRandomSeq.ExtractRandomDouble(minAngOffsetOutlier, maxAngOffsetOutlier)*(cmnPI / 180.0);
      vctAxAnRot3 Rrod(rndAx, rndAn);
      noisySampleNorms.at(k) = vctRot3(Rrod)*sampleNorms.at(k);
      // set L to any set of axis perpendicular to the sample
      vct3 a = vctCrossProduct(z, sampleNorms.at(k));
      if (a.Norm() > 0.001)
      {
        noiseL.at(k).Column(0) = a.NormalizedSelf();
        noiseL.at(k).Column(1) = vctCrossProduct(a, sampleNorms.at(k)).Normalized();
      }
      else
      {
        noiseL.at(k).Column(0) = vct3(1.0, 0.0, 0.0);
        noiseL.at(k).Column(1) = vct3(0.0, 1.0, 0.0);
      }
    }
    else
    {
      noisySampleNorms.at(k) = sampleNorms.at(k);
      noiseL.at(k) = vct3x2(0.0);
    }
  }
  randSeqPos = cisstRandomSeq.GetSequencePosition();

  // save noisy samples
  if (SavePath_NoisySamples)
  {
    std::string noisySampsSavePath = *SavePath_NoisySamples;
    if (cisstPointCloud::WritePointCloudToFile(noisySampsSavePath, noisySamples, noisySampleNorms) < 0)
    {
      std::cout << "ERROR: Samples save failed" << std::endl;
      assert(0);
    }
  }
  // save cov
  if (SavePath_Cov)
  {
    WriteToFile_Cov(sampleCov, *SavePath_Cov);
  }
  // save L
  if (savePath_L)
  {
    WriteToFile_L(noiseL, *savePath_L);
  }
}

#if 0   // TODO: code needs to be updated for changes to library

// Generate noise having different variance perpendicular to its mesh
//  triangle than noise in-plane
void GenerateNoisyMesh(cisstMesh &mesh,
  cisstMesh &noisyMesh,
  std::ifstream &randnStream,
  //unsigned int randSeed, unsigned int &randSeqPos
  double noiseStdDev,
  std::string *SavePath_NoisyMesh)
{
  ////initialize random numbers
  //cmnRandomSequence &cisstRandomSeq = cmnRandomSequence::GetInstance();
  //cisstRandomSeq.SetSeed( randSeed );
  //cisstRandomSeq.SetSequencePosition( 0 );

  unsigned int numTriangles = mesh.NumTriangles();
  unsigned int numVertices = mesh.NumVertices();

  noisyMesh.faces.SetSize(numTriangles);
  noisyMesh.vertices.SetSize(numVertices);
  noisyMesh.TriangleCov.SetSize(numTriangles);
  noisyMesh.TriangleCovEig.SetSize(numTriangles);

  // add isotropic random noise to the mesh vertices
  for (unsigned int k = 0; k < numVertices; k++)
  {
    //=== generate noisy vertices ===//

    // Generate an isotropic Gaussian noise scaled by the noise standard deviation
    double p1, p2, p3;
    //p1 = GenerateGaussianRV(cisstRandomSeq);
    //p2 = GenerateGaussianRV(cisstRandomSeq);
    //p3 = GenerateGaussianRV(cisstRandomSeq);
    p1 = ExtractGaussianRVFromStream(randnStream);
    p2 = ExtractGaussianRVFromStream(randnStream);
    p3 = ExtractGaussianRVFromStream(randnStream);
    vct3 p(p1, p2, p3);
    p.Multiply(noiseStdDev);

    // Assign noisy vertex to noisy mesh
    noisyMesh.vertices.Element(k) = mesh.vertices.Element(k) + p;
  }

  // construct triangles for noisy mesh
  cisstTriangle T(&noisyMesh);
  for (unsigned int i = 0; i < numTriangles; i++)
  {
    // get triangle vertex indices from old mesh
    int v0, v1, v2;
    mesh.TriangleVertexIndexes(i, v0, v1, v2);
    // copy indices to new triangle and compute normal vector
    T.SetVertexIndexes(v0, v1, v2);
    T.UpdatePlaneFromVertices();
    // assign triangle to noisy mesh
    noisyMesh.Triangles.Element(i) = T;
  }

  double var = noiseStdDev * noiseStdDev;
  vct3x3 M(0.0);
  M.Element(0, 0) = var;
  M.Element(1, 1) = var;
  M.Element(2, 2) = var;
  noisyMesh.TriangleCov.SetAll(M);
  noisyMesh.TriangleCovEig.SetAll(vct3(var));

  //randSeqPos = cisstRandomSeq.GetSequencePosition();

  // save noisy target mesh
  if (SavePath_NoisyMesh)
  {
    if (noisyMesh.SaveMeshFile((*SavePath_NoisyMesh).append(".mesh")) < 0)
    {
      std::cout << "ERROR: Save mesh failed" << std::endl;
      assert(0);
    }
  }
}


// Defines the measurement noise for each triangle in a mesh
//  (not the measurement noise of samples, but of the mesh itself)
void SetMeshTriangleCovariances(cisstMesh &mesh,
  double stdDevPerpPlane, double stdDevInPlane)
{
  unsigned int numTriangles = mesh.NumTriangles();
  vct3x3 M;
  vct3 n;
  vctRot3 R;
  vct3 z(0.0, 0.0, 1.0);
  double EigMax = stdDevPerpPlane > stdDevInPlane ? stdDevPerpPlane : stdDevInPlane;
  EigMax = EigMax*EigMax;

  // ensure normals have been computed for this mesh
  if (mesh.Triangles[0].norm.NormSquare() < 0.1)
  {
    std::cout << "ERROR: it appears normals were not computed for this mesh!" << std::endl;
    assert(0);
  }

  // initialize size of triangle covariances
  mesh.TriangleCov.SetSize(numTriangles);
  mesh.TriangleCovEig.SetSize(numTriangles);
  mesh.TriangleCovEig.SetAll(vct3(EigMax));

  // set eigenvalues of noise covariance
  M.SetAll(0.0);
  M.Element(0, 0) = stdDevInPlane*stdDevInPlane;
  M.Element(1, 1) = stdDevInPlane*stdDevInPlane;
  M.Element(2, 2) = stdDevPerpPlane*stdDevPerpPlane;     // z-axis chosen for perpendicular noise component

  for (unsigned int i = 0; i < numTriangles; i++)
  {
    // get normal to triangle plane
    n = mesh.Triangles[i].norm;

    // find rotation to align triangle normal with the z-axis
    vct3 xProd = vctCrossProduct(n, z);
    if (xProd.Norm() <= 1e-6)  // protect from divide by zero
    { // norm is already oriented with z-axis
      R = vctRot3::Identity();
    }
    else
    {
      // the norm of the cross product is the same for angles of x deg & x+180 deg
      //  between two vectors => use dot product to determine the angle
      //   NOTE: the angle corresponding to the cross product axis is always > 0;
      //         acos of the dot product gives the correct form
      //   NOTE: the problem with using norm of cross product isn't that we aren't
      //         going the right direction, but rather that we don't rotate far enough
      //         if A & B are seperated by more than 90 degrees.  I.e. if angular
      //         seperation between A & B is 100 degrees, then asin(norm(AxB)) gives
      //         the same angle as if A & B are seperated by 80 degrees => we don't
      //         know if the actual angle is X or 90+X using the norm of cross product.
      vct3 ax = xProd.Normalized();
      double an = acos(vctDotProduct(n, z));
      //double an = asin(t.Norm());
      vctAxAnRot3 R_AxAn(ax, an);
      R = vctRot3(R_AxAn);
    }

    // set noise covariance of this triangle
    //  rotate to align normal with perpendicular, apply noise covariance, rotate back
    mesh.TriangleCov[i] = R.Transpose()*M*R;
  }
}
#endif

// TODO: enable
#if 0

void Callback_TrackRegPath_Utility( cisstICP::CallbackArg &arg, void *userData )
{
  // cast to norm callback arg type
  cisstICP::CallbackArgNormals *argp;
  argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  if (argp==0)
  {
    std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
    assert(argp);   // terminate program
  }

  // Save to file:
  //  - error function (-loglik)
  //  - incremental transform
  //  - IMLP params
  //  - residual match errors
  // output format:
  //  err r00 r01 r02 r10 r11 r12 r20 r21 r22 tx ty tz normWeight posWeight avgNormError avgPosError 
  std::ofstream *fs = (std::ofstream *)(userData);
  (*fs) << argp->E << " " << argp->dF.Rotation().Row(0) << " " << argp->dF.Rotation().Row(1) << " " 
    << argp->dF.Rotation().Row(2) << " " << argp->dF.Translation() << " " 
    << argp->normWeight << " " << argp->posWeight << " "
    << argp->MatchPosErrorAvg << " " << argp->MatchPosErrorSD << " "
    << argp->MatchNormErrorAvg << " " << argp->MatchNormErrorSD << std::endl;
}

void Callback_SaveIterationsToFile_Utility( cisstICP::CallbackArg &arg, void *userData )
{
  // cast to norm callback arg type
  cisstICP::CallbackArgNormals *argp;
  argp = dynamic_cast<cisstICP::CallbackArgNormals*>(&arg);
  if (argp==0)
  {
    std::cerr << "ERROR: cannot cast callback argument to cisstICP type" << std::endl;
    assert(argp);   // terminate program
  }

  std::ofstream *fs = (std::ofstream *)(userData);
  vctRodRot3 dR(argp->dF.Rotation());
  std::stringstream ss;
  ss << cmnPrintf("iter=%u E=%.1f tolE=%.4f  (dAng/dPos)=%.2f/%.2f  t=%.3f  nW/pW=%.4f/%.4f  (circSD/posSD)=%.2f/%.2f  NNodes=%u/%u/%u NOut=%u/%u  PosErr=%.2f/%.2f  NormErr=%.2f/%.2f")
    //  dTheta=%.2f/%.2f/%.2f  
    << argp->iter 
    << argp->E
    << argp->tolE
    << dR.Norm()*180.0/cmnPI << arg.dF.Translation().Norm()
    << argp->time
    << argp->normWeight << argp->posWeight
    << argp->circSD*180.0/cmnPI << sqrt(argp->posVar)
    //<< argp->dThetaMin*180.0/cmnPI << argp->dThetaMax*180.0/cmnPI << argp->dThetaAvg*180.0/cmnPI
    << argp->maxNodesSearched << argp->avgNodesSearched << argp->minNodesSearched
    << argp->nOutliersPos << argp->nOutliersNorm
    << argp->MatchPosErrorAvg << argp->MatchPosErrorSD
    << argp->MatchNormErrorAvg << argp->MatchNormErrorSD;
  (*fs) << ss.str();
  if (arg.isAccelStep)
  {
    (*fs) << "\t-- accel step --";
  }
  (*fs) << std::endl;
}
#endif

// compute rotation of minimal angle that aligns vectors a & b:  b = R*a
vctRot3 XProdRotation(const vct3 &a, const vct3 &b)
{
  double EPS = 1.0e-6;

  // Find rotation that aligns a to b
  vct3 xProd = vctCrossProduct(a, b);
  double dProd = vctDotProduct(a, b);
  if (xProd.Norm() <= EPS)  // protect from divide by zero
  { // vectors are either already aligned or pointing opposite
    if (dProd > 0)
    {
      return vctRot3::Identity();
    }
    else
    {
      // need to rotate a by 180 degrees towards b
      int dim = 0;
      vct3 v = a;
      while (true)
      {
        // find any vector perpendicular to a & b
        v(dim) += 1.0;  // find a vector not parrallel to vector "a"
        v.NormalizedSelf();
        vct3 xProd2 = vctCrossProduct(a, v);
        if (xProd2.Norm() > EPS)
        {
          vct3 ax = xProd2.Normalized();
          vctAxAnRot3 R_AxAn(ax, cmnPI);
          return vctRot3(R_AxAn);
        }
      }
    }
  }
  else
  {
    // the norm of the cross product is the same for angles of x deg & x+180 deg
    //  between two vectors => use dot product to determine the angle
    //   NOTE: the angle corresponding to the cross product axis is always > 0;
    //         acos of the dot product gives the correct form
    //   NOTE: the problem with using norm of cross product isn't that we aren't
    //         going the right direction, but rather that we don't rotate far enough
    //         if A & B are seperated by more than 90 degrees.  I.e. if angular
    //         seperation between A & B is 100 degrees, then asin(norm(AxB)) gives
    //         the same angle as if A & B are seperated by 80 degrees => we don't
    //         know if the actual angle is X or 90+X using the norm of cross product.
    vct3 ax = xProd.Normalized();
    double an = acos(vctDotProduct(a, b));
    //double an = asin(t.Norm());
    vctAxAnRot3 R_AxAn(ax, an);
    return vctRot3(R_AxAn);
  }
}

#if 0
bool vct3HasValue(int x, vctFixedSizeVector<int, 3> vj)
{
  if (x == vj[0] || x == vj[1] || x == vj[2])
    return true;
  else 
    return false;
}

bool TrianglesAreNeighbors(vctFixedSizeVector<int, 3> vi, vctFixedSizeVector<int, 3> vj)
{
  // Two triangles are neighbors if any two vertices are
  //  held in common

  if (vct3HasValue(vi[0],vj))
  { // vertex 0 in common
    if (vct3HasValue(vi[1], vj))
    { // vertex 0 & 1 in common
      return true;
    }
    else if (vct3HasValue(vi[2], vj))
    { // vertex 0 & 2 in common
      return true;
    }
  }
  else if (vct3HasValue(vi[1], vj))
  { // vertex 1 in common
    if (vct3HasValue(vi[2], vj))
    { // vertex 1 & 2 in common
      return true;
    }
  }
  return false;
}


double ComputeAvgNeighborDistance(std::string meshFile)
{
  cisstMesh mesh;
  mesh.LoadMeshFile(meshFile);
  return ComputeAvgNeighborDistance(mesh);
}

// Compute average distance between neighboring triangle
//  center-points in a mesh
double ComputeAvgNeighborDistance( cisstMesh mesh )
{
  unsigned int numT = mesh.NumTriangles();
  double distSum = 0.0;
  unsigned int numNbrs = 0;

  // compute distances to triangle neighbors
  vctFixedSizeVector<int,3> vi,vj;  // vertex indices of triangles
  for (unsigned int i = 0; i < numT; i++)
  {
    mesh.TriangleVertexIndexes(i, vi[0], vi[1], vi[2]);
    unsigned int numLocalNbrs = 0;

    // find neighbors of triangle i
    for (unsigned int j = i + 1; j < numT; j++)
    {
      mesh.TriangleVertexIndexes(j, vj[0], vj[1], vj[2]);
      if (TrianglesAreNeighbors(vi, vj))
      {
        numNbrs++;
        distSum += (mesh.Triangles(i).Midpoint() - mesh.Triangles(j).Midpoint()).Norm();

        numLocalNbrs++;
        if (numLocalNbrs >= 3)
          break;  // a triangle cannot have more than 3 neighbors
      }
    }
  }

  return distSum / (double)numNbrs;
}
#endif