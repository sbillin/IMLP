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

#include <algorithm>
#include <limits>

#include "algDirICP_PIMLOP.h"
#include "DirPDTreeNode.h"
#include "utilities.h"


#include <assert.h>
#undef NDEBUG       // enable assert in release mode

#define EPS  1.0e-14
//#define DEBUG_DISABLE_ORIENTATION_OPT


// Constructor
algDirICP_PIMLOP::algDirICP_PIMLOP(
  DirPDTreeBase *pDirTree,
  const vctDynamicVector<vct3> &samplePts,
  const vctDynamicVector<vct2> &sampleNorms2d,
  const vctDynamicVector<vctRot3> &Rx_pln,
  const vctDynamicVector<double> &sample_k,
  const vctDynamicVector<vct3x3> &sample_M)
  : algDirICP(pDirTree, samplePts, vctDynamicVector<vct3>(samplePts.size(),vct3(0.0))),
  //algDirPDTree_vonMisesPrj(pDirTree),
  dlib(this)
{
  SetSamples(samplePts, sampleNorms2d, Rx_pln, sample_k, sample_M);
}

double algDirICP_PIMLOP::ICP_EvaluateErrorFunction()
{
  // Just return match error for now
  //  optionally, could also add the log term out front
  double Error = 0.0;
  for (unsigned int i = 0; i < nSamples; i++)
  {
    Error += MatchError(samplePtsXfmd[i], sampleNorms2d[i],
      matchPts[i], matchNorms[i], Ry_pln[i],
      k[i], R_invM_Rt[i]);
  }
  return Error;
}

vctFrm3 algDirICP_PIMLOP::ICP_RegisterMatches()
{
  //vctRodRot3 R0( Fact.Rotation() );
  //debugStream << "F0:  rod: " << R0 << " trans: " << Fact.Translation() << std::endl;

  //for (unsigned int i=0; i<2; i++)
  //{
  //  debugStream << " " << i << " Xp: " << pICP->TransformedSamplePtsSets[0](i) <<
  //    " Xn: " << pICP->TransformedSampleNormSets[0](i) << std::endl;
  //  debugStream << " " << i << " Yp: " << pICP->ClosestPointSets[0](i) <<
  //    " Yn: " << pICP->ClosestPointNormSets[0](i) << std::endl;
  //}

  vctFrm3 dF;
  vct6 x0(0.0);
  vct6 x;

  // x_prev must begin at a different value than x0
  x_prev.SetAll(std::numeric_limits<double>::max());

  //x = gsl.ComputeRegistration( x0 );
  x = dlib.ComputeRegistration(x0);
  //debugStream << "x0: " << x0 << std::endl;
  //debugStream << "xf: " << x << std::endl;

  // update transform
  vctFixedSizeVectorRef<double, 3, 1> alpha(x, 0);
  vctFixedSizeVectorRef<double, 3, 1> t(x, 3);
  dF.Rotation() = vctRot3(vctRodRot3(alpha));
  dF.Translation() = t;
  Freg = dF * Freg;

  return Freg;
}

void algDirICP_PIMLOP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
  //std::string debugFile = "C:/Code/Repos_Git/RegistrationTools/ICP_TestData/LastRun/KentDebug.txt";
  //debugStream.open( debugFile.c_str() );

  // initialize base class
  algDirICP::ICP_InitializeParameters(FGuess);

  // monitoring variables
  errFuncNormWeight = 0.0;
  errFuncPosWeight = 0.0;

  if (k.size() != nSamples || M.size() != nSamples || invM.size() != nSamples 
    || N.size() != nSamples || invN.size() != nSamples || Dmin.size() != nSamples
    || R_invM_Rt.size() != nSamples || N_Rt.size() != nSamples || inv_N_Rt.size() != nSamples
    || Yp_t.size() != nSamples || Rat_Yp_RaXp_t.size() != nSamples //|| Yp_RaXp_t.size() != nSamples 
    || invM_Rat_Yp_RaXp_t.size() != nSamples
    || Yprj.size() != nSamples || Ypln.size() != nSamples || Ynorm.size() != nSamples
    || sampleNorms2d.size() != nSamples || Rx_pln.size() != nSamples 
    || Ry_pln.size() != nSamples)
  {
    std::cout << "ERROR: noise parameters not properly initialized" << std::endl;
    assert(0);
  }

  // initialize rotation dependent parameters to prepare for initial match
  UpdateNoiseModel_SamplesXfmd(FGuess);
}


void algDirICP_PIMLOP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg)
{
  // base class
  algDirICP::ICP_UpdateParameters_PostRegister(Freg);

  UpdateNoiseModel_SamplesXfmd(Freg);
}

unsigned int algDirICP_PIMLOP::ICP_FilterMatches()
{
  return 0;
}


// Helper Methods

void algDirICP_PIMLOP::UpdateNoiseModel_SamplesXfmd(vctFrm3 &Freg)
{
  // update noise models of the transformed sample points
  vctRot3 R(Freg.Rotation());
  for (unsigned int s = 0; s < nSamples; s++)
  {
    Ry_pln[s] = R * Rx_pln[s];

    R_invM_Rt[s] = R * invM[s] * R.Transpose();
    N_Rt[s] = N[s] * R.Transpose();
    inv_N_Rt[s] = R*invN[s];
  }
}

void algDirICP_PIMLOP::SetSamples(
  const vctDynamicVector<vct3> &samplePts,
  const vctDynamicVector<vct2> &sampleNorms2d,
  const vctDynamicVector<vctRot3> &Rx_pln,
  const vctDynamicVector<double> &sample_k,
  const vctDynamicVector<vct3x3> &sample_M)
{
  // base class
  nSamples = samplePts.size();
  vctDynamicVector<vct3> dummy3dNorms(nSamples, vct3(0.0));
  algDirICP::SetSamples(samplePts, dummy3dNorms);
  //this->samplePts = samplePts;

  // other sample data
  this->sampleNorms2d = sampleNorms2d;
  this->Rx_pln = Rx_pln;

  Ry_pln.SetSize(nSamples);

  // sample noise model  
  this->k = sample_k;
  this->M = sample_M;
  invM.SetSize(nSamples);
  N.SetSize(nSamples);
  invN.SetSize(nSamples);
  Dmin.SetSize(nSamples);
  R_invM_Rt.SetSize(nSamples);
  N_Rt.SetSize(nSamples);
  inv_N_Rt.SetSize(nSamples);
  //Emin.SetSize(nSamples);

  // optimizer variables
  Yp_t.SetSize(nSamples);
  Rat_Yp_RaXp_t.SetSize(nSamples);
  //Yp_RaXp_t.SetSize(nSamples);
  invM_Rat_Yp_RaXp_t.SetSize(nSamples);
  Yprj.SetSize(nSamples);
  Ypln.SetSize(nSamples);
  Ynorm.SetSize(nSamples);

  if (nSamples > 0)
  {
    ComputeNoiseModelDecompositions();
  }
}

// Note: this function depends on the SamplePreMatch() function
//       to set the noise model of the current transformed sample
//       point before this function is called
void algDirICP_PIMLOP::ComputeNoiseModelDecompositions()
{
  if (nSamples <= 0 || M.size() != nSamples)
  {
    std::cout << "ERROR: noise model has different # of elements than samples or no samples exist"
      << std::endl;
    assert(0);
  }

  for (unsigned int i = 0; i < nSamples; i++)
  {
    algDirPDTree_vonMisesPrj::ComputeCovDecomposition(
      M[i], invM[i], N[i], invN[i], Dmin[i]);
  }
}


void algDirICP_PIMLOP::UpdateOptimizerCalculations(const vct6 &x)
{
  // TODO: speed-up by checking for identity matrix (i.e. zero-vector of x)
  a.Assign(x[0], x[1], x[2]);
  t.Assign(x[3], x[4], x[5]);

  // incremental rotation matrix
  Ra = vctRot3(vctRodRot3(a));

  vctDynamicVectorRef<vct3>   Xp_xfm(samplePtsXfmd);
  vctDynamicVectorRef<vct3>   Yp(matchPts);  
  vctDynamicVectorRef<vct3>   Yn(matchNorms);  
  vctRot3 Ra_pln;
  for (unsigned int i = 0; i < nSamples; i++)
  {    
    //--- orientation ---//
    Ra_pln = Ra * Freg.Rotation() * Rx_pln.Element(i);
    //RaR = RaRreg * Rx_pln.Element(i);
    
    Yprj.Element(i) = (Ra_pln.TransposeRef() * Yn.Element(i)).XY();
    Ynorm.Element(i) = Yprj.Element(i).Norm();
    if (Ynorm.Element(i) < 0.001)
    {
      Ypln.Element(i).Assign(vct2(0.0));  // makes Jacobian for this term = 0 (rather than approaching inf)
      Ynorm.Element(i) = 0.001;           // protects from division by zero
    }
    else
    {
      Ypln.Element(i) = Yprj.Element(i) / Ynorm.Element(i);
    }

    //--- position ---//
    Yp_t.Element(i) = Yp.Element(i) - t;
    Rat_Yp_RaXp_t.Element(i) = Ra.TransposeRef() * Yp_t.Element(i) - Xp_xfm.Element(i);
    invM_Rat_Yp_RaXp_t.Element(i) = R_invM_Rt.Element(i) * Rat_Yp_RaXp_t.Element(i);
    //Yp_RaXp_t.Element(i) = Yp.Element(i) - Ra*Xp_xfm.Element(i) - t;
    //invM_Rat_Yp_RaXp_t.Element(i) = R_invM_Rt.Element(i) * (Ra.TransposeRef() * Yp_RaXp_t.Element(i));
    //invM_Yp_RaXp_t.Element(i) = R_invM_Rt.Element(i) * Yp_RaXp_t.Element(i);
  }
  x_prev = x;
}

double algDirICP_PIMLOP::CostFunctionValue(const vct6 &x)
{
  // don't recompute these if already computed for gradient
  if (x.NotEqual(x_prev))
  {
    UpdateOptimizerCalculations(x);
  }

  double f = 0.0;
  vctDynamicVectorRef<vct2> Xpln(sampleNorms2d);  
  for (unsigned int i = 0; i < nSamples; i++)
  {
#ifndef DEBUG_DISABLE_ORIENTATION_OPT
    // NOTE: Ypln has already been transformed to local Xpln coordinates
    f += k[i] * (1.0 - Ypln[i] * Xpln[i])
      + (Rat_Yp_RaXp_t[i] * invM_Rat_Yp_RaXp_t[i]) / 2.0;
      //+ (Yp_RaXp_t[i] * Ra * invM_Rat_Yp_RaXp_t[i]) / 2.0;
#else
    f += (Yp_RaXp_t[i] * Ra * invM_Rat_Yp_RaXp_t[i]) / 2.0;
#endif
  }

  return f;
}


void algDirICP_PIMLOP::CostFunctionGradient(const vct6 &x, vct6 &g)
{
  vct3 gTheta;
  vctFixedSizeVector<vctRot3, 3> dRa;  // Rodrigues Jacobians of R(a) wrt ax,ay,az

  // don't recompute these if already computed for cost function value
  if (x.NotEqual(x_prev))
  {
    UpdateOptimizerCalculations(x);
  }

  ComputeRodriguesJacobians(a, dRa);
  
  // form the cost function gradient
  g.SetAll(0.0);
  vctFixedSizeVectorRef<double, 3, 1> ga(g, 0);
  vctFixedSizeVectorRef<double, 3, 1> gt(g, 3);
  vctDynamicVectorRef<vct3>   Xp_xfm(samplePtsXfmd);
  vctDynamicVectorRef<vct2>   Xpln(sampleNorms2d);
  vctDynamicVectorRef<vct3>   Yp(matchPts);
  vctDynamicVectorRef<vct3>   Yn(matchNorms);
  vct2x2 Jpln_prj;
  vct2x3 Jprj_a;
  vct2x2 YYt;
  vct3 dRaxYn, dRayYn, dRazYn;
  //vct3 dRaxXp, dRayXp, dRazXp;
  //vct3 dRatx_Yp_RaXp_t, dRaty_Yp_RaXp_t, dRatz_Yp_RaXp_t;
  vct3x3 Jz_a;
  for (unsigned int s = 0; s < nSamples; s++)
  {
#ifndef DEBUG_DISABLE_ORIENTATION_OPT
    // von-Mises gradient term (orientation)
    //  rotational effect
    YYt.OuterProductOf(Yprj[s], Yprj[s]);
    Jpln_prj = vct2x2::Eye() / Ynorm[s] - YYt / pow(Ynorm[s], 3.0);
    dRaxYn = dRa[0].TransposeRef() * Yn[s];
    dRayYn = dRa[1].TransposeRef() * Yn[s];
    dRazYn = dRa[2].TransposeRef() * Yn[s];
    Jprj_a.Column(0) = (Ry_pln[s].TransposeRef() * dRaxYn).XY();
    Jprj_a.Column(1) = (Ry_pln[s].TransposeRef() * dRayYn).XY();
    Jprj_a.Column(2) = (Ry_pln[s].TransposeRef() * dRazYn).XY();
    ga -= k[s] * Xpln[s] * (Jpln_prj * Jprj_a);
#endif

    // Gaussian gradient term (positions)
    //  rotational effect
    Jz_a.Column(0) = dRa[0].TransposeRef() * Yp_t[s];
    Jz_a.Column(1) = dRa[1].TransposeRef() * Yp_t[s];
    Jz_a.Column(2) = dRa[2].TransposeRef() * Yp_t[s];
    ga += invM_Rat_Yp_RaXp_t[s] * Jz_a;
    //dRaxXp = dRa[0] * Xp_xfm[s];
    //dRayXp = dRa[1] * Xp_xfm[s];
    //dRazXp = dRa[2] * Xp_xfm[s];
    //dRatx_Yp_RaXp_t = dRa[0].TransposeRef() * Yp_RaXp_t[s];
    //dRaty_Yp_RaXp_t = dRa[1].TransposeRef() * Yp_RaXp_t[s];
    //dRatz_Yp_RaXp_t = dRa[2].TransposeRef() * Yp_RaXp_t[s];
    //Jz_a.Column(0) = -Ra.TransposeRef() * dRaxXp + dRatx_Yp_RaXp_t;
    //Jz_a.Column(1) = -Ra.TransposeRef() * dRayXp + dRaty_Yp_RaXp_t;
    //Jz_a.Column(2) = -Ra.TransposeRef() * dRazXp + dRatz_Yp_RaXp_t;
    //ga += invM_Rat_Yp_RaXp_t[s] * Jz_a;
    // translational effect
    gt -= Ra * invM_Rat_Yp_RaXp_t[s];
  }
}

// Xp    ~ transformed source point
// Xpln  ~ non-transformed source orientation in 2d plane coordinates
//         where Xn (3d orientation in local coors) = Rx_pln * [Xpln; 0]
// Yp    ~ target point
// Yn    ~ target 3d orientation
// Ry_pln  ~ Yprj = Pxy( Ry_pln' * Yn)  (projection of Yn to plane containing Xpln)
// k     ~ orientation concentration
// invM  ~ positional covariance for transformed sample
double algDirICP_PIMLOP::MatchError(
  const vct3 &Xp, const vct2 &Xpln,
  const vct3 &Yp, const vct3 &Yn,
  const vctRot3 &Ry_pln,
  double k, const vct3x3 &invM)
{
  // project y orientation to 2d plane coords
  vct2 Ypln;
  vct2 Yprj = (Ry_pln.TransposeRef() * Yn).XY();
  double Ynorm = Yprj.Norm();
  if (Ynorm < 0.001)
  {
    Ypln.SetAll(0.0);   // effectively sets orientation error to 90 deg.
  }
  else
  {
    Ypln = Yprj / Ynorm;
  }

  // add extra k*1 term so that match error >= 0
#ifndef DEBUG_DISABLE_ORIENTATION_OPT
  return k*(1 - Ypln*Xpln) + ((Yp - Xp)*invM*(Yp - Xp)) / 2.0;
#else
  return ((Yp - Xp)*invM*(Yp - Xp)) / 2.0;
#endif
}
