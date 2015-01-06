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
#include "cisstTriangleClosestPointSolver.h"

#include <fstream>


vct3 cisstTriangleClosestPointSolver::ProjectOnSegment(const vct3& c,const vct3& p, const vct3& r) 
{
	vct3 pc=c-p;
	vct3 pr=r-p;
	double lam = (pc*pr)/(pr*pr);
	if (lam<=0.0) { return p;};
	if (lam>1.0) {return r;};
	return p+pr*lam;
};

cisstTriangleClosestPointSolver::cisstTriangleClosestPointSolver() 
			: A(3,2,VCT_COL_MAJOR), h(3), P(3), g(3), b(3), x(3), 
			  LeastSquaresSolver(3,2,1,VCT_COL_MAJOR), 
			  B(3,1,VCT_COL_MAJOR) 
{};


void cisstTriangleClosestPointSolver::SolveLamMuNu(const vct3& a, 
													const vct3& p,
													const vct3& q,
													const vct3& r,
													double &lambda,
													double &mu,
													double &nu)
{ vct3 pa=a-p; // b(0)=pa.x; b(1)=pa.y; b(2)=pa.z;
  vct3 pq=q-p; // A(0,0)=pq.x; A(1,0)=pq.y; A(2,0)=pq.z;
  vct3 pr=r-p; // A(0,1)=pr.x; A(1,1)=pr.y; A(2,1)=pr.z;
  for (int i=0;i<3;i++) 
	{ b(i)=pa(i); A(i,0)=pq(i); A(i,1)=pr(i); 
	  B(i,0)=b(i);  // because using dumb interface to solver
	};
  LeastSquaresSolver.Solve(A,B);		   // replace with HFTI call when can do so
  lambda=B(0,0); mu=B(1,0);				   // ditto
  // nmrAlgorithmHFTI(A,h,P,Tau,g,b,x,0);  // x := lambda and mu
  // lambda = x(0);
  // mu = x(1);
  nu = 1.-lambda-mu;
};

// TODO: can this algorithm fail for obtuse triangle if both mu & lambda < 0?
//
// rht this function was inlined
int  cisstTriangleClosestPointSolver::FindClosestPointOnTriangle(	const vct3& a, 
																	const vct3& p,
																	const vct3& q,
																	const vct3& r,
																	double distBound, // -1 if ignore
																	vct3& ret)
	{ 
		vct3 pa=a-p; // b(0)=pa.x; b(1)=pa.y; b(2)=pa.z;
		vct3 pq=q-p; // A(0,0)=pq.x; A(1,0)=pq.y; A(2,0)=pq.z;
		vct3 pr=r-p; // A(0,1)=pr.x; A(1,1)=pr.y; A(2,1)=pr.z;
		for (int i=0;i<3;i++) 
			{ b(i)=pa(i); A(i,0)=pq(i); A(i,1)=pr(i); 
			B(i,0)=b(i);  // because using dumb interface to solver
			};		
		LeastSquaresSolver.Solve(A,B);		     // replace with HFTI call when can do so
		double lambda=B(0,0); double mu=B(1,0);	 // ditto
		// nmrAlgorithmHFTI(A,h,P,Tau,g,b,x,0);  // x := lambda and mu		double lambda = x(0);
		// double lambda=x(0);
		// double mu = x(1);
		vct3 c = p+pq*lambda+pr*mu;
		if (distBound >= 0.0) 
		{ double dist = (c-a).Norm();
			if (dist>distBound) 
			{ ret=c; return 0; };
		};
		if (lambda<=0.0)		{ ret =  ProjectOnSegment(c,r,p); return 1; };
		if (mu <=0.0)			{ ret =  ProjectOnSegment(c,p,q); return 2; };
		if (lambda+mu> 1.0 )	{ ret =  ProjectOnSegment(c,q,r); return 3; };
		ret= c; return 4;
  };

