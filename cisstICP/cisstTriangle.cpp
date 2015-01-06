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


#include <cisstNumerical/nmrLSSolver.h>

#include "cisstTriangle.h"
#include "cisstMesh.h"


// Get vertex coordinates
vct3 cisstTriangle::VertexCoord(int i) const
{ return myMesh->VertexCoord(Vx[i]);
}

// Set triangle normal direction using the vertex coordinate
//   Note: assumes face direction follows right hand rule
//         with curl direction defined by Vx[0]->Vx[1]->Vx[2]
void cisstTriangle::SetNormalFromVertices()
{
  assert( Vx[0]>=0 && Vx[1]>=0 && Vx[2]>=0 );
	vct3 a = VertexCoord(0);
	vct3 b = VertexCoord(1);
	vct3 c = VertexCoord(2);
	vct3 n;
  n.CrossProductOf((b-a),(c-a));
  norm.Assign(n.Normalized());
}

// update plane params
//  (i.e. normal vector and distance from origin)
//  assumes face direction follows right hand rule
//   with curl direction defined by Vx[0]->Vx[1]->Vx[2]
void cisstTriangle::UpdatePlaneFromVertices()
{ 
	//if (Vx[0]<0||Vx[1]<0||Vx[2]<0) return;
  assert( Vx[0]>=0 && Vx[1]>=0 && Vx[2]>=0 );
	vct3 a = myMesh->VertexCoord(Vx[0]);
	vct3 b = myMesh->VertexCoord(Vx[1]);
	vct3 c = myMesh->VertexCoord(Vx[2]);
	vct3 n;
  n.CrossProductOf((b-a),(c-a));
  n.NormalizedSelf();
  norm.Assign(n);
  d = a.DotProduct(n);
  //vct3 n = (b-a)%(c-a);
	//n=n/n.Norm();
  //SetD(a*n);
}

// get triangle midpoint
vct3 cisstTriangle::Midpoint() const 
{
  return (myMesh->VertexCoord(Vx[0]) +
					myMesh->VertexCoord(Vx[1]) +
					myMesh->VertexCoord(Vx[2]))/3.0;
}

// get an interpolated position on the triangle face interpolated
//  from the vertex positions
vct3 cisstTriangle::Interpolate(double lam0, double lam1, double lam2) const
{ 
  return (myMesh->VertexCoord(Vx[0])*lam0) +
			   (myMesh->VertexCoord(Vx[1])*lam1) +
			   (myMesh->VertexCoord(Vx[2])*lam2);
} 

// get the center of a bouding sphere
//  SDB TODO: should also return sphere radius if actually using this.
vct3 cisstTriangle::EnclosingSphereCenter() const
{ 
	vct3 V0 = myMesh->VertexCoord(Vx[0]);
	vct3 q = myMesh->VertexCoord(Vx[1])-V0;
	vct3 r = myMesh->VertexCoord(Vx[2])-V0;
	vct3 qr = q-r;
	double q2=q*q;
	double r2=r*r;
	double qr2 = qr*qr;
	// pick the longest side
	// initially assume it is q

	if (q2>r2)
		{ // either q or qr
			if (q2<qr2)
				{ // it is qr.  So set "origin at vertex 1
					V0 =	myMesh->VertexCoord(Vx[1]);
					 q = myMesh->VertexCoord(Vx[2])-V0;
					 r = myMesh->VertexCoord(Vx[0])-V0;
					 r2 = r*r;
				};
		}
	else
		{ // either r or qr
			if (q2<qr2)
				{ // it is qr.  So set "origin at vertex 1
					V0 =	myMesh->VertexCoord(Vx[1]);
					 q = myMesh->VertexCoord(Vx[2])-V0;
					 r = myMesh->VertexCoord(Vx[0])-V0;
					 r2 = r*r;
				}
			else
				{ // it is r.  So set "origin at vertex 2
					V0 =	myMesh->VertexCoord(Vx[2]);
					 q = myMesh->VertexCoord(Vx[0])-V0;
					 r = myMesh->VertexCoord(Vx[1])-V0;
					 r2 = r*r;
				};
		};
	// here q = longest side (relocated to V0 origin)
	
	vct3 w = q/2.0;
	double w2=w*w;
	double wr=w*r;
	double alpha = wr/w2;
	double lambda = (r2-wr*alpha-2*wr)/2.0;
	if (lambda<0) 
		{ return V0+w;
		}
	else
		{ vct3 n = r-w*alpha;
			vct3 c = w+n*lambda;
			//JA debug
			vct3 tmp = V0+c;
			return tmp;
		}
}

//void cisstTriangle::Print(FILE* chan) const
//{ 
//  fprintf(chan,"(%4d,%4d,%4d) ",Vx[0],Vx[1],Vx[2]);
//	for (int i=0;i<3;i++) fprintfVct3Bracketed(chan,myMesh->VertexCoord(Vx[i]));
//}
