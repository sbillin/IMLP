
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

#ifndef _cisstTriangle_h_
#define _cisstTriangle_h_

#include <stdio.h>
#include <cisstVector.h>
#include <cisstCommon.h>

#include "cisstBoundingBox.h"

class	cisstMesh;  // forward decleration for mutual dependency


class cisstTriangle
{
	friend class cisstMesh;

private:

  cisstMesh *myMesh;

public:

	int     Vx[3];  // vertex indices (use int type so we have -1 for unitialized index)
	vct3    norm;   // normal to triangle plane

  //  NOTE: referencing this bounding box takes more time than building
  //        a new one from scratch in the node proximity tests for Std ICP
  cisstBoundingBox BB;  // bounding box around triangle

  // constructors
  //
  //   NOTE: as the mesh may still be under construction, we can't 
  //         be certain that the vertex coordinates have been initialized 
  //         yet => don't set plane parmeters (norm & d) using vertex
  //         coordinates in any of the constructors
  //       
  cisstTriangle( cisstMesh *M )
  { assert(M);
    myMesh = M;
    Vx[0]=Vx[1]=Vx[2] = -1;
    norm.SetAll(0.0);
  }
	cisstTriangle( cisstMesh *M, 
							   int v0, int v1, int v2 ) 
	{ assert(M);
    myMesh = M;
    SetVertexIndexes(v0,v1,v2);
    norm.SetAll(0.0);
  }
	cisstTriangle( cisstMesh *M, 
				         int v0, int v1, int v2,
                 const vct3 &normal) 
	{ assert(M);
    myMesh = M;
    SetVertexIndexes(v0,v1,v2);
    norm = normal.Normalized();
  }
  // don't use this constructor
  //  (it is here only to support creating dynamic vectors of this class type)
  cisstTriangle()
  { myMesh = 0;
    Vx[0]=Vx[1]=Vx[2] = -1;
    norm.SetAll(0.0);
  }

  // Get vertex coordinates
  inline vct3 VertexCoord(int i) const;

  // set vertex indexes
  void SetVertexIndexes(int vx0, int vx1, int vx2)
  { Vx[0] = vx0; Vx[1] = vx1; Vx[2] = vx2;
    ComputeBoundingBox();
  }

  // compute bounding box for the current vertices
  void ComputeBoundingBox()
  { BB.Include(VertexCoord(0));
    BB.Include(VertexCoord(1));
    BB.Include(VertexCoord(2));
  }

  // Plane params
  void SetNormalFromVertices();

	vct3 Midpoint() const;
	vct3 Interpolate(double lam0, double lam1, double lam2) const;
	vct3 EnclosingSphereCenter() const;


  //-- Legacy Code --//

  int Nx[3];    // neighbor triangle indices
	double  d;	  // distance from triangle plane to origin
                //  vector plane equation: N*x = d

	int VertexIndex(int i) const
	{	assert(i>=0&&i<3);
		return Vx[i];
	}
	int& VertexIndex(int i)
	{	assert(i>=0&&i<3);
		return Vx[i];
	};

  double FaceD() const	{ return d; };

  cisstTriangle( cisstMesh *M, 
    				           int v0, int v1, int v2,
                       int n0, int n1, int n2 )
  { Vx[0]=v0; Vx[1]=v1; Vx[2]=v2;
    UpdatePlaneFromVertices();
    Nx[0]=n0; Nx[1]=n1; Nx[2]=n2;
    myMesh = M;
  };
  int& NeighborIndex (int i)
  {	assert(i>=0&&i<3);
	  // cisstGENERIC_ASSERT_WITH_MSG((i>=0&&i<3), "SubscriptError in indexing cisstTriangle");
	  return Nx[i];
  };
  // get neighbor index
  int NeighborIndex (int i) const
  {	assert(i>=0&&i<3);
    return Nx[i];
  };
  void UpdatePlaneFromVertices();

  //void Print(FILE* chan) const;

};

#endif // _cisstTriangle_h_

// ****************************************************************************
//                              Change History
// ****************************************************************************
//
//
// ****************************************************************************
