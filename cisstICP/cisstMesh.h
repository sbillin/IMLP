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
#ifndef _cisstMesh_h_
#define _cisstMesh_h_

#include <stdio.h>
#include <cisstVector.h>
#include <cisstCommon.h>
#include "cisstTriangle.h"

class cisstMesh 
{

protected: 

	void Reset();

	int cisstError(const char* Msg) const 
	{ assert(0);
	  // cisstGENERIC_ASSERT_WITH_MSG(0,Msg); 
    std::cerr << Msg << std::endl;  // SDB
	  return 0;
  };

public:

	vctDynamicVector<vct3>          VertexCoordinates;
	vctDynamicVector<cisstTriangle> Triangles;

  // Noise model of the mesh
  //  NOTE: if used, this must be set manually by the user AFTER loading the mesh file
  //        (defaults to all zeroes, i.e. zero measurement noise on the mesh)
  vctDynamicVector<vct3x3>  TriangleCov;        // triangle covariances
  vctDynamicVector<vct3>    TriangleCovEig;     // triangle covariance eigenvalues (in decreasing size)


  // constructors
  cisstMesh() {};
  cisstMesh(const char *fn) { ReadMeshFile(fn); }

  // destructor
  ~cisstMesh() { Reset(); };

	inline int NumTriangles() const {return Triangles.size();}
	inline int NumVertices() const  {return VertexCoordinates.size();}

  void ComputeTriangleNoiseModels(
    double noiseInPlaneVar,
    double noisePerpPlaneVar);

  void SaveTriangleCovariances(std::string &filePath);

  //// get noise model for a given triangle index
  //inline vct3x3& TriangleCov( int ti )
  //  { return TriangleCov[ti];
  //  }

  // get triangle at a given triangle index
	inline cisstTriangle& Triangle(int ti)
    { return Triangles[ti]; 
    }

  // get mesh vertex index for the given triangle index and vertex
	inline int TriangleVertexIndex(int ti, int v)
    { return Triangles[ti].VertexIndex(v); 
    }

	// get mesh vertex indexes for the given triangle index
	inline virtual void TriangleVertexIndexes(int ti, int &v0, int &v1, int &v2) const
		{ v0=Triangles[ti].Vx[0];
			v1=Triangles[ti].Vx[1];
			v2=Triangles[ti].Vx[2];
		}

  // get triangle normal for the given triangle index
  inline vct3 TriangleNorm(int ti)
    { return Triangles[ti].norm;
    }

  // get coordinates of all vertices for a given triangle index  (optimized)
	inline void VerticesCoords(int ti, vct3 &v0, vct3 &v1, vct3 &v2) const
	{ v0 = VertexCoordinates.Element(Triangles[ti].Vx[0]);
    v1 = VertexCoordinates.Element(Triangles[ti].Vx[1]);
    v2 = VertexCoordinates.Element(Triangles[ti].Vx[2]);
	}
  // get vertex coordinate for a given vertex index
	inline vct3 VertexCoord(int vx) const
  { return VertexCoordinates.at(vx);
  }

  // set vertex coords for a given vertex index
  //  Note: changing vertex coordinates changes the normal vectors
  //        of all associated triangles
  //        it is responsibility of the user to keep these values
  //        consistent
	inline void SetVertexCoord(int vx, vct3 v)	{VertexCoordinates[vx].Assign(v);}
	inline void SetVertexCoord(int vx, float x, float y, float z) {VertexCoordinates[vx].Assign(x,y,z);}
	inline void SetVertexCoord(int vx, double x, double y, double z) {VertexCoordinates[vx].Assign(x,y,z);}

  // Mesh I/O
  int  LoadMesh(
    const vctDynamicVector<vct3> &V,
    const vctDynamicVector<vctInt3> &T,
    const vctDynamicVector<vct3> &N);
  int  LoadMeshFromSTLFile(const std::string &stlFilePath);
  //int  LoadMeshFromLegacyVTKFile(const std::string &vtkFilePath);
  int  LoadMeshFile( const std::string &meshFilePath );
  int  LoadMeshFileMultiple( const std::vector<std::string> &meshFilePaths );  
  int  SaveMeshFile( const std::string &filePath );
private:
  int  AddMeshFile(const std::string &meshFilePath);
public:

  // Legacy Mesh I/O
	void ReadMeshFile(const char *fn);
  void ReadMeshMeshFile( const std::string &meshFilePath );
	void ReadSURMeshFile(const char *fn);
	void ReadSFCMeshFile(const char *fn);
	void WriteMeshFile(const char *fn);
	void WriteSURMeshFile(const char *fn);

	//void Print(FILE* chan);

};

#endif // _cisstMesh_h_