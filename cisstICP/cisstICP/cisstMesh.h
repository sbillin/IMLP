// ****************************************************************************
//
//    Copyright (c) 2015, Seth Billings, Russell Taylor, Johns Hopkins University
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
#include <ply_io.h>


class cisstMesh 
{

public:

  //--- Variables ---//

	vctDynamicVector<vct3>      vertices;       // the coordinates for each vertex in the mesh
  vctDynamicVector<vctInt3>   faces;          // the vertex indices for each triangle in the mesh
  vctDynamicVector<vct3>      faceNormals;    // the face normal for each triangle in the mesh

  // optional mesh properties
  vctDynamicVector<vct3>      vertexNormals;  // a normal orientation associated with each vertex  
  vctDynamicVector<vctInt3>   faceNeighbors;  // the face indices for the neighbors of each triangle
                                              //  in the mesh

  // mesh noise model
  //  NOTE: if used, this must be set manually by the user AFTER loading the mesh file
  //        (defaults to all zeroes, i.e. zero measurement noise on the mesh)
  vctDynamicVector<vct3x3>  TriangleCov;        // triangle covariances
  vctDynamicVector<vct3>    TriangleCovEig;     // triangle covariance eigenvalues (in decreasing size)

private:

  ply_io ply_obj; 

public:

  //--- Methods ---//

  // constructor
  cisstMesh() {};

  // destructor
  ~cisstMesh() {}

  // initializes all mesh properties to empty (default initializer);
  //  this is a useful routine to use while building a mesh,
  //  since some mesh properties are optional and may not be
  //  initialized by the data used to build the mesh; calling this 
  //  ensures that unused properties are emptied rather than left with
  //  possibly invalid values
  void ResetMesh();

  inline int NumVertices() const { return vertices.size(); }
  inline int NumTriangles() const { return faces.size(); }

  // initializes triangle noise models to zero (default initializer)
  void InitializeNoiseModel();

  // computes noise model covariances for each triangle in the mesh
  //  such that the in-plane and perpendicular-plane directions have
  //  the specified variance
  void InitializeNoiseModel( double noiseInPlaneVar, double noisePerpPlaneVar);

  void SaveTriangleCovariances(std::string &filePath);

 // // get mesh vertex index for the given triangle index and vertex
	//inline int TriangleVertexIndex(int ti, int v)
 //   { return faces[ti][v]; 
 //   }

	//// get mesh vertex indexes for the given triangle index
	//inline virtual void TriangleVertexIndexes(int ti, int &v0, int &v1, int &v2) const
	//	{ v0=Triangles[ti].Vx[0];
	//		v1=Triangles[ti].Vx[1];
	//		v2=Triangles[ti].Vx[2];
	//	}

  // get coordinates of all three vertices for a given face index
	inline void FaceCoords(int ti, vct3 &v0, vct3 &v1, vct3 &v2) const
	{ v0 = vertices[faces[ti][0]];
    v1 = vertices[faces[ti][1]];
    v2 = vertices[faces[ti][2]];
	}
  // get vertex coordinate for a given face/vertex index
	inline vct3& FaceCoord(int ti, int vi)
  { return vertices[faces[ti][vi]];
  }

  // assumes vertex order follows right-hand rule with curl v1->v2->v3
  void ComputeFaceNormalsFromVertices();

 // // set vertex coords for a given vertex index
 // //  Note: changing vertex coordinates changes the normal vectors
 // //        of all associated triangles
 // //        it is responsibility of the user to keep these values
 // //        consistent
	//inline void SetVertexCoord(int vx, vct3 v)	{VertexCoordinates[vx].Assign(v);}
	//inline void SetVertexCoord(int vx, float x, float y, float z) {VertexCoordinates[vx].Assign(x,y,z);}
	//inline void SetVertexCoord(int vx, double x, double y, double z) {VertexCoordinates[vx].Assign(x,y,z);}


  // Mesh I/O

  // Build mesh from data arrays
  int  LoadMesh(
    const vctDynamicVector<vct3> *vertices,
    const vctDynamicVector<vctInt3> *faces,
    const vctDynamicVector<vct3> *face_normals = NULL, 
    const vctDynamicVector<vctInt3> *face_neighbors = NULL,
    const vctDynamicVector<vct3> *vertex_normals = NULL
    );

  // Load mesh from PLY file
  void LoadPLY(const std::string &input_file);

  // Save mesh to PLY fle
  void SavePLY(const std::string &output_file);


//  // -- Deprecated I/O --
//
//  // Build mesh from an array of vertices, faces, and face normals
//  int  LoadMesh(
//    const vctDynamicVector<vct3> &V,
//    const vctDynamicVector<vctInt3> &T,
//    const vctDynamicVector<vct3> &N);
//
//  // Build mesh from an array of vertices and faces;
//  // face normals are computed from the vertex positions assuming
//  // vertex order follows the right-hand rule relative to the face normal
//  int  LoadMesh(
//    const vctDynamicVector<vct3> &V,
//    const vctDynamicVector<vctInt3> &T);
//
//  // Build new mesh from a single .mesh file
//  int  LoadMeshFile( const std::string &meshFilePath );
//
//  // Build new mesh from multiple .mesh files
//  int  LoadMeshFileMultiple( const std::vector<std::string> &meshFilePaths ); 
//
//  // Save mesh to .mesh file
//  int  SaveMeshFile( const std::string &filePath );
//
//  // Build new mesh from .stl file
//  int  LoadMeshFromSTLFile(const std::string &stlFilePath);
//
//private:
//
//  // Load .mesh file, adding it to the current mesh while preserving all
//  //  data currently existing in the mesh
//  int  AddMeshFile(const std::string &meshFilePath);
//
//public:
//
//  // Legacy Mesh I/O
//	void ReadMeshFile(const char *fn);
//  void ReadMeshMeshFile( const std::string &meshFilePath );
//	void ReadSURMeshFile(const char *fn);
//	void ReadSFCMeshFile(const char *fn);
//	void WriteMeshFile(const char *fn);
//	void WriteSURMeshFile(const char *fn);
};

#endif // _cisstMesh_h_
