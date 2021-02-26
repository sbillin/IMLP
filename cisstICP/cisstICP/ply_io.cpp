// ****************************************************************************
//
//    Copyright (c) 2016, Seth Billings, Johns Hopkins University
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

#include "ply_io.h"
#include "rply.h"

// temp global variables used when loading a ply file
int _nVertices;
int _nFaces;
int _nPoints;
int _idx_v;
int _idx_vn;
int _idx_f;
int _idx_fn;
int _idx_fnbr;
int _idx_p;
int _idx_pn;
int _listLength;
int _v[3];
vctDynamicVector<vct3>      _vertices;
vctDynamicVector<vct3>      _vertex_normals;
vctDynamicVector<vctInt3>   _faces;
vctDynamicVector<vct3>      _face_normals;
vctDynamicVector<vctInt3>   _face_neighbors;
vctDynamicVector<vct3>      _points_3d;
vctDynamicVector<vct3>      _point_normals_3d;
vctDynamicVector<vct2>      _points_2d;
vctDynamicVector<vct2>      _point_normals_2d;


static void init_temp_variables_mesh(
  int num_vertices, 
  int num_faces,
  int num_face_normals,
  int num_face_neighbors,
  int num_vertex_normals)
{
  _nVertices = num_vertices;
  _nFaces = num_faces;

  _idx_v = 0;
  _idx_f = 0;
  _idx_fn = 0;
  _idx_fnbr = 0;
  _idx_vn = 0;

  _vertices.SetSize(num_vertices);
  _faces.SetSize(num_faces);
  _face_normals.SetSize(num_face_normals);
  _face_neighbors.SetSize(num_face_neighbors);
  _vertex_normals.SetSize(num_vertex_normals);
}

static void init_temp_variables_pointcloud3d(
  int num_points, 
  int num_point_normals)
{
  _nPoints = num_points;

  _idx_p = 0;
  _idx_pn = 0;

  _points_3d.SetSize(num_points);
  _point_normals_3d.SetSize(num_point_normals);
}

static void init_temp_variables_pointcloud2d(
  int num_points,
  int num_point_normals)
{
  _nPoints = num_points;

  _idx_p = 0;
  _idx_pn = 0;

  _points_2d.SetSize(num_points);
  _point_normals_2d.SetSize(num_point_normals);
}

static int vertex_cb( p_ply_argument argument ) 
{
  long param_idx;

  // get whether this parameter value is x, y, or z
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_v >= _nVertices) {
    std::cout << "ERROR: vertex count exceeded" << std::endl;
    return 0;
  }
  _vertices(_idx_v)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 2) _idx_v += 1;

  return 1;
}

static int vertex_normal_cb(p_ply_argument argument)
{
  long param_idx;

  // and whether this element value is x, y, or z
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_vn >= _nVertices) {
    std::cout << "ERROR: vertex count exceeded while loading vertex normals" << std::endl;
    return 0;
  }
  _vertex_normals(_idx_vn)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 2) _idx_vn += 1;

  return 1;
}

static int face_cb( p_ply_argument argument ) 
{
  long length, param_idx, vertex_value;

  // get whether this parameter value is list length or vertex index
  ply_get_argument_property(argument, NULL, &length, &param_idx);

  switch (param_idx) {
  
  case -1:  // list length
    _listLength = (int)ply_get_argument_value(argument);
    break;
  
  default:  // vertex index
    if (param_idx >= _listLength) {
      std::cout << "ERROR: face parameter exceeds the list length" << std::endl;
      return 0;
    }
    if (_idx_f >= _faces.size()) {
      // face count may be exceeded if more than one face is provided per face list
      _faces.resize(_faces.size()*2);
      std::cout << "WARNING: face count exceeded, resizing face vector" << std::endl;
    }
    // convert vertex index list to triangle connectivity
    vertex_value = (int)ply_get_argument_value(argument);
    if (param_idx <= 2) {
      _faces(_idx_f)(param_idx) = vertex_value;
      _v[param_idx] = vertex_value;
      if (param_idx == 2) _idx_f += 1;  // next face
    } else {
      // following the first face, each additional vertex specifies a complete face
      _v[1] = _v[2];
      _v[2] = vertex_value;
      _faces(_idx_f)(0) = _v[0];
      _faces(_idx_f)(1) = _v[1];
      _faces(_idx_f)(2) = _v[2];
      _idx_f += 1; // next face
      _nFaces += 1;
    }

    break;
  }

  return 1;
}

static int face_normal_cb( p_ply_argument argument ) 
{
  long param_idx;

  // and whether this element value is x, y, or z
  ply_get_argument_user_data(argument, NULL, &param_idx);
  
  if (_idx_fn >= _nFaces) {
    std::cout << "ERROR: face count exceeded while loading face normals" << std::endl;
    return 0;
  }
  _face_normals(_idx_fn)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 2) _idx_fn += 1;

  return 1;
}

static int face_neighbor_cb(p_ply_argument argument) 
{
  long param_idx;

  // and whether this element value is vertex 1, 2, or 3
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_fnbr >= _nFaces) {
    std::cout << "ERROR: face count exceeded while loading face neighbors" << std::endl;
    return 0;
  }
  _face_neighbors(_idx_fnbr)(param_idx) = (int)ply_get_argument_value(argument);
  if (param_idx == 2) _idx_fnbr += 1;

  return 1;
}

static int point_3d_cb(p_ply_argument argument)
{
  long param_idx;

  // get whether this parameter value is x, y, or z
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_p >= _nPoints) {
    std::cout << "ERROR: point count exceeded" << std::endl;
    return 0;
  }
  _points_3d(_idx_p)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 2) _idx_p += 1;

  return 1;
}

static int point_normal_3d_cb(p_ply_argument argument)
{
  long param_idx;

  // get whether this parameter value is x, y, or z
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_pn >= _nPoints) {
    std::cout << "ERROR: point count exceeded" << std::endl;
    return 0;
  }
  _point_normals_3d(_idx_pn)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 2) _idx_pn += 1;

  return 1;
}

static int point_2d_cb(p_ply_argument argument)
{
  long param_idx;

  // get whether this parameter value is x or y
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_p >= _nPoints) {
    std::cout << "ERROR: point count exceeded" << std::endl;
    return 0;
  }
  _points_2d(_idx_p)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 1) _idx_p += 1;

  return 1;
}

static int point_normal_2d_cb(p_ply_argument argument)
{
  long param_idx;

  // get whether this parameter value is x or y
  ply_get_argument_user_data(argument, NULL, &param_idx);

  if (_idx_pn >= _nPoints) {
    std::cout << "ERROR: point count exceeded" << std::endl;
    return 0;
  }
  _point_normals_2d(_idx_pn)(param_idx) = ply_get_argument_value(argument);
  if (param_idx == 1) _idx_pn += 1;

  return 1;
}

int ply_io::read_ply_mesh(
  const std::string &input_ply,
  vctDynamicVector<vct3>    *vertices,
  vctDynamicVector<vctInt3> *faces,
  vctDynamicVector<vct3>    *face_normals,
  vctDynamicVector<vctInt3> *face_neighbors,
  vctDynamicVector<vct3>    *vertex_normals)
{
  long nVertices, nFaceLists, nFaceNormals, nFaceNeighbors, nVertexNormals;

  if (!vertices) {
    std::cout << "ERROR: vertices argument must be a valid pointer" << std::endl;
    return 0;
  }
  if (!faces) {
    std::cout << "ERROR: faces argument must be a valid pointer" << std::endl;
    return 0;
  }

  p_ply ply = ply_open(input_ply.c_str(), NULL, 0, NULL);
  if (!ply) {
    std::cout << "ERROR: failed to load ply file " << input_ply << std::endl;
    return 0;
  }

  if (!ply_read_header(ply)) {
    std::cout << "ERROR: failed to read header for ply file " << input_ply << std::endl;
    return 0;
  }
    
  nVertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 1);
  ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 2);

  nVertexNormals = ply_set_read_cb(ply, "vertex", "nx", vertex_normal_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "ny", vertex_normal_cb, NULL, 1);
  ply_set_read_cb(ply, "vertex", "nz", vertex_normal_cb, NULL, 2);

  nFaceLists = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0);
  
  nFaceNormals = ply_set_read_cb(ply, "face_normal", "nx", face_normal_cb, NULL, 0);
  ply_set_read_cb(ply, "face_normal", "ny", face_normal_cb, NULL, 1);
  ply_set_read_cb(ply, "face_normal", "nz", face_normal_cb, NULL, 2);

  nFaceNeighbors = ply_set_read_cb(ply, "face_neighbor", "nb1", face_neighbor_cb, NULL, 0);
  ply_set_read_cb(ply, "face_neighbor", "nb2", face_neighbor_cb, NULL, 1);
  ply_set_read_cb(ply, "face_neighbor", "nb3", face_neighbor_cb, NULL, 2);

  init_temp_variables_mesh(nVertices, nFaceLists, nFaceNormals, nFaceNeighbors, nVertexNormals);

  if (!ply_read(ply)) return 1;
  ply_close(ply);

  // faces array may be larger than the number of faces due to dynamic resizing
  if (_nFaces < _faces.size()) {
    _faces.resize(_nFaces);
  }

  // check that all values were loaded
  if (_idx_v != _vertices.size() ||
    _idx_f != _faces.size() ||
    (nFaceNormals > 0 && _idx_fn != _face_normals.size()) ||
    (nFaceNeighbors > 0 && _idx_fnbr != _face_neighbors.size()) ||
    (nVertexNormals > 0 && _idx_vn != _vertex_normals.size()) )
  {
    std::cout << "ERROR: PLY values did not load properly" << std::endl;
  }

  *vertices = _vertices;
  *faces = _faces;
  if (face_normals) *face_normals = _face_normals;
  if (face_neighbors) *face_neighbors = _face_neighbors;
  if (vertex_normals) *vertex_normals = _vertex_normals;

  return 1;
}

int ply_io::write_ply_mesh(const std::string &output_ply,
  const vctDynamicVector<vct3>    *vertices,
  const vctDynamicVector<vctInt3> *faces,
  const vctDynamicVector<vct3>    *face_normals,
  const vctDynamicVector<vctInt3> *face_neighbors,
  const vctDynamicVector<vct3>    *vertex_normals
  )
{
  e_ply_type type;
  e_ply_type length_type = PLY_UINT;
  e_ply_type list_type = PLY_UINT;

  if (!vertices) {
    std::cout << "ERROR: vertices argument must be a valid pointer" << std::endl;
    return 0;
  }
  if (!faces) {
    std::cout << "ERROR: faces argument must be a valid pointer" << std::endl;
    return 0;
  }

  // create PLY object
  p_ply oply = ply_create(output_ply.c_str(), PLY_ASCII, NULL, 0, NULL);
  if (!oply) {
    std::cout << "Unable to create file " << output_ply << std::endl;
    ply_close(oply);
    return 0;
  }

  std::string msg_element_error = "error adding element";
  std::string msg_property_error = "error adding property";
  std::string msg_comment_error = "Failed adding comments";

  // add elements and properties
  // vertices
  long nVertices;
  if (!vertices->empty()) {
    nVertices = vertices->size();
    if (vertex_normals && !vertex_normals->empty()) {
      if (vertices->size() != vertex_normals->size()) {
        std::cout << "ERROR: vertices and vertex normals must have same number of elements" << std::endl;
        return 0;
      }
    }
  }
  else if (vertex_normals) {
    nVertices = vertex_normals->size();
  }
  else {
    nVertices = 0;
  }
  if (nVertices > 0) {
    // add element
    ply_add_element(oply, "vertex", nVertices);

    // add element properties
    if (!vertices->empty()) {
      ply_add_property(oply, "x", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "y", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "z", PLY_FLOAT32, length_type, list_type);
    }
    if (!vertex_normals->empty()) {
      ply_add_property(oply, "nx", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "ny", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "nz", PLY_FLOAT32, length_type, list_type);
    }
  }

  // faces
  if (!faces->empty()) {
    // add element
    ply_add_element(oply, "face", faces->size());

    // add element properties
    ply_add_property(oply, "vertex_indices", PLY_LIST, length_type, list_type);
  }

  // face normals
  if (face_normals && !face_normals->empty()) {
    // add element
    ply_add_element(oply, "face_normal", face_normals->size());

    // add element properties
    if (!face_normals->empty()) {
      ply_add_property(oply, "nx", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "ny", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "nz", PLY_FLOAT32, length_type, list_type);
    }
  }

  // face neighbors
  if (face_neighbors && !face_neighbors->empty()) {
    // add element
    ply_add_element(oply, "face_neighbor", face_neighbors->size());

    // add element properties
    if (!face_normals->empty()) {
      ply_add_property(oply, "nb1", PLY_UINT, length_type, list_type);
      ply_add_property(oply, "nb2", PLY_UINT, length_type, list_type);
      ply_add_property(oply, "nb3", PLY_UINT, length_type, list_type);
    }
  }

  // add comments
  ply_add_comment(oply, "created by cisstICP");

  // write header
  ply_write_header(oply);

  // write data
  // vertices
  if (nVertices > 0) {
    for (long i = 0; i < nVertices; i++) {
      if (!vertices->empty()) {
        ply_write(oply, (*vertices)(i)[0]);
        ply_write(oply, (*vertices)(i)[1]);
        ply_write(oply, (*vertices)(i)[2]);
      }
      if (vertex_normals && !vertex_normals->empty()) {
        ply_write(oply, (*vertex_normals)(i)[0]);
        ply_write(oply, (*vertex_normals)(i)[1]);
        ply_write(oply, (*vertex_normals)(i)[2]);
      }
    }
  }

  // faces
  if (!faces->empty()) {
    long nFaces = faces->size();
    for (long i = 0; i < nFaces; i++) {
      ply_write(oply, 3.0); // list length
      ply_write(oply, (double)(*faces)(i)[0]);
      ply_write(oply, (double)(*faces)(i)[1]);
      ply_write(oply, (double)(*faces)(i)[2]);
    }
  }

  // face normals
  if (face_normals && !face_normals->empty()) {
    long nFaces = face_normals->size();
    for (long i = 0; i < nFaces; i++) {
      ply_write(oply, (double)(*face_normals)(i)[0]);
      ply_write(oply, (double)(*face_normals)(i)[1]);
      ply_write(oply, (double)(*face_normals)(i)[2]);
    }
  }

  // face neighbors
  if (face_neighbors && !face_neighbors->empty()) {
    long nFaces = face_neighbors->size();
    for (long i = 0; i < nFaces; i++) {
      ply_write(oply, (double)(*face_neighbors)(i)[0]);
      ply_write(oply, (double)(*face_neighbors)(i)[1]);
      ply_write(oply, (double)(*face_neighbors)(i)[2]);
    }
  }

  // close output file
  if (!ply_close(oply)) {
    std::cout << "ERROR: failed to close PLY output file" << std::endl;
    return 0;
  }
  return 1;
}

// read ply file containing a 3D point cloud
// returns  1: success  0: error
int ply_io::read_ply_pointcloud(const std::string &input_ply,
  vctDynamicVector<vct3> *points,
  vctDynamicVector<vct3> *point_normals
  )
{
  long nPoints, nPointNormals;

  p_ply ply = ply_open(input_ply.c_str(), NULL, 0, NULL);
  if (!ply) {
    std::cout << "ERROR: failed to load ply file " << input_ply << std::endl;
  }

  if (!ply_read_header(ply)) {
    std::cout << "ERROR: failed to read header for ply file " << input_ply << std::endl;
  }

  nPoints = ply_set_read_cb(ply, "vertex", "x", point_3d_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "y", point_3d_cb, NULL, 1);
  ply_set_read_cb(ply, "vertex", "z", point_3d_cb, NULL, 2);

  nPointNormals = ply_set_read_cb(ply, "vertex", "nx", point_normal_3d_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "ny", point_normal_3d_cb, NULL, 1);
  ply_set_read_cb(ply, "vertex", "nz", point_normal_3d_cb, NULL, 2);

  init_temp_variables_pointcloud3d(nPoints, nPointNormals);

  if (!ply_read(ply)) return 1;
  ply_close(ply);

  // check that all values were loaded
  if (_idx_p != _points_3d.size() || _idx_pn != _point_normals_3d.size())
  {
    std::cout << "ERROR: PLY values did not load properly" << std::endl;
  }

  if (points) *points = _points_3d;
  if (point_normals) *point_normals = _point_normals_3d;
  return 1;
}


// write 3D point cloud to ply file
// returns  1: success  0: error
int ply_io::write_ply_pointcloud(const std::string &output_ply,
  const vctDynamicVector<vct3> *points,
  const vctDynamicVector<vct3> *point_normals
  )
{
  e_ply_type type;
  e_ply_type length_type = PLY_UINT;
  e_ply_type list_type = PLY_UINT;

  // create PLY object
  p_ply oply = ply_create(output_ply.c_str(), PLY_ASCII, NULL, 0, NULL);
  if (!oply) {
    std::cout << "Unable to create file " << output_ply << std::endl;
    ply_close(oply);
    return 0;
  }

  std::string msg_element_error = "error adding element";
  std::string msg_property_error = "error adding property";
  std::string msg_comment_error = "Failed adding comments";

  // add elements and properties
  // points
  long nPoints;
  if (!points->empty()) {
    nPoints = points->size();
    if (point_normals && !point_normals->empty()) {
      if (points->size() != point_normals->size()) {
        std::cout << "ERROR: points and point normals must have same number of elements" << std::endl;
        return 0;
      }
    }
  }
  else if (point_normals) {
    nPoints = point_normals->size();
  }
  else {
    nPoints = 0;
  }

  if (nPoints > 0) {
    // add element
    ply_add_element(oply, "vertex", nPoints);

    // add element properties
    if (points && !points->empty()) {
      ply_add_property(oply, "x", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "y", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "z", PLY_FLOAT32, length_type, list_type);
    }
    if (point_normals && !point_normals->empty()) {
      ply_add_property(oply, "nx", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "ny", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "nz", PLY_FLOAT32, length_type, list_type);
    }
  }

  // add comments
  ply_add_comment(oply, "created by cisstICP");

  // write header
  ply_write_header(oply);

  // write data
  // points
  if (nPoints > 0) {
    for (long i = 0; i < nPoints; i++) {
      if (points && !points->empty()) {
        ply_write(oply, (*points)(i)[0]);
        ply_write(oply, (*points)(i)[1]);
        ply_write(oply, (*points)(i)[2]);
      }
      if (point_normals && !point_normals->empty()) {
        ply_write(oply, (*point_normals)(i)[0]);
        ply_write(oply, (*point_normals)(i)[1]);
        ply_write(oply, (*point_normals)(i)[2]);
      }
    }
  }

  // close output file
  if (!ply_close(oply)) {
    std::cout << "ERROR: failed to close PLY output file" << std::endl;
    return 0;
  }
  return 1;
}

// read ply file containing a 2D point cloud
// returns  1: success  0: error
int ply_io::read_ply_pointcloud(const std::string &input_ply,
  vctDynamicVector<vct2> *points,
  vctDynamicVector<vct2> *point_normals
  )
{
  long nPoints, nPointNormals;

  p_ply ply = ply_open(input_ply.c_str(), NULL, 0, NULL);
  if (!ply) {
    std::cout << "ERROR: failed to load ply file " << input_ply << std::endl;
  }

  if (!ply_read_header(ply)) {
    std::cout << "ERROR: failed to read header for ply file " << input_ply << std::endl;
  }

  nPoints = ply_set_read_cb(ply, "vertex", "x", point_2d_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "y", point_2d_cb, NULL, 1);

  nPointNormals = ply_set_read_cb(ply, "vertex", "nx", point_normal_2d_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "ny", point_normal_2d_cb, NULL, 1);

  init_temp_variables_pointcloud2d(nPoints, nPointNormals);

  if (!ply_read(ply)) return 1;
  ply_close(ply);

  // check that all values were loaded
  if (_idx_p != _points_2d.size() || _idx_pn != _point_normals_2d.size())
  {
    std::cout << "ERROR: PLY values did not load properly" << std::endl;
  }

  if (points) *points = _points_2d;
  if (point_normals) *point_normals = _point_normals_2d;
  return 1;
}

// write 2D point cloud to ply file
// returns  1: success  0: error
int ply_io::write_ply_pointcloud(const std::string &output_ply,
  const vctDynamicVector<vct2> *points,
  const vctDynamicVector<vct2> *point_normals
  )
{
  e_ply_type type;
  e_ply_type length_type = PLY_UINT;
  e_ply_type list_type = PLY_UINT;

  // create PLY object
  p_ply oply = ply_create(output_ply.c_str(), PLY_ASCII, NULL, 0, NULL);
  if (!oply) {
    std::cout << "Unable to create file " << output_ply << std::endl;
    ply_close(oply);
    return 0;
  }

  std::string msg_element_error = "error adding element";
  std::string msg_property_error = "error adding property";
  std::string msg_comment_error = "Failed adding comments";

  // add elements and properties
  // points
  long nPoints;
  if (!points->empty()) {
    nPoints = points->size();
    if (point_normals && !point_normals->empty()) {
      if (points->size() != point_normals->size()) {
        std::cout << "ERROR: points and point normals must have same number of elements" << std::endl;
        return 0;
      }
    }
  }
  else if (point_normals) {
    nPoints = point_normals->size();
  }
  else {
    nPoints = 0;
  }

  if (nPoints > 0) {
    // add element
    ply_add_element(oply, "vertex", nPoints);

    // add element properties
    if (points && !points->empty()) {
      ply_add_property(oply, "x", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "y", PLY_FLOAT32, length_type, list_type);
    }
    if (point_normals && !point_normals->empty()) {
      ply_add_property(oply, "nx", PLY_FLOAT32, length_type, list_type);
      ply_add_property(oply, "ny", PLY_FLOAT32, length_type, list_type);
    }
  }

  // add comments
  ply_add_comment(oply, "created by cisstICP");

  // write header
  ply_write_header(oply);

  // write data
  // points
  if (nPoints > 0) {
    for (long i = 0; i < nPoints; i++) {
      if (points && !points->empty()) {
        ply_write(oply, (*points)(i)[0]);
        ply_write(oply, (*points)(i)[1]);
      }
      if (point_normals && !point_normals->empty()) {
        ply_write(oply, (*point_normals)(i)[0]);
        ply_write(oply, (*point_normals)(i)[1]);
      }
    }
  }

  // close output file
  if (!ply_close(oply)) {
    std::cout << "ERROR: failed to close PLY output file" << std::endl;
    return 0;
  }
  return 1;
}