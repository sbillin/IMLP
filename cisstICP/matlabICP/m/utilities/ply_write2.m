function ply_write2( data, path, format )
% This function writes mesh / point cloud data stored in a struct 
% to a PLY file with the following optional fields:
%  
%  struct field name           PLY field.property name
% -------------------         -------------------------
%   vertices                   vertex.{x,y,z}
%   vertex_normals             vertex.{nx,ny,nz}
%   faces                      face.vertex_indices [list]
%   face_normals               face_normal.{nx,ny,nz}
%   face_neighbors             face_neighbor.{nb1,nb2,nbe3}
%     
%  This function is dependent on the PLY_IO library, and assumes that
%  the "ply_write" function from this library is on the path;
%  some sections of code below have also been extracted and modified
%  from the "tri_mesh_to_ply.m" function from this library.
%  
%    PLY_IO Library:  http://people.sc.fsu.edu/~jburkardt/m_src/ply_io/ply_io.html
%
% inputs
%
%  data         struct containing the mesh / point cloud data
%
%  path         path to outpu PLY file
%
%  format
%      'ascii'                  ASCII text data (default)
%      'binary_little_endian'   binary data, little endian
%      'binary_big_endian'      binary data, big endian
%

if ~exist('format','var') || isempty(format)
  format = 'ascii';
end


% transfer the data to a PLY-formatted data structure
ply_data = [];

%% vertex positions
if isfield(data,'vertices')
  ply_data.vertex.x = data.vertices(:,1);
  ply_data.vertex.y = data.vertices(:,2);
  ply_data.vertex.z = data.vertices(:,3);
end

%% vertex normals
if isfield(data,'vertex_normals')
  ply_data.vertex.nx = data.vertex_normals(:,1);
  ply_data.vertex.ny = data.vertex_normals(:,2);
  ply_data.vertex.nz = data.vertex_normals(:,3);
end

%% faces
if isfield(data,'faces')
  numFaces = size(data.faces,1);
  ply_data.face.vertex_indices = cell ( numFaces, 1 );
  for k = 1 : numFaces
    % subtract one to convert from 1-base to 0-base indexing
    ply_data.face.vertex_indices{k} = data.faces(k,:) - 1;
  end
end

%% face normals
if isfield(data,'face_normals')
  ply_data.face_normal.nx = data.face_normals(:,1);
  ply_data.face_normal.ny = data.face_normals(:,2);
  ply_data.face_normal.nz = data.face_normals(:,3);  
end

%% face neighbors
if isfield(data,'face_neighbors')
  % subtract one to convert from 1-base to 0-base indexing
  ply_data.face_neighbor.nb1 = data.face_neighbors(:,1) - 1;
  ply_data.face_neighbor.nb2 = data.face_neighbors(:,2) - 1;
  ply_data.face_neighbor.nb3 = data.face_neighbors(:,3) - 1;  
end


% write to PLY file
ply_write( ply_data, path, format )

end
