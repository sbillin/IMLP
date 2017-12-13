%clear, clc

addpath('..\..\..\dependencies\Matlab Dependencies\PLY_IO');

inputDir = '..\..\..\test_data';
plyFile = 'ProximalFemur.ply';
plyPath = [inputDir,'\',plyFile];

bComputeVertexNormals = false;
bComputeFaceNormals = true;
bComputeFaceNeighbors = true;

% load standard PLY file
mesh = ply_read2( plyPath );

if bComputeFaceNormals
  mesh.face_normals = MeshComputeTriangleNormals( mesh.vertices, mesh.faces );
end

if bComputeFaceNeighbors
  % compute neighbors
  [mesh.faces, mesh.face_normals, mesh.face_neighbors] = ...
    MeshComputeNeighbors( mesh.faces, mesh.face_normals );
end

if bComputeVertexNormals
  % compute vertex normals
  [mesh.vertex_normals] = ...
    MeshComputeVertexNormals( mesh.vertices, mesh.faces, mesh.face_normals );
end

% save back to PLY format
[dir, name, ext] = fileparts( plyPath );
newPLYFile = [inputDir,'/',name,'_new.ply'];
ply_write2( mesh, newPLYFile );
