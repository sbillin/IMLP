function [ data ] = ply_read2 ( Path )
% This function loads a ply file to a struct with the following optional fields:
%  
%  struct field name           PLY field.property name
% -------------------         -------------------------
%   vertices                   {'vertex','Vertex','point','Point','pts','Pts'}.{x,y,z}
%   vertex_normals             {'vertex','Vertex','point','Point','pts','Pts'}.{nx,ny,nz}
%   faces                      {'face','Face','poly','Poly','tri','Tri'}.{'vertex_indices','vertex_indexes','vertex_index','indices','indexes'}
%   face_normals               face_normal.{nx,ny,nz}
%   face_neighbors             face_neighbor.{nb1,nb2,nbe3}
%     
%  This function is dependent on the PLY_IO library, and assumes that
%  the "ply_read" function from this library is on the path;
%  some sections of code below have also been extracted and modified
%  from the "ply_read" function.
%  
%    PLY_IO Library:  http://people.sc.fsu.edu/~jburkardt/m_src/ply_io/ply_io.html
%

% load the ply file
[Elements, Comments] = ply_read( Path );


%% process the returned PLY data and return the processed data in a struct
data = [];

%% vertex field
Names = {'vertex','Vertex','point','Point','pts','Pts'};
ind_arr = isfield(Elements,Names);
if any(ind_arr)
  if sum(ind_arr) > 1
    error('multiple fields found for vertex (or point) type')
  end
  
  fieldName = Names{ind_arr};
  field = Elements.(fieldName);
  prop_names = fieldnames(field);

  % vertices
  if ( any(strcmp(prop_names,'x')) & any(strcmp(prop_names,'y')) & any(strcmp(prop_names,'z')) )
    data.vertices = [field.x, field.y, field.z];
  else
    warning('no vertex position property found for vertex field')
  end
  
  % vertex normals
  if ( any(strcmp(prop_names,'nx')) & any(strcmp(prop_names,'ny')) & any(strcmp(prop_names,'nz')) )
    data.vertex_normals = [field.nx, field.ny, field.nz];
  end
end


%% face field
Names = {'face','Face','poly','Poly','tri','Tri'};
ind_arr = isfield(Elements,Names);
if any(ind_arr)
  if sum(ind_arr) > 1
    error('multiple fields found for triangle face type')
  end
  
  fieldName = Names{ind_arr};
  field = Elements.(fieldName);

  Names = {'vertex_indices','vertex_indexes','vertex_index','indices','indexes'};
  ind_arr = isfield(field,Names);
  if any(ind_arr)
    if sum(ind_arr) > 1
      error('multiple vertex index property types found for face field')
    end
  else
    error('no vertex index property found for face field');
  end
  
  propertyName = Names{ind_arr};
  vertex_indices = field.(propertyName);
  
  % vertex indices
  % convert vertex index lists to triangle connectivity
  N = length(vertex_indices);
  indices = zeros(N*2,3);  % add extra length in case of multiple faces per row
  Extra = 0;
  for k = 1 : N
    indices(k,:) = vertex_indices{k}(1:3);
    
    for j = 4 : length(vertex_indices{k})
      Extra = Extra + 1;
      indices(N + Extra,1) = vertex_indices{k}(1);
      indices(N + Extra,2) = vertex_indices{k}(j-1);
      indices(N + Extra,3) = vertex_indices{k}(j);
    end

  end  
  % Add 1 to convert from zero-based (PLY) to one-based (MATLAB) indexing
  indices = indices(1:N+Extra,:) + 1;       
  data.faces = indices;
end

%% face normals field
Names = {'face_normal'};
ind_arr = isfield(Elements,Names);
if any(ind_arr)
  if sum(ind_arr) > 1
    error('multiple fields found for face normal type')
  end
  
  fieldName = Names{ind_arr};
  field = Elements.(fieldName);
  prop_names = fieldnames(field);

  % face normals
  if ( any(strcmp(prop_names,'nx')) & any(strcmp(prop_names,'ny')) & any(strcmp(prop_names,'nz')) )
    data.face_normals = [field.nx, field.ny, field.nz];
  else
    error('properties for face normals field are unrecognized')
  end
end

%% face neighbors field
Names = {'face_neighbor'};
ind_arr = isfield(Elements,Names);
if any(ind_arr)
  if sum(ind_arr) > 1
    error('multiple fields found for face neighbor type')
  end
  
  fieldName = Names{ind_arr};
  field = Elements.(fieldName);
  prop_names = fieldnames(field);

  % face neighbors
  if ( any(strcmp(prop_names,'nb1')) & any(strcmp(prop_names,'nb2')) & any(strcmp(prop_names,'nb3')) )
    % add one for base-1 indexing
    data.face_neighbors = [field.nb1, field.nb2, field.nb3] + 1;
  else
    error('properties for face neighbors field are unrecognized')
  end
end


end
