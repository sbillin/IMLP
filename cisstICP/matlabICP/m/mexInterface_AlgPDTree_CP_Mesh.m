classdef mexInterface_AlgPDTree_CP_Mesh < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_AlgPDTree_CP_Mesh()
      this.h = mexAlgPDTree_CP_Mesh('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlgPDTree_CP_Mesh('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %   V     ~ target vertices   Nv x 3
    %   T,Tn  ~ target triangles and triangle normals   Nt x 3 
    function Initialize( this, V,T,Tn)
      % convert triangles to C++ indexing (base-0) and to integer data type
      Tcpp = int32(T - ones(size(T)));
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  mesh
      % //    V ~ 3D vertex positions     (3 x Nv double)
      % //    T ~ triangle vertex indices (3 x Nt integer)
      % //    N ~ 3D triangle normals     (3 x Nt double)
      % //
      mexAlgPDTree_CP_Mesh('Initialize', this.h,...
        V',Tcpp',Tn');
    end

    %% Compute Matches
    % inputs:
    %  samplePts        ~ (Ns x 3)  sample points        
    %  matchDatumsInit  ~ (Ns x 1)  initializer for match datums [optional]
    % outputs:
    %  matchPts      ~ (Ns x 3)  match points
    %  matchNorms    ~ (Ns x 3)  match norms
    %  matchDatums   ~ (Ns x 1)  match datums
    function [matchPts, matchNorms, matchDatums] = ...
      ComputeMatches(this, samplePts, matchDatumsInit)
     
      if exist('matchDatumsInit','var') && ~isempty(matchDatumsInit)
        [matchPts, matchNorms, matchDatums] = ...
          mexAlgPDTree_CP_Mesh('ComputeMatches', this.h, ...
          samplePts', matchDatumsInit);
      else
        [matchPts, matchNorms, matchDatums] = ...
          mexAlgPDTree_CP_Mesh('ComputeMatches', this.h, samplePts');      
      end
      
      matchPts = matchPts';
      matchNorms = matchNorms';
    end    
  end
end