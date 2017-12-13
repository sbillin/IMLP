classdef mexInterface_AlgPDTree_MLP_Mesh < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_AlgPDTree_MLP_Mesh()
      this.h = mexAlgPDTree_MLP_Mesh('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlgPDTree_MLP_Mesh('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %   V     ~ mesh vertices   Nv x 3
    %   T,Tn  ~ mesh triangles and triangle normals   Nt x 3 
    %   Tcov  ~ mesh triangle covariances 3 x 3 x Nt  [optional] 
    function Initialize( this, V,T,Tn,Tcov)
      % convert triangles to C++ indexing (base-0) and to integer data type
      Tcpp = int32(T - ones(size(T)));
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  mesh
      % //    V ~ 3D vertex positions     (3 x Nv double)
      % //    T ~ triangle vertex indices (3 x Nt integer)
      % //    N ~ triangle normals        (3 x Nt double)
      % //  mesh noise model (optional)
      % //    Tcov ~ triangle covariances (3 x 3 x Nt double)      
      % //
      if ~exist('Tcov','var') || isempty(Tcov)
        mexAlgPDTree_MLP_Mesh('Initialize', this.h,...
          V',Tcpp',Tn');
      else
        mexAlgPDTree_MLP_Mesh('Initialize', this.h,...
          V',Tcpp',Tn',Tcov);        
      end
    end

    %% Compute Matches
    % inputs:
    %  samplePts        ~ (Ns x 3)      sample points        
    %  sampleCov        ~ (3 x 3 x Ns)  sample covariances
    %  matchDatumsInit  ~ (Ns x 1)      initializer for match datums [optional]
    % outputs:
    %  matchPts      ~ (Ns x 3)  match points
    %  matchNorms    ~ (Ns x 3)  match norms
    %  matchDatums   ~ (Ns x 1)  match datums
    function [matchPts, matchNorms, matchDatums] = ...
      ComputeMatches(this, samplePts, sampleCov, matchDatumsInit)
     
      if exist('matchDatumsInit','var') && ~isempty(matchDatumsInit)
        [matchPts, matchNorms, matchDatums] = ...
          mexAlgPDTree_MLP_Mesh('ComputeMatches', this.h, samplePts', sampleCov, matchDatumsInit);
      else
        [matchPts, matchNorms, matchDatums] = ...
          mexAlgPDTree_MLP_Mesh('ComputeMatches', this.h, samplePts', sampleCov);
      end
      
      matchPts = matchPts';
      matchNorms = matchNorms';
    end    
  end
end