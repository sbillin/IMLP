classdef mexInterface_Alg2D_DirPDTree_CP_PointCloud < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_Alg2D_DirPDTree_CP_PointCloud()
      this.h = mexAlg2D_DirPDTree_CP_PointCloud('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlg2D_DirPDTree_CP_PointCloud('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %  points    ~ 2D model shape points   (N x 2)
    function Initialize( this, points )
                       
      mexAlg2D_DirPDTree_CP_PointCloud('Initialize',this.h, ...
        points');
    end
        
    %% Compute Matches
    % inputs:
    %  samplePts     ~ (Ns x 2)  2D sample points        
    %  sampleNorms   ~ (Ns x 2)  2D sample orientations  
    %  matchDatumsInit  ~ (Ns x 1)  initializer for match datums [optional]
    % outputs:    
    %  matchPts2D    ~ (Ns x 2)
    %  matchNorms2D  ~ (Ns x 2) 
    %  matchDatums   ~ (Ns x 1)  match datums
    function [matchPts2D, matchNorms2D, matchDatums] = ...
      ComputeMatches(this, samplePts, sampleNorms, matchDatumsInit)
     
      if exist('matchDatumsInit','var') && ~isempty(matchDatumsInit)
        [matchPts2D, matchNorms2D, matchDatums] = ...
          mexAlg2D_DirPDTree_CP_PointCloud('ComputeMatches', this.h, ...
          samplePts', sampleNorms', matchDatumsInit);
      else
        [matchPts2D, matchNorms2D, matchDatums] = ...
          mexAlg2D_DirPDTree_CP_PointCloud('ComputeMatches', this.h, ...
          samplePts', sampleNorms');
      end

      matchPts2D = matchPts2D';
      matchNorms2D = matchNorms2D';
    end    
  end
end