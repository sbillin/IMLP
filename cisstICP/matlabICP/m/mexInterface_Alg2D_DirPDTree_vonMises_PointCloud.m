classdef mexInterface_Alg2D_DirPDTree_vonMises_PointCloud < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods

    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_Alg2D_DirPDTree_vonMises_PointCloud()
      this.h = mexAlg2D_DirPDTree_vonMises_PointCloud('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlg2D_DirPDTree_vonMises_PointCloud('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %  points       ~ 2D model shape points                (N x 2)
    %  norms        ~ 2D model shape point orientations    (N x 2)
    function Initialize( this, points, norms )
      mexAlg2D_DirPDTree_vonMises_PointCloud('Initialize',this.h, ...
        points', norms');
    end
        
    %% Compute Matches
    % inputs:
    %  samplePts        ~ (Ns x 2)  2D sample points        
    %  sampleNorms      ~ (Ns x 2)  2D sample orientations  
    %  sigma2           ~ (Ns x 1)  2D sample point variance
    %  k                ~ (Ns x 1)  2D sample orientation concentration    
    %  matchDatumsInit  ~ (Ns x 1)  initializer for match datums          [optional]
    %  match_ThetaMax   ~ (Ns x 1)  maximum orientation error for a match (degrees) [optional]
    %  
    % outputs:    
    %  matchPts2D    ~ (Ns x 2)  (in image coords)
    %  matchNorms2D  ~ (Ns x 2)  (in image coords)     
    %  matchDatums   ~ (Ns x 1)  match datums
    %  matchPermitted ~ (Ns x 1)  permitted matches  (i.e. true if theta < thetaMax)
    function [matchPts2D, matchNorms2D, matchDatums, matchPermitted] = ...
      ComputeMatches(this, samplePts, sampleNorms, sigma2, k, matchDatumsInit, match_ThetaMax)

      if ~exist('matchDatumsInit','var')
        matchDatumsInit = [];
      end
    
      if ~exist('match_ThetaMax', 'var') || isempty(match_ThetaMax)
        match_ThetaMax = 180;  % default
      end
      match_ThetaMax = match_ThetaMax * pi/180;    
    
      [matchPts2D, matchNorms2D, matchDatums, matchPermitted] = ...
        mexAlg2D_DirPDTree_vonMises_PointCloud('ComputeMatches', this.h, ...
        samplePts', sampleNorms', sigma2, k, matchDatumsInit, match_ThetaMax);

      matchPts2D = matchPts2D';
      matchNorms2D = matchNorms2D';
      matchPermitted = matchPermitted';
      matchDatums = bsxfun(@plus,matchDatums,1);  % add 1 to convert base-0 (C++) to base-1 (Matlab) indexing
    end

  end
end