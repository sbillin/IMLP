classdef mexInterface_Alg2D_DirPDTree_vonMises_Edges < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods

    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_Alg2D_DirPDTree_vonMises_Edges()
      this.h = mexAlg2D_DirPDTree_vonMises_Edges('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlg2D_DirPDTree_vonMises_Edges('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %  edgesV1       ~ first vertices of 2D edges    (Ne x 2)  (in image coords)
    %  edgesV2       ~ second vertices of 2D edges   (Ne x 2)  (in image coords)
    %  edgesNorm     ~ 2D edge normals               (Ne x 2)  (in image coords)
    function Initialize( this, edgesV1, edgesV2, edgesNorm)
      mexAlg2D_DirPDTree_vonMises_Edges('Initialize',this.h, ...
        edgesV1', edgesV2', edgesNorm');
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
        mexAlg2D_DirPDTree_vonMises_Edges('ComputeMatches', this.h, ...
        samplePts', sampleNorms', sigma2, k, matchDatumsInit, match_ThetaMax);

      matchPts2D = matchPts2D';
      matchNorms2D = matchNorms2D';
      matchPermitted = matchPermitted';
      matchDatums = bsxfun(@plus,matchDatums,1);  % add 1 to convert base-0 (C++) to base-1 (Matlab) indexing
    end
    
    % %% Set Noise Model
    % % inputs:
    % %   k       ~ (Ns x 1)  or  scalar
    % %   sigma2  ~ (Ns x 1)  or  scalar
    % function SetNoiseModel( k, sigma2 )
    %   mexAlg2D_DirPDTree_vonMises_Edges('SetNoiseModel', this.h, k, sigma2);
    % end
  end
end