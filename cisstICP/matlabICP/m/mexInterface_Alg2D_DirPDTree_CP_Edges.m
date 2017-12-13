classdef mexInterface_Alg2D_DirPDTree_CP_Edges < handle
% MATLAB class wrapper for C++ mex library

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_Alg2D_DirPDTree_CP_Edges()
      this.h = mexAlg2D_PDTree_CP_Edges('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlg2D_PDTree_CP_Edges('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %  edgesV1       ~ first vertices of 2D edges    (Ne x 2)  (in image coords)
    %  edgesV2       ~ second vertices of 2D edges   (Ne x 2)  (in image coords)
    %  edgesNorm     ~ 2D edge normals               (Ne x 2)  (in image coords)
    %  edgesV1_3D    ~ first vertices of 3D edges    (Ne x 3)  (in camera coords)
    %  edgesV2_3D    ~ second vertices of 3D edges   (Ne x 3)  (in camera coords)
    function Initialize( this, edgesV1, edgesV2, edgesNorm, ...
                         edgesV1_3D, edgesV2_3D)
                       
      mexAlg2D_PDTree_CP_Edges('Initialize',this.h, ...
        edgesV1', edgesV2', edgesNorm', ...
        edgesV1_3D', edgesV2_3D');
    end
        
    %% Compute Matches
    % inputs:
    %  samplePts     ~ (Ns x 2)  2D sample points        
    %  sampleNorms   ~ (Ns x 2)  2D sample orientations  
    %  matchDatumsInit  ~ (Ns x 1)  initializer for match datums [optional]
    % outputs:    
    %  matchPts2D    ~ (Ns x 2)  (in image coords)
    %  matchNorms2D  ~ (Ns x 2)  (in image coords)     
    %  matchPts3D    ~ (Ns x 3)  (in camera coords)
    %  matchDatums   ~ (Ns x 1)  match datums
    function [matchPts2D, matchNorms2D, matchPts3D, matchDatums] = ...
      ComputeMatches(this, samplePts, sampleNorms, matchDatumsInit)
     
      if exist('matchDatumsInit','var') && ~isempty(matchDatumsInit)
        [matchPts2D, matchNorms2D, matchPts3D, matchDatums] = ...
          mexAlg2D_PDTree_CP_Edges('ComputeMatches', this.h, ...
          samplePts', sampleNorms', matchDatumsInit);
      else
        [matchPts2D, matchNorms2D, matchPts3D, matchDatums] = ...
          mexAlg2D_PDTree_CP_Edges('ComputeMatches', this.h, ...
          samplePts', sampleNorms');
      end

      matchPts2D = matchPts2D';
      matchNorms2D = matchNorms2D';
      matchPts3D = matchPts3D';
    end    
  end
end