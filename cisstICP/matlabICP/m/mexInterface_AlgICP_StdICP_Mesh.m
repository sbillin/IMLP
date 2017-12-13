classdef mexInterface_AlgICP_StdICP_Mesh < handle
% MATLAB class wrapper to an underlying C++ class

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_AlgICP_StdICP_Mesh()
      this.h = mexAlgICP_StdICP_Mesh('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlgICP_StdICP_Mesh('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %   X     ~ source points     N x 3
    %   V     ~ target vertices   Nv x 3
    %   T,Tn  ~ target triangles and triangle normals   Nt x 3
    function Initialize(this, X, V,T)
      % convert triangles to C++ indexing (base-0) and to integer data type
      Tcpp = int32(T - ones(size(T)));
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  sample points                 (3 x Ns double)
      % //  mesh
      % //    V ~ 3D vertex positions     (3 x Nv double)
      % //    T ~ triangle vertex indices (3 x Nt integer)
      mexAlgICP_StdICP_Mesh('Initialize',this.h, X', V',Tcpp');
    end
    
    %% Set Samples
    % inputs:
    %   X     ~ source points     N x 3
    function SetSamples(this, X)

      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  sample points                 (3 x Ns double)
      mexAlgICP_StdICP_Mesh('SetSamples',this.h, X');  
    end
        
    %% ICP: Initialize Parameters
    % inputs:
    %   [R0,t0]   ~ initial registration guess   (3 x 4 double)
    function ICP_InitializeParameters(this, F0)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  Freg ~ [R,t]      (3 x 4 double)
      mexAlgICP_StdICP_Mesh('ICP_InitializeParameters',this.h, F0(1:3,:));
    end    
        
    %% ICP: Compute Matches
    % computes matches and post-match parameters
    %
    % outputs:
    %   Y   ~ matches (inliers & outliers)   N x 3
    function Y = ICP_ComputeMatches(this)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  SampleMatchPts ~ matches for 3D sample pts    (3 x numPts double)  [optional]      
      Y = mexAlgICP_StdICP_Mesh('ICP_ComputeMatches', this.h);
      Y = Y';
    end    
    
    %% ICP: Register Matches
    % registers matches and computes post-registration parameters
    %
    % outputs:
    %   F   ~ registration parameters
    function [F] = ICP_RegisterMatches(this)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  Freg ~ [R,t]      (3 x 4 double)    (Y = Freg * X)
      % //
      Freg = mexAlgICP_StdICP_Mesh('ICP_RegisterMatches', this.h);
      F = getFrm3(Freg(1:3,1:3),Freg(1:3,4));
    end        

    %% ICP: Update Parameters (Post Register)
    %   use this function when explicitly computing parameters following 
    %   registration, (i.e. for example if computing the registration in
    %   Matlab rather than using the mex routine for this)
    % inputs:
    %   [R0,t0]   ~ initial registration guess   (3 x 4 double)    
    function ICP_UpdateParameters_PostRegister( this, Freg )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  Freg ~ [R,t]      (3 x 4 double)
      mexAlgICP_StdICP_Mesh('ICP_UpdateParameters_PostRegister',this.h, Freg(1:3,:));      
    end
    
    function fval = ICP_EvaluateErrorFunction( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  fval  ~ value of error function (double)
      % //      
      fval = mexAlgICP_StdICP_Mesh('ICP_EvaluateErrorFunction', this.h);
    end
    
    function bTerminate = ICP_Terminate( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  bTerminate  ~ boolean specifying algorithm-specific termination (logical scalar)
      % //
      bTerminate = mexAlgICP_StdICP_Mesh('ICP_Terminate', this.h);
    end
  end
end