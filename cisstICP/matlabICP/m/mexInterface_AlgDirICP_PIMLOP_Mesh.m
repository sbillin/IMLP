classdef mexInterface_AlgDirICP_PIMLOP_Mesh < handle
% MATLAB class wrapper to an underlying C++ class

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexInterface_AlgDirICP_PIMLOP_Mesh()
      this.h = mexAlgDirICP_PIMLOP_Mesh('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlgDirICP_PIMLOP_Mesh('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %   V     ~ target vertices   Nv x 3
    %   T,Tn  ~ target triangles and triangle normals   Nt x 3
    %   Xp    ~ source points             Ns x 3
    %   Xpln  ~ source 2d normals         Ns x 2
    %   Rx_plane ~ 3d source normal coords: Xn3d = Rx_plane * [Xn2d;0]  3 x 3 x Ns
    %   k ~ concentration of source orientations    Ns x 1
    %   M ~ covariance matrices of soure points     3 x 3 x Ns
    function Initialize(this, V,T,Tn, Xp,Xpln,Rx_plane,k,M)
      % convert triangles to C++ indexing (base-0) and to integer data type
      Tcpp = int32(T - ones(size(T)));
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  mesh
      % //    V ~ 3D vertex positions     (3 x Nv double)
      % //    T ~ triangle vertex indices (3 x Nt integer)
      % //    N ~ 3D triangle normals     (3 x Nt double)
      % //  sample points                 (3 x Ns double)
      % //  sample 2d normals             (2 x Ns double)
      % //  Rx_plane ~ 3d sample normal coords: Xn3d = Rx_plane * [Xn2d;0]  (3 x 3 x Ns double)
      % //  k ~ concentration of sample orientations    (Ns x 1 double)
      % //  M ~ covariance matrices of sample points    (3 x 3 x Ns double)
      mexAlgDirICP_PIMLOP_Mesh('Initialize',this.h,...
        V',Tcpp',Tn',...
        Xp',Xpln',Rx_plane, k,M);  
    end
    
    %% Set Samples
    % inputs:
    %   source points             Ns x 3
    %   source 2d normals         Ns x 2
    %   Rx_plane ~ 3d source normal coords: Xn3d = Rx_plane * [Xn2d;0]  3 x 3 x Ns
    %   k ~ concentration of source orientations    Ns x 1
    %   M ~ covariance matrices of soure points     3 x 3 x Ns
    function SetSamples(this, Xp,Xpln,Rx_plane,k,M)

      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  sample points                 (3 x Ns double)
      mexAlgDirICP_PIMLOP_Mesh('SetSamples',this.h,...
        Xp',Xpln',Rx_plane, k,M);
    end
        
    %% ICP: Initialize Parameters
    % inputs:
    %   [R0,t0]   ~ initial registration guess   (3 x 4 double)
    function ICP_InitializeParameters(this, F0)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  Freg ~ [R,t]      (3 x 4 double)
      mexAlgDirICP_PIMLOP_Mesh('ICP_InitializeParameters',this.h, F0(1:3,:));
    end    
        
    %% ICP: Compute Matches
    % computes matches and post-match parameters
    %
    % outputs:
    %   Yp   ~ match positions    N x 3
    %   Yn   ~ match orientations N x 3
    function [Yp, Yn] = ICP_ComputeMatches(this)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  SampleMatchPts ~ matches for 3D sample pts    (3 x numPts double)  [optional]      
      [Yp, Yn] = mexAlgDirICP_PIMLOP_Mesh('ICP_ComputeMatches', this.h);
      Yp = Yp';
      Yn = Yn';
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
      Freg = mexAlgDirICP_PIMLOP_Mesh('ICP_RegisterMatches', this.h);
      F = getFrm3(Freg(1:3,1:3),Freg(1:3,4));
    end        

    function fval = ICP_EvaluateErrorFunction( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  fval  ~ value of error function (double)
      % //      
      fval = mexAlgDirICP_PIMLOP_Mesh('ICP_EvaluateErrorFunction', this.h);
    end
    
    function bTerminate = ICP_Terminate( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  bTerminate  ~ boolean specifying algorithm-specific termination (logical scalar)
      % //
      bTerminate = mexAlgDirICP_PIMLOP_Mesh('ICP_Terminate', this.h);
    end
  end
end