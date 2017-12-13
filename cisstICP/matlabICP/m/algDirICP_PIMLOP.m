classdef algDirICP_PIMLOP < algICP

properties

  bEnableDebug;
  bTargetAsMesh;
  
  % sample data
  Xp;          % source pts        N x 3
  %Xp_xfm;      % xfmd source pts   N x 3
  Yp;          % match pts         N x 3
  
  Xpln;        % source normals 2d (in-plane)   N x 2  
  Yn;          % match normals 3d  N x 3
  Rx_plane;    % rotation to align in-plane normals with 3d coordinates of source data  3 x 3 x N
  
  k;    % concentration of source orientations    N x 1
  M;    % covariance of source positions          3 x 3 x N
  
  % target shape
  % mesh and point cloud target:
  V;        % vertices            Nv x 3
  % mesh target only:
  T;        % triangle faces (vertex indices)   Nt x 3
  Tn;       % triangle normals                  Nt x 3
  
  % member objects
  mexAlg_PIMLOP;   % mex object implementing this algorithm
  
end

methods
  
  %% Constructor
  function [this] = algDirICP_PIMLOP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlg_PIMLOP.delete();
  end  
  
  %% Initialize Algorithm
  % V     ~ target vertices   Nv x 3
  % T,Tn  ~ target triangles and triangle normals   Nt x 3
  % Xp    ~ source points             Ns x 3
  % Xpln  ~ source 2d normals         Ns x 2
  % Rx_plane  ~ 3d source normal coords: Xn3d = Rx_plane * [Xn2d;0]  3 x 3 x Ns
  % k   ~ concentration of source orientations    Ns x 1
  % M   ~ covariance matrices of soure points     3 x 3 x Ns
  % bTargetAsMesh   ~ set point cloud or mesh target
  % bEnableDebug    ~ enable debug messages
  function Initialize( this, ...
    V,T,Tn, ...
    Xp,Xpln,Rx_plane,k,M, ...    
    bTargetAsMesh,bEnableDebug )
    
    % initialize superclass
    Initialize@algICP( this );
    
    % debug enable
    if ~exist('bEnableDebug','var') || isempty(bEnableDebug)
      bEnableDebug = 0;
    end
    this.bEnableDebug = bEnableDebug;    
    
    % source shape
    if ~exist('Xp','var') || isempty(Xp)
      warning('Missing input Xp')
      this.Xp = [];
    else      
      if size(Xp,2) ~= 3
        error('Invalid size of Xp')
      end
      this.Xp = Xp;
    end
    if ~exist('Xpln','var') || isempty(Xpln)
      warning('Missing input Xpln')
      this.Xpln = [];
    else
      this.Xpln = Xpln;
    end
    if ~exist('Rx_plane','var') || isempty(Rx_plane)
      warning('Missing input Rx_plane')
      this.Rx_plane = [];
    else
      this.Rx_plane = Rx_plane;
    end
    if ~exist('k','var') || isempty(k)
      warning('Missing input k')
      this.k = [];
    end
    if ~exist('M','var') || isempty(M)
      warning('Missing input M')
      this.M = [];
    else
      this.M = M;
    end    
        
    % target shape
    if ~exist('bEnableDebug','var') || isempty(bTargetAsMesh)
      bTargetAsMesh = 1;
    end
    if size(V,2) ~= 3
      error('Invalid size of target vertices (V)')
    end
    this.bTargetAsMesh = bTargetAsMesh;
    this.V = V;
    if (bTargetAsMesh)
      this.T = T;
      this.Tn = Tn;
    else
      this.T = [];
      this.Tn = [];
    end
        
    % Mex objects
    if (bTargetAsMesh)
      utlDebugMsg( this.bEnableDebug,'Creating mex object: mexInterface_AlgICP_vonMisesPrj_Mesh...\n')
      this.mexAlg_PIMLOP = mexInterface_AlgDirICP_PIMLOP_Mesh();
      utlDebugMsg( this.bEnableDebug,' ...done\n');
      
      utlDebugMsg( this.bEnableDebug,'Initializing mexAlg_PIMLOP...\n')
      this.mexAlg_PIMLOP.Initialize( ...
        this.V, this.T, this.Tn,...
        this.Xp, this.Xpln, this.Rx_plane, this.k, this.M);
      utlDebugMsg( this.bEnableDebug,' ...done\n');
    else
      error('Point Cloud version of mexAlgICP_vonMisesPrj not ready');
    end    
  end  
  
  %% Set Samples
  % Xp    ~ source points             Ns x 3
  % Xpln  ~ source 2d normals         Ns x 2
  % Rx_plane  ~ 3d source normal coords: Xn3d = Rx_plane * [Xn2d;0]  3 x 3 x Ns
  % k   ~ concentration of source orientations    Ns x 1
  % M   ~ covariance matrices of soure points     3 x 3 x Ns
  function SetSamples( this, Xp,Xpln,Rx_plane,k,M )  
    this.Xp = Xp;
    this.Xpln = Xpln;
    this.Rx_plane = Rx_plane;
    this.k = k;
    this.M = M;
    this.mexAlg_PIMLOP.SetSamples( Xp,Xpln,Rx_plane,k,M );
  end
  
  %% ICP: Initialize Registration
  function ICP_InitializeParameters( this, F0, algPayload )
    
    % initial guess
    if ~exist('F0','var') || isempty(F0)
      F0 = getFrm3(eye(3),zeros(3,1));
    end
    
    % initialize superclass
    ICP_InitializeParameters@algICP(this, F0);
    
    utlDebugMsg( this.bEnableDebug,'initializing ICP parameters...\n')
    this.mexAlg_PIMLOP.ICP_InitializeParameters( F0 );
    utlDebugMsg( this.bEnableDebug,' ...done\n');
  end 
  
  %% ICP: Compute Matches
  function [ Yp, Yn ] = ICP_ComputeMatches( this )
    % compute matches
    %  also computes post-match parameters
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    [this.Yp, this.Yn] = this.mexAlg_PIMLOP.ICP_ComputeMatches();
    Yp = this.Yp;
    Yn = this.Yn;
    utlDebugMsg( this.bEnableDebug,' ...done\n' );  
  end

  %% ICP: Register Matches
  function [F,dRrod,dt] = ICP_RegisterMatches( this )
   
    % save prior registration
    Fprior = this.F;
    
    % register matches
    %  also computes post-registration parameters
    utlDebugMsg( this.bEnableDebug,'registering matches...\n' );
    [this.F] = this.mexAlg_PIMLOP.ICP_RegisterMatches();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );    
    
    % compute change in registration
    dF = this.F * invFrm3(Fprior);
    this.dRrod = rot2rodrigues(getRot(dF));
    this.dt = getPos(dF);
    
    % return values
    F = this.F;
    dRrod = this.dRrod;
    dt = this.dt;
    
  end

  
  %% ICP: Evaluate Error Function
  function fval = ICP_EvaluateErrorFunction( this )
    utlDebugMsg( this.bEnableDebug,'evaluating error function...\n' );
    this.fval = this.mexAlg_PIMLOP.ICP_EvaluateErrorFunction();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    % return values
    fval = this.fval;
  end
  
end

end
