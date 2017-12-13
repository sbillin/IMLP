classdef algDirICP_StdICP < algICP

properties

  bEnableDebug;
  bTargetAsMesh;
  
  % sample data
  Xp;          % source pts        N x 3
  Xn;          % source norms      N x 3
  %Xp_xfm;      % xfmd source pts  N x 3
  
  Yp;          % match pts         N x 3
  Yn;          % match normals 3d  N x 3
  
  % target shape
  % mesh and point cloud target:
  V;        % vertices            Nv x 3
  % mesh target only:
  T;        % triangle faces (vertex indices)   Nt x 3
  Tn;       % triangle normals                  Nt x 3
  
  % member objects
  mexAlg_StdICP;   % mex object implementing this algorithm
  
end

methods
  
  %% Constructor
  function [this] = algDirICP_StdICP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlg_StdICP.delete();
  end  
  
  %% Initialize Algorithm
  %   V     ~ target vertices   Nv x 3
  %   T,Tn  ~ target triangles and triangle normals   Nt x 3
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  % bTargetAsMesh   ~ set point cloud or mesh target
  % bEnableDebug    ~ enable debug messages
  function Initialize( this, ...
    V,T,Tn, ...
    Xp,Xn,...
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
    if ~exist('Xn','var') || isempty(Xn)
      warning('Missing input Xn')
      this.Xn = [];
    else
      this.Xn = Xn;
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
      utlDebugMsg( this.bEnableDebug,'Creating mex object: mexInterface_AlgDirICP_StdICP_Mesh...\n')
      this.mexAlg_StdICP = mexInterface_AlgDirICP_StdICP_Mesh();
      utlDebugMsg( this.bEnableDebug,' ...done\n');
      
      utlDebugMsg( this.bEnableDebug,'Initializing mexAlg_StdICP...\n')
      this.mexAlg_StdICP.Initialize( ...
        this.V, this.T, this.Tn,...
        this.Xp, this.Xn );
      utlDebugMsg( this.bEnableDebug,' ...done\n');
    else
      error('Point Cloud version of mexalgDirICP_StdICP not ready');
    end    
  end  
  
  %% Set Samples
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  function SetSamples( this, Xp,Xn )  
    this.Xp = Xp;
    this.Xn = Xn;
    this.mexAlg_StdICP.SetSamples( Xp,Xn );
  end
  
  %% ICP: Initialize Parameters
  function ICP_InitializeParameters( this, F0 )
    
    % initial guess
    if ~exist('F0','var') || isempty(F0)
      F0 = eye(4);
    end
    
    % initialize superclass
    ICP_InitializeParameters@algICP(this, F0);
    
    utlDebugMsg( this.bEnableDebug,'initializing ICP parameters...\n')
    this.mexAlg_StdICP.ICP_InitializeParameters( F0 );
    utlDebugMsg( this.bEnableDebug,' ...done\n');
  end 
  
  %% ICP: Compute Matches
  function [ Yp, Yn ] = ICP_ComputeMatches( this )
    % compute matches
    %  also computes post-match parameters
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    [this.Yp, this.Yn] = this.mexAlg_StdICP.ICP_ComputeMatches();
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
    [this.F] = this.mexAlg_StdICP.ICP_RegisterMatches();
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
    this.fval = this.mexAlg_StdICP.ICP_EvaluateErrorFunction();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    % return values
    fval = this.fval;
  end
  
end

end
