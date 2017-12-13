classdef algDirICP_IMLOP < algICP

properties

  bEnableDebug;
  bTargetAsMesh;
  
  % sample data
  Xp;          % source pts        N x 3
  Xn;          % source norms      N x 3
  %Xp_xfm;      % xfmd source pts  N x 3
  
  Yp;          % match pts         N x 3
  Yn;          % match normals 3d  N x 3
  
  sigma2_init; % variance of source positions         (double)
  kinit;       % concentration of source orientations (double)
  wRpos;       % weight of positions in estimating param k
  bDynamicParamEst;   % dynamic estimation of params
  
  % monitor variables for dynamic noise model estimates
  sigma2;
  k;
  
  % target shape
  % mesh and point cloud target:
  V;        % vertices            Nv x 3
  % mesh target only:
  T;        % triangle faces (vertex indices)   Nt x 3
  Tn;       % triangle normals                  Nt x 3
  
  % member objects
  mexAlg_IMLOP;   % mex object implementing this algorithm
  
end

methods
  
  %% Constructor
  function [this] = algDirICP_IMLOP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlg_IMLOP.delete();
  end  
  
  %% Initialize Algorithm
  %   V     ~ target vertices   Nv x 3
  %   T,Tn  ~ target triangles and triangle normals   Nt x 3
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  %   sigma2_init       double
  %   kinit             double
  %   wRpos             double
  %   bDynamicParamEst  bool
  % bTargetAsMesh   ~ set point cloud or mesh target
  % bEnableDebug    ~ enable debug messages
  function Initialize( this, ...
    V,T,Tn, ...
    Xp,Xn, ...
    sigma2_init, kinit, wRpos, bDynamicParamEst, ...
    bTargetAsMesh, bEnableDebug )
    
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
      if size(Xn,2) ~= 3
        error('Invalid size of Xn')
      end      
      this.Xn = Xn;
    end
    if ~exist('sigma2_init','var') || isempty(sigma2_init)      
      this.sigma2_init = 1.0;
      warning(['Missing input sigma2_init; setting default: ',...
        num2str(this.sigma2_init)])
    else
      this.sigma2_init = sigma2_init;
    end     
    if ~exist('kinit','var') || isempty(kinit)
      this.kinit = 0;
      warning(['Missing input kinit; setting default: ',...
        num2str(this.kinit)])      
    else
      this.kinit = kinit;
    end
    if ~exist('wRpos','var') || isempty(wRpos)
      this.wRpos = 0.5;
      warning(['Missing input wRpos; setting default: ',...
        num2str(this.wRpos)])         
    else
      this.wRpos = wRpos;
    end 
    if ~exist('bDynamicParamEst','var') || isempty(bDynamicParamEst)
      this.bDynamicParamEst = true;
      warning(['Missing input bDynamicParamEst; setting default: ',...
        num2str(this.bDynamicParamEst)])           
    else
      this.bDynamicParamEst = bDynamicParamEst;
    end
        
    % target shape
    if ~exist('bTargetAsMesh','var') || isempty(bTargetAsMesh)
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
      utlDebugMsg( this.bEnableDebug,'Creating mex object: mexInterface_AlgDirICP_IMLOP_Mesh...\n')
      this.mexAlg_IMLOP = mexInterface_AlgDirICP_IMLOP_Mesh();
      utlDebugMsg( this.bEnableDebug,' ...done\n');
      
      utlDebugMsg( this.bEnableDebug,'Initializing mexAlg_IMLOP...\n')
      this.mexAlg_IMLOP.Initialize( ...
        this.V, this.T, this.Tn,...
        this.Xp, this.Xn,...
        this.sigma2_init, this.kinit, this.wRpos, this.bDynamicParamEst);
      utlDebugMsg( this.bEnableDebug,' ...done\n');
    else
      error('Point Cloud version of mexalgDirICP_vMFG not ready');
    end    
  end  
  
  %% Set Samples
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  %   sigma2_init       double
  %   kinit             double
  %   wRpos             double
  %   bDynamicParamEst  bool
  function SetSamples( this, Xp,Xn, sigma2_init,kinit,wRpos,bDynamicParamEst )  
    this.Xp = Xp;
    this.Xn = Xn;
    if exist('sigma2_init','var') && ~isempty(sigma2_init)
      this.sigma2_init = sigma2_init;
    end
    if exist('kinit','var') && ~isempty(kinit)
      this.kinit = kinit;
    end
    if exist('wRpos','var') && ~isempty(wRpos)
      this.wRpos = wRpos;
    end
    if exist('bDynamicParamEst','var') && ~isempty(bDynamicParamEst)
      this.bDynamicParamEst = bDynamicParamEst;
    end
    this.mexAlg_IMLOP.SetSamples( this.Xp, this.Xn,...
      this.sigma2_init, this.kinit, this.wRpos, this.bDynamicParamEst );
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
    this.mexAlg_IMLOP.ICP_InitializeParameters( F0 );
    utlDebugMsg( this.bEnableDebug,' ...done\n');
  end 
  
  %% ICP: Compute Matches
  function [ Yp, Yn ] = ICP_ComputeMatches( this )
    % compute matches
    %  also computes post-match parameters
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    [this.Yp, this.Yn] = this.mexAlg_IMLOP.ICP_ComputeMatches();
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
    [this.F, this.sigma2, this.k] = this.mexAlg_IMLOP.ICP_RegisterMatches();
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
    this.fval = this.mexAlg_IMLOP.ICP_EvaluateErrorFunction();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    % return values
    fval = this.fval;
  end
  
  %% Iteration Callback
  function ICP_IterationCallback( this, iter, fid )
    
    % save iteration data
    dAng = norm(this.dRrod);
    dPos = norm(this.dt);
    this.IterData.AddIterData( iter, this.fval, this.dRrod, this.dt,...
      dAng, dPos);

    % form iteration message
    if iter == 0
      str = sprintf('%45s\n','(sigma2/k)');
      str = [str,sprintf('iter %2u   %.3f   %.3f %.3f\n',...
        iter, this.fval, dAng*180/pi, dPos )];
    else
      str = sprintf('iter %2u   %.3f   %.3f %.3f   %.2f %.0f\n',...
        iter, this.fval, dAng*180/pi, dPos, ...
        this.sigma2, this.k );      
    end
    
    % display message
    bWriteToFile = 1;
    if ~exist('fid','var') || isempty(fid)
      fid = [];
      bWriteToFile = 0;
    end    
    printfCapture(str,fid,bWriteToFile, this.bPrintToScreen);
    
    % %  1) to terminal
    % fprintf('%s',str);
    % %  2) to file
    % if exist('fid','var') && ~isempty(fid)      
    %   fprintf(fid,'%s',str);
    % end
  end      
  
end

end
