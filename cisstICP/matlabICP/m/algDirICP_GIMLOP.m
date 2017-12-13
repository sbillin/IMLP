classdef algDirICP_GIMLOP < algICP

properties

  bEnableDebug;
  bTargetAsMesh;
  
  % sample data
  Xp;          % source pts        N x 3
  Xn;          % source norms      N x 3
  %Xp_xfm;      % xfmd source pts  N x 3
  
  Yp;          % match pts         N x 3
  Yn;          % match normals 3d  N x 3
  
  M;    % covariance of source positions          3 x 3 x N
  k;    % concentration of source orientations    N x 1
  E;    % eccentricity of source orientations     N x 1
        %  E ~ [0,1);  B = E * k/2
  L;    % Kent major/minor axis                   2 x 3 x N
  dynamicParamEst;    % 0 ~ fixed params (int)
                      % 1 ~ dynamic threshold params (int)
                      % 2 ~ fixed params; 1st iter positions only (int)
  
  % monitor variables for dynamic noise model updates
  meanSigma2;
  meanK;
  meanE;
  
  % target shape
  % mesh and point cloud target:
  V;        % vertices            Nv x 3
  % mesh target only:
  T;        % triangle faces (vertex indices)   Nt x 3
  Tn;       % triangle normals                  Nt x 3
  
  % member objects
  mexAlg_GIMLOP;   % mex object implementing this algorithm
  
end

methods
  
  %% Constructor
  function [this] = algDirICP_GIMLOP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlg_GIMLOP.delete();
  end  
  
  %% Initialize Algorithm
  %   V     ~ target vertices   Nv x 3
  %   T,Tn  ~ target triangles and triangle normals   Nt x 3
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  %   M ~ covariance matrices of soure points     3 x 3 x Ns
  %   k ~ concentration of source orientations    Ns x 1
  %   E ~ eccentricities of source orientations   Ns x 1
  %   L ~ Kent major/minor axis                   2 x 3 x Ns
  %   dynamicParamEst (int)  0 ~ fixed params
  %                          1 ~ dynamic threshold params  
  % bTargetAsMesh   ~ set point cloud or mesh target
  % bEnableDebug    ~ enable debug messages
  function Initialize( this, ...
    V,T,Tn, ...
    Xp,Xn,M,k,E,L, dynParamEst,...
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
    if ~exist('M','var') || isempty(M)
      warning('Missing input M')
      this.M = [];
    else
      this.M = M;
    end     
    if ~exist('k','var') || isempty(k)
      warning('Missing input k')
      this.k = [];
    else
      this.k = k;
    end
    if ~exist('E','var') || isempty(E)
      warning('Missing input E')
      this.E = [];
    else
      if E < 0 | E >= 1
        error(['invalid value E: ', num2str(E)]);
      end
      this.E = E;
    end         
    if ~exist('L','var') || isempty(L)
      warning('Missing input L')
      this.L = [];
    else
      nPts = size(L,3);
      tol = 1e-8;
      for i = 1:nPts
        d1 = Xn(i,:)*L(1,:,i)';
        d2 = Xn(i,:)*L(2,:,i)';
        d3 = L(1,:,i)*L(2,:,i)';
        if d1 > tol || d2 > tol || d3 > tol
          error(['L matrix not perpendicular for index: ',num2str(i)]);
        end
      end      
      this.L = L;
    end  
    if ~exist('dynParamEst','var') || isempty(dynParamEst)
      warning('Missing input dynamicParamEst defaulting to: 0')
      this.dynamicParamEst = int32(0);
    else
      this.dynamicParamEst = int32(dynParamEst);
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
      utlDebugMsg( this.bEnableDebug,'Creating mex object: mexInterface_AlgDirICP_GIMLOP_Mesh...\n')
      this.mexAlg_GIMLOP = mexInterface_AlgDirICP_GIMLOP_Mesh();
      utlDebugMsg( this.bEnableDebug,' ...done\n');
      
      utlDebugMsg( this.bEnableDebug,'Initializing mexAlg_GIMLOP...\n')
      this.mexAlg_GIMLOP.Initialize( ...
        this.V, this.T, this.Tn,...
        this.Xp, this.Xn, this.M, this.k, this.E, this.L, this.dynamicParamEst);
      utlDebugMsg( this.bEnableDebug,' ...done\n');
    else
      error('Point Cloud version of mexAlgDirICP_Kent not ready');
    end    
  end  
  
  %% Set Samples
  %   Xp    ~ source points             Ns x 3
  %   Xn    ~ source normals            Ns x 3
  %   M ~ covariance matrices of soure points     3 x 3 x Ns
  %   k ~ concentration of source orientations    Ns x 1
  %   E ~ eccentricities of source orientations   Ns x 1
  %   L ~ Kent major/minor axis                   2 x 3 x Ns
  %   dynamicParamEst (int)  0 ~ fixed params
  %                          1 ~ dynamic threshold params  
  function SetSamples( this, Xp,Xn, M,k,E,L, dynParamEst )  
    if E < 0 | E >= 1
      error(['invalid value E: ', num2str(E)]);
    end    
    nPts = size(L,3);
    tol = 1e-8;
    for i = 1:nPts
      d1 = Xn(i,:)*L(1,:,i)';
      d2 = Xn(i,:)*L(2,:,i)';
      d3 = L(1,:,i)*L(2,:,i)';
      if d1 > tol || d2 > tol || d3 > tol
        error(['L matrix not perpendicular for index: ',num2str(i)]);
      end
    end
    
    this.Xp = Xp;
    this.Xn = Xn;
    this.M = M;
    this.k = k;
    this.E = E;
    this.L = L;
    this.dynamicParamEst = int32(dynParamEst);
    this.mexAlg_GIMLOP.SetSamples( Xp,Xn, M,k,E,L, dynParamEst);
  end
  
  %% ICP: Initialize Registration
  function ICP_InitializeParameters( this, F0 )
    
    % initial guess
    if ~exist('F0','var') || isempty(F0)
      F0 = getFrm3(eye(3),zeros(3,1));
    end
    
    % initialize superclass
    ICP_InitializeParameters@algICP(this, F0);
    
    utlDebugMsg( this.bEnableDebug,'initializing ICP parameters...\n')
    this.mexAlg_GIMLOP.ICP_InitializeParameters( F0 );
    utlDebugMsg( this.bEnableDebug,' ...done\n');
  end 
  
  %% ICP: Compute Matches
  function [ Yp, Yn ] = ICP_ComputeMatches( this )
    % compute matches
    %  also computes post-match parameters
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    [this.Yp, this.Yn] = this.mexAlg_GIMLOP.ICP_ComputeMatches();
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
    [this.F, this.meanSigma2, this.meanK, this.meanE] = ...
      this.mexAlg_GIMLOP.ICP_RegisterMatches();
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
    this.fval = this.mexAlg_GIMLOP.ICP_EvaluateErrorFunction();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    % return values
    fval = this.fval;
  end
  
  %% Debug: Set Matches
  function Debug_SetMatches( this, Yp, Yn )
    utlDebugMsg( this.bEnableDebug,'debug: forcing matches...\n' );
    this.mexAlg_GIMLOP.Debug_SetMatches(Yp,Yn);
    utlDebugMsg( this.bEnableDebug,' ...done\n' );    
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
      str = sprintf('%45s\n','(sigma2/k/e)');
      str = [str,sprintf('iter %2u   %.3f   %.3f %.3f\n',...
        iter, this.fval, dAng*180/pi, dPos )];
    else
      str = sprintf('iter %2u   %.3f   %.3f %.3f   %.2f %.0f %0.2f\n',...
        iter, this.fval, dAng*180/pi, dPos, ...
        this.meanSigma2, this.meanK, this.meanE );      
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
