classdef algICP_TrimmedICP < algICP
  
properties

  bEnableDebug;
  bTargetAsMesh;
  bEstimateScale;
  bPlotIterations;
    
  % sample data
  X;          % source pts        N x 3
  Y;          % matches 3D pts    N x 3
  Xxfm;       % xfmd source pts   N x 3
  matchDatums;  % covariance tree datums of matches N x 3
  
  % target shape
  % mesh and point cloud target:
  V;        % vertices            Nv x 3
  % mesh target only:
  T;        % triangle faces (vertex indices)         Nt x 3
  Tn;       % triangle normals                        Nt x 3
  
  % target shape for plotting
  Vplot;
  Tplot;
  Tnplot;
  
  % outliers
  trimRate;   % percent of matches in decimal form (i.e. 0.3 for 30%)
  outliers;   % logical array of outliers
  
  % scale estimation
  scale;
  prevScale;
  scaleTermThresh;
  
  % member objects
  mexAlgPDTreeCP;   % mex object for StdICP algorithm
  
end

methods
  
  %% Constructor
  function [this] = algICP_TrimmedICP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlgPDTreeCP.delete();
  end  
  
  %% Initialize Algorithm
  % V       ~ target vertices     Nv x 3
  % T,Tn    ~ target triangles and triangle normals (for mesh target only)  Nt x 3
  % X       ~ source points       N x 3  
  % options
  %   bEnableDebug    ~ enable debug messages
  %   bTargetAsMesh   ~ set point cloud or mesh target
  %   bEstimateScale
  %   trimRate
  function Initialize( this, ...
    V,T,Tn, X, options, Vplot, Tplot, Tnplot )
    
    % initialize superclass
    Initialize@algICP( this );
        
    % debug enable
    if ~exist('options','var') || ~isfield(options,'bEnableDebug') || isempty(options.bEnableDebug)
      options.bEnableDebug = 0;
    end
    this.bEnableDebug = options.bEnableDebug;
    
    % trim rate
    if ~exist('options','var') || ~isfield(options,'trimRate') || isempty(options.trimRate)
      options.trimRate = 0;
    end
    this.trimRate = options.trimRate;
    
    % scale estimation
    if ~exist('options','var') || ~isfield(options,'bEstimateScale') || isempty(options.bEstimateScale)
      options.bEstimateScale = false;
    end
    this.bEstimateScale = options.bEstimateScale;      
    
    % plotting
    if ~exist('options','var') || ~isfield(options,'bPlotIterations') || isempty(options.bPlotIterations)
      options.bPlotIterations = false;
    end
    this.bPlotIterations = options.bPlotIterations;      
    
    % source shape
    if ~exist('X','var') || isempty(X)
      %error('Missing source shape input X')
      warning('Missing source shape input X')
      this.X = [];
    else      
      if size(X,2) ~= 3
        error('Invalid size of source shape (X)')
      end
      this.X = X;
    end
        
    % target shape
    if ~exist('options','var') || ~isfield(options,'bTargetAsMesh') || isempty(options.bTargetAsMesh)
      options.bTargetAsMesh = true;
    end
    if size(V,2) ~= 3
      error('Invalid size of target vertices (V)')
    end
    this.bTargetAsMesh = options.bTargetAsMesh;
    this.V = V;
    if (this.bTargetAsMesh)
      this.T = T;
      this.Tn = Tn;
    else
      this.T = [];
      this.Tn = [];
    end
    
    % target shape for plotting
    if ~exist('Vplot','var') || ~exist('Tplot','var') || ~exist('Tnplot','var') || ...
        isempty(Vplot) || isempty(Tplot) || isempty(Tnplot)
      Vplot = V;
      Tplot = T;
      Tnplot = Tn;
    end
    this.Vplot = Vplot;
    this.Tplot = Tplot;
    this.Tnplot = Tnplot;
        
    % Mex objects
    if (this.bTargetAsMesh)
      utlDebugMsg( this.bEnableDebug,'creating 3D point matcher algorithm...\n')
      this.mexAlgPDTreeCP = mexInterface_AlgPDTree_CP_Mesh();
      utlDebugMsg( this.bEnableDebug,' ...done\n')

      utlDebugMsg( this.bEnableDebug,'initializing 3D point matcher algorithm...\n')
      this.mexAlgPDTreeCP.Initialize( this.V, this.T, this.Tn);
      utlDebugMsg( this.bEnableDebug,' ...done\n')
    else
      error('point cloud target not supported')
    end    
  end  
  
  %% Set Samples
  % X  ~ source points       N x 3
  function SetSamples( this, X )  
    this.X = X;
  end
  
  %% ICP: Initialize Registration
  function ICP_InitializeParameters( this, F0, algPayload )       
    
    % initial alignment
    if ~exist('F0','var') || isempty(F0)
      F0 = getFrm3(eye(3),zeros(3,1));
    end
    
    % initial scale
    if ~exist('algPayload','var') || ~isfield(algPayload,'scale0') || isempty(algPayload.scale0)
      algPayload.scale0 = 1;
    end
    this.scale = algPayload.scale0;
    
    % scale termination threshold
    if ~exist('algPayload','var') || isempty(algPayload) ...
      || ~isfield(algPayload,'scaleTermThresh')    
      algPayload.scaleTermThresh = 0.0001;  % one-hundredth of a percent
    end
    this.scaleTermThresh = algPayload.scaleTermThresh;    
    
    % initialize superclass
    ICP_InitializeParameters@algICP(this, F0, algPayload);
    
    this.ComputeXfmdSamples( F0, this.scale );
    
    this.matchDatums = [];
  end 
  
  %% ICP: Compute Matches
  function [ Y ] = ICP_ComputeMatches( this )
    
    % compute matches
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    [this.Y, ~, this.matchDatums] = ...
      this.mexAlgPDTreeCP.ComputeMatches( this.Xxfm, this.matchDatums );
    Y = this.Y;
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    this.ComputeOutliers();
  end

  %% ICP: Register Matches
  function [F,dRrod,dt,dScale] = ICP_RegisterMatches( this )
   
    % save prior registration
    Fprior = this.F;
    scalePrior = this.scale;    

    % register matches
    %  also computes post-registration parameters
    utlDebugMsg( this.bEnableDebug,'registering matches...\n' );
    inliers = ~this.outliers;
    [this.F, ~,~, scaleEst] = Register_P2P( ...
      this.scale * this.X(inliers,:), this.Y(inliers,:), false, this.bEstimateScale );
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    if this.bEstimateScale
      this.scale = scaleEst * this.scale;
    end
    
    % compute change in registration
    dF = this.F * invFrm3(Fprior);
    this.dRrod = rot2rodrigues(getRot(dF));
    this.dt = getPos(dF);
    
    % update xfmd sample positions
    this.ComputeXfmdSamples( this.F, this.scale );
    
    % return values
    R = getRot(this.F);
    t = getPos(this.F);
    F = getFrm3(this.scale*R, t, false);
    dRrod = this.dRrod;
    dt = this.dt;
    dScale = (this.scale-scalePrior)/scalePrior;    
  end

  
  %% ICP: Evaluate Error Function
  function fval = ICP_EvaluateErrorFunction( this )    
    this.fval = 0;
    fval = this.fval;
  end
  
  %% ICP: Iteration Callback
  % iter ~ iteration number
  % fid  ~ handle for output file
  function ICP_IterationCallback( this, iter, fid )
    
    % save iteration data
    dAng = norm(this.dRrod);
    dPos = norm(this.dt);
    this.IterData.AddIterData( iter, this.fval, this.dRrod, this.dt,...
      dAng, dPos);

    % form iteration message
    str = sprintf('iter %2u   %.3f   %.3f %.3f  scale=%.3f\n',...
      iter, this.fval, dAng*180/pi, dPos, this.scale);
    
    % display message
    bWriteToFile = true;
    if ~exist('fid','var') || isempty(fid)
      fid = [];
      bWriteToFile = false;
    end    
    printfCapture(str,fid,bWriteToFile, this.bPrintToScreen);
    
    if this.bPlotIterations
      this.PlotIteration( iter )
    end        
  end  
  
  
  %% Helper Functions
    
  function ComputeXfmdSamples( this, F, scale )
    if ~exist('scale','var') || isempty(scale)
      scale = 1;
    end
    this.Xxfm = applyFrm3( F, scale * this.X );
  end
  
  % Outlier Handling
  function ComputeOutliers( this )
    % compute match distances
    res = this.Y - this.Xxfm;
    d = sqrt( sum(res.^2,2) );
    
    % sort match distances
    [dSort, dSortIdx] = sort(d,'descend');
    
    % remove worst matches as a percentage of all matches
    if (this.trimRate > 0)
      outlierIdx = dSortIdx( 1:ceil(this.trimRate*length(d)) );
    else
      outlierIdx = [];
    end
    this.outliers = false(length(d),1);
    this.outliers(outlierIdx) = true;
  end
  
  function PlotIteration( this, iter )
    if iter > 0
      [AZ,EL] = view(); % save current viewpoint (in case user changed it)
    end
        
    % plot mesh
    figure(1);
    PlotMesh(this.Vplot,this.Tplot,0.5,0.5);
    hold on
    % plot matches
    plot3(this.Y(:,1),this.Y(:,2),this.Y(:,3),'bo');
    % plot inliers
    plot3(this.Xxfm(~this.outliers,1),this.Xxfm(~this.outliers,2),this.Xxfm(~this.outliers,3),'r.','MarkerSize',12);
    % plot outliers
    if any(this.outliers)
      plot3(this.Xxfm(this.outliers,1),this.Xxfm(this.outliers,2),this.Xxfm(this.outliers,3),'m*','MarkerSize',4);
    end
    hold off
    axis equal;
    if exist('AZ','var')
      view([AZ,EL]);  % restore viewpoint
    end
    drawnow;
  end
  
end

end
