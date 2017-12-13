classdef algOptimizer_Matlab_SE3Rodrigues < handle
% algorithm class for using Matlab optimizer to solve
%  a rigid body (SE3) registration using a rodrigues
%  parameterization of rotation; estimation of scale
%  is also optional

properties
  
  options;    % the standard Matlab optimization options (optimoptions)
  
  bEstScale;  % enable / disable scale estimation
  bConstrainScale;  % enable / disable scale constraints
  %scaleConstraints; % min / max bounds on scale
  
end

methods
  %% constructor
  function this = algOptimizer_Matlab_SE3Rodrigues( options, bEstScale, bConstrainScale ) %scaleConstraints )
    if ~exist('bEstScale','var') || isempty(bEstScale)
      bEstScale = [];
    end
    if ~exist('bConstrainScale','var') || isempty(bConstrainScale)
      bConstrainScale = [];
    end    
    % if ~exist('scaleConstraints','var') || isempty(scaleConstraints)
    %   scaleConstraints = [];
    % end
    this.Initialize(options,bEstScale,bConstrainScale);    
  end
  
  %% initialize optimization options
  function Initialize( this, options, bEstScale, bConstrainScale) %scaleConstraints )    
    
    this.bEstScale = bEstScale;    
    this.bConstrainScale = bConstrainScale;
    % this.scaleConstraints = scaleConstraints;
    % if ~isempty(scaleConstraints) && length(scaleConstraints) ~= 2
    %     error('invalid scale constraints')
    % end
    
    if ~exist('options','var') || isempty(options)
      if bConstrainScale
        % define default unconstrained options
        options = optimoptions('fminunc',     ...
          'GradObj',          'on',             ...
          'Hessian',          'off',            ...
          'DerivativeCheck',  'off',            ...
          'FunValCheck',      'off',            ...
          'FinDiffType',      'central',        ...
          'Display',          'off',            ...
          'MaxIter',          20,               ...      
          'Algorithm',        'quasi-newton');
          %'Algorithm',        'trust-region');
          %'TypicalX',         [pi/180,pi/180,pi/180,1,1,1], ...    
          %'PlotFcns',         @optimplotfval,      ...
          %'Display',          'off',  ...
          %'Display',          'iter-detailed',  ...
          %'Display',          'final',  ...
          %'Algorithm',        'quasi-newton')
          %'Algorithm',        'trust-region')
          %'TolFun',           1e-6,             ...
          %'TolX',             1e-6,             ...
          %'FinDiffType',      'forward',         ...
          %'FinDiffType',      'central',         ...
          %'Diagnostics',      'on',             ...
          %'DiffMinChange',    1e-7,             ...
          %'TypicalX',         [1,1,1]*0.1*pi/180, ...    
      else
        % define default constrained options
        options = optimoptions('fmincon',     ...
          'GradObj',          'on',             ...
          'DerivativeCheck',  'off',            ...
          'FunValCheck',      'off',            ...
          'FinDiffType',      'central',        ...
          'Display',          'off',            ...
          ...%'Hessian',          'off',            ...
          ...%'MaxIter',          20,               ...      
          'Algorithm',        'active-set');
          % 'Algorithm',        'trust-region-reflective');          
          % 'Algorithm',        'sqp');          
          % 'Algorithm',        'interior-point');          
      end
    end
    
    % set optimizer options
    this.options = options;
  end
  
  %% run optimizer
  %
  % input:
  %   func  ~ cost function to be optimized
  %   x0    ~ initial guess = [alpha;t;s]
  %           alpha = Rodrigues rotation vector
  %           t = translation vector
  %           s = scale [optional]
  %
  % output:
  %   F     ~ final xfm (homogeneous form)
  %   alpha ~ final rotation (Rodrigues form)
  %   fval  ~ final cost function value
  %   grad  ~ final gradient value
  %   scale ~ scale
  %
  function [F,alpha,fval,grad,scale] = Register( this, func, x0, scaleConstraints )
    
    if this.bConstrainScale && ~this.bEstScale
      error('Scale constraints are active, but scale estimation is inactive')
    end
    if ~exist('scaleConstraints','var') || isempty(scaleConstraints)
      if this.bConstrainScale
        error('Scale constraints are active, but none were provided')
      end
    else
      if ~this.bConstrainScale
        error('Scale constraints are not active, but yet were provided')
      end
      if any(size(scaleConstraints) ~= [2,1]) && any(size(scaleConstraints) ~= [1,2])
        error('Invalid scale constraints')
      end
    end
    
    if ~exist('func','var') || isempty(func)
      error('ERROR: 1st argument must be pointer to cost function');
    end
      
    if ~exist('x0','var') || isempty(x0)
      % begin registration from zero offset
      if this.bEstScale
        x0 = [zeros(6,1);1];
      else
        x0 = zeros(6,1);
      end
    end    
    
    if ~this.bConstrainScale
      % compute new x to minimize the unconstrained cost function
      [x,fval,exitflag,output,grad] = ...
        fminunc( func, x0, this.options );
    else
      % define lower / upper bounds on the variable values
      lb = -Inf(size(x0));
      ub = Inf(size(x0));
      lb(7) = scaleConstraints(1);
      ub(7) = scaleConstraints(2);
      % compute new x to minimize the constrained cost function
      [x,fval,exitflag,output,labmda,grad] = ...
        fmincon( func, x0, [],[],[],[],lb,ub,[], this.options );
    end
    
    % convert result to matrix form F = [R,t]
    alpha = x(1:3);
    R = rodrigues2rot(x(1:3));
    t = x(4:6);
    F = getFrm3(R,t);
    
    if this.bEstScale
      scale = x(7);
    else
      scale = [];
    end    
  end
end

end
