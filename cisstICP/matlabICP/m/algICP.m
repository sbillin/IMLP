classdef algICP < handle
% algorithm class for using Matlab optimizer

properties
  
  IterData;  % per iteration data  (type objIterData)
  
  % registration parameters
  F0;     % initial guess of registration
  
  F;      % current registration estimate
  dRrod;  % incremental change in rotation (Rodrigues vector)
  dt;     % incremental change in translation
  
  fval;   % current value of cost function

  % for callbacks
  bPrintToScreen;
  
end

methods (Abstract)
  
  %--- ICP Loop Routines ---%
  
  % inputs:
  %  [R0,t0] ~ initial registration parameter estimates  (3 x 4 double)
  %ICP_InitializeParameters( this, R0, t0 );
  
  % Post-Match parameter update and post-match filtering
  %  must take place in this call as well
  % output:
  %  Yp ~ match position (Nx3) (for default plotting feature)
  ICP_ComputeMatches( this );
  
  % Post-Register parameter update must take place in this call as well
  % output:
  %  Freg = [R,t]'   (3 x 4 double)
  ICP_RegisterMatches( this );
  
  % output:
  %  fval ~ error function value  (double scalar)
  ICP_EvaluateErrorFunction( this );
  
end % abstract methods

methods
  
  %% Constructor
  function this = algICP()
    this.bPrintToScreen = true;   % default    
  end        
  
  %% Initialize this base class
  function Initialize( this )
    % create object for storing per iteration data
    this.IterData = objIterData();
  end
  
  %% ICP: Initialize Registration
  % algPayload ~ optional algorithm-dependent parameters
  function ICP_InitializeParameters( this, F0, algPayload )
    this.F0 = F0;
    this.F = F0;
    this.dRrod = rot2rodrigues(getRot(F0));
    this.dt = getPos(F0);
  end
  
  %% Iteration Callback
  % iter ~ object of type objIterData
  function ICP_IterationCallback( this, iter, fid )
    
    % save iteration data
    dAng = norm(this.dRrod);
    dPos = norm(this.dt);
    this.IterData.AddIterData( iter, this.fval, this.dRrod, this.dt,...
      dAng, dPos);

    % form iteration message
    str = sprintf('iter %2u   %.3f   %.3f %.3f\n',...
      iter, this.fval, dAng*180/pi, dPos);
    
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
  
  
  %% ICP: algorithm-specific termination test (default method)
  % output:
  %  term   scalar:  1 ~ terminate 
  %                 -1 ~ do not terminate
  %                  0 ~ nothing special   
  function term = ICP_Terminate( this )
    term = 0;   % do nothing special
  end  
  
  %ICP_IterationCallback( this );
  %ComputeMatchDistance( this );  
end

end
