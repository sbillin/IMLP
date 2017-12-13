classdef objOptICP < handle
% handle class for storing registration options

properties
 
  % termination condition
  maxIter;          % maximum iteration count
  term_holdIter;    % minimum number of iterations to maintain thresholds
  term_dAng;        % minimum threshold on change in rotation angle (radians)
  term_dPos;        % minimum threshold on change in translation    
 
end

methods  
  
  %% Constructor
  function this = objOptICP()
    % defaults
    this.maxIter = 100;
    this.term_dAng = 0.001*pi/180;
    this.term_dPos = 0.001;
    this.term_holdIter = 2;
  end
  
  %% Methods

end
end
