classdef objIterData < handle
% handle class for storing iteration data from a registration

properties

  % per iteration data
  fval;   % cost function value                               (iter + 1 x 1)
  dRrod;  % incremental rotation (Rodrigues form)             (iter + 1 x 3)
  dt;     % incremental translation                           (iter + 1 x 3)
  dAng;   % magnitude of incremental change in rotation       (iter + 1 x 1)
  dPos;   % magnitude of incremental change in translation    (iter + 1 x 1)

end

methods  
  
  %% Constructor
  function this = objIterData()
  end
  
  %% Methods
  function AddIterData( this, iter,fval,dRrod,dt,dAng,dPos )
    index = iter + 1;
    if ~isempty(fval)
      this.fval(index,1) = fval;
    end
    this.dRrod(index,:) = dRrod;
    this.dt(index,:) = dt;
    this.dAng(index,1) = dAng;
    this.dPos(index,1) = dPos;
  end
end
end
