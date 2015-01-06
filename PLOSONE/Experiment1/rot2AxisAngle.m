function [axis, angle] = rot2AxisAngle(rot)

rod = rot2rodrigues(rot);

% protect for divide by zero
if (norm(rod) >= 1e-10)
  axis = rod / norm(rod);
  angle = norm(rod);    
else
  axis = [0;0;1];
  angle = 0;
end

end