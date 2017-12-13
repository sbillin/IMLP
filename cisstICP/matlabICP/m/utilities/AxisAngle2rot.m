function [ R ] = AxisAngle2rot( axis, angle, bRenormalize )
     
if ~exist('bRenormalize','var') || isempty(bRenormalize)
  bRenormalize = false;
end

if bRenormalize
  axis = axis / norm(axis);
end

% if rotation less than very small number, return identity
%  to prevent NaN's in calculation
if (angle < eps)
    R = eye(3);
else
    I = eye(3);
    kSkew = [...
        0       -axis(3)  axis(2)
        axis(3)  0        -axis(1)
        -axis(2) axis(1)  0];
    R = I + kSkew*sin(angle) + (1-cos(angle))*(kSkew*kSkew);
    %rot = I*cos(angle) + kSkew*sin(angle) + (1-cos(angle))*(axis*axis')
end

end