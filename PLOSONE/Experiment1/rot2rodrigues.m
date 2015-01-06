function [rod] = rot2rodriguez(rot)

theta = acos((rot(1,1) + rot(2,2) + rot(3,3) - 1)/2);

% protect from divide by zero
% (happens when rot = I)
rod = zeros(3,1);
if (theta >= 1e-10)
    rod(1) = (theta / (2*sin(theta))) * (rot(3,2) - rot(2,3));
    rod(2) = (theta / (2*sin(theta))) * (rot(1,3) - rot(3,1));
    rod(3) = (theta / (2*sin(theta))) * (rot(2,1) - rot(1,2));
end

end