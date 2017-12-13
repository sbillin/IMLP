function [rot] = euler2rot(yaw,pitch,roll)
% Converts from yaw/pitch/roll euler angles (azimuth/elevation/roll) to a
%  rotation matrix
% yaw/pitch/roll = 1xm
% rot = 3x3xm
sz = size(yaw);
rot = zeros(3,3,sz(2));
for i=1:sz(2)
% Multiply using a local frame-of-reference
%  (i.e. first rotation goes furthest to the left)
rot(:,:,i) = rotz(yaw(i))*roty(pitch(i))*rotx(roll(i));
end