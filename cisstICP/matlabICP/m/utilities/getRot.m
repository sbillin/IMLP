function [ R ] = getRot( F )
% Returns the rotation component of a 3D homogeneous transformation matrix
%  F = 4x4
R = F(1:3,1:3);
end

