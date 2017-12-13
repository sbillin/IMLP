function [ P ] = getPos( F )
% Returns the position component of a 3D homogeneous transformation matrix
%  F = 4x4
P = F(1:3,4);
end