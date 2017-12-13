function [ F ] = getFrm3(R, t, renormalize)
% Creates a 3D homogeneous coordinate transformation
%  from a rotation matrix and position vector
% rot 3x3
% pos 3x1

if ~exist('renormalize','var') || isempty(renormalize)
  renormalize = false;
end

if renormalize
  R = RenormalizeRotation(R);
end

F = [[R; 0 0 0], [t; 1]];
end

