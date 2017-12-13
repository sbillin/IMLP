function [ Rn ] = RenormalizeRotation( R )
%
% renormalizes rotation to ensure properly constrained matrix
%

% renormalize rotation
[u,~,v] = svd(R);
Rn = u*v';

end