function [R] = Compute3dRotation_Iso(X, Y)
% Arun's method for rotation
%

H = X'*Y;
% compute decomposition of H = U*S*V'
[U, S, V] = svd(H);
R = V*diag([1, 1, det(V*U)])*U';

end