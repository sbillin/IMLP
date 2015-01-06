function [R, t] = ComputeRigidBodyXfm_Iso(X, Y)
% Arun's method for rigid body translation
%
numPts = size(X,1);

Xbar = mean(X,1)';
Ybar = mean(Y,1)';
Xp = X'-repmat(Xbar,1,numPts);
Yp = Y'-repmat(Ybar,1,numPts);
H = Xp*Yp';
% compute decomposition of H = U*S*V'
[U, S, V] = svd(H);
R = V*diag([1, 1, det(V*U)])*U';
t = Ybar - R*Xbar;

end