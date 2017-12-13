function [F, avgResidualDist, residualDist, scale] = ...
  Register_P2P(X, Y, bComputeResiduals, bEstimateScale)
% Arun's method for rigid body translation
%  X ~ nx3
%  Y ~ nx3
%
%  Solves yi = F*[xi;1] = R*xi + t
%  

if ~exist('bComputeResiduals','var') || isempty(bComputeResiduals)
  bComputeResiduals = false;
end
  
if ~exist('bEstimateScale','var') || isempty(bEstimateScale)
  bEstimateScale = false;
end

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
F = getFrm3(R,t);

if bEstimateScale
  % compute scale
  Hx = Xp*Xp'/numPts;
  Hy = Yp*Yp'/numPts;  
  Lx = eig(Hx);
  Ly = eig(Hy);
  scale = sqrt((Lx'*Ly)/(Lx'*Lx));
else
  scale = 1;
end

avgResidualDist = [];
residualDist = [];
if (bComputeResiduals)
  % compute average residual distance
  Xxfm = bsxfun(@plus, X*R', t');
  residualDist = sqrt(sum((Y - Xxfm).^2,2));
  avgResidualDist = mean(residualDist);
end

end