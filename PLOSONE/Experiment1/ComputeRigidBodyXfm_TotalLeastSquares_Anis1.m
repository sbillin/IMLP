function [data] = ComputeRigidBodyXfm_TotalLeastSquares_Anis1(xobs,yobs,M,maxIter,threshAng,threshPos,isoInit,loopCovDecomp)
% Algorithm: Total-Least-Squares registration of corresponding point sets
%
%  ** assumes anisotropic noise in only the My point set; assumes the
%     inverse covariance has been pre-computed before calling this method
%
%  xobs ~ Nx3
%  yobs ~ Nx3
%  Mx ~ 3x3xN
%  My ~ 3x3xN
%

numPts = size(xobs,1);
Minv = zeros(size(M));

J = zeros(3*numPts,6);
Jt_Minv = zeros(6,3*numPts);

% Step 1: initialize
R = eye(3);
t = zeros(3,1);
J(:,4:6) = repmat(-eye(3),[numPts,1]);

if (isoInit)
  % initialize transform parameters to isotropic values
  [R,t] = ComputeRigidBodyXfm_Iso(xobs,yobs);
end

% compute inverse covariances
if (loopCovDecomp)
  % slow method
  for i=1:numPts
    Minv(:,:,i) = inv(M(:,:,i));
  end  
else
  % fast method
  idxI = repmat(reshape(1:3*numPts,3,1,numPts),[1 3 1]);
  idxJ = repmat(reshape(1:3*numPts,1,3,numPts),[3 1 1]);
  sparseM = sparse(idxI(:),idxJ(:),M(:));
  Minv = reshape((sparseM \ repmat(eye(3),numPts,1))', 3,3,numPts);
end

numIter = 0;
dTheta = 1;
dt_norm = 1;

while ( ((dTheta > threshAng) | (dt_norm > threshPos)) & (numIter < maxIter))

  % Step 2: F0
  Rxobs = xobs*R';
  Txobs = Rxobs + repmat(t',[numPts,1]);
  Txobs_t = Txobs';
  yobs_t = yobs';
  F0 = yobs_t(:) - Txobs_t(:);  % stack into a single vector [f1x f1y f1z f2x...]'  
  
  % Step 3: J & Fx
  %   Note: Fx is block diagonal with all sub-blocks being equal
  %         to each other (equal to Fxi)
  J(:,1:3) = skewSet(Rxobs);
  Fxi = -R;
  
  % Step 4: solve dP by least squares
  
  % These are fixed for noise in 1 point set
  % % update covariances
  % M = mtimesx(Fxi,mtimesx(Mx,Fxi,'t')) + My;
  % % compute inverse covariances
  % idxI = repmat(reshape(1:3*numPts,3,1,numPts),[1 3 1]);
  % idxJ = repmat(reshape(1:3*numPts,1,3,numPts),[3 1 1]);
  % sparseM = sparse(idxI(:),idxJ(:),M(:));
  % Minv = reshape((sparseM \ repmat(eye(3),numPts,1))', 3,3,numPts);
  
  Jt_Minv(1:3,:) = reshape(mtimesx(reshape(J(:,1:3)',3,3,numPts),Minv),3,3*numPts);
  Jt_Minv(4:6,:) = -reshape(Minv,3,3*numPts);
  Jt_Minv_J = Jt_Minv * J;
  neg_Jt_Minv_F0 = Jt_Minv * (-F0);
  
  % Solve: Ax = b
  % Cholesky
  %  A = C'C
  %  C'Cx=b  =>  solve C'y=b  then  Cx=y
  C = chol(Jt_Minv_J);
  y = C'\neg_Jt_Minv_F0;  % forward substitution
  dP = C\y;               % back substitution
  dAlpha = dP(1:3);
  dt = dP(4:6);

  % Step 5: [R,t]
  dR = rodrigues2rot(dAlpha);
  R = dR*R;
  t = t + dt;
  
  [~, dTheta] = rot2AxisAngle(dR);
  dTheta = dTheta*180/pi;
  dt_norm = norm(dt);
  
  numIter = numIter + 1;
end 

data{1} = [R t];
data{2} = numIter;

end







