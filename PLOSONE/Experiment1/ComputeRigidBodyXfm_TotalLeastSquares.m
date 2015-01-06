function [data] = ComputeRigidBodyXfm_TotalLeastSquares(xobs,yobs,Mx,My,maxIter,threshAng,threshPos,isoInit,loopCovDecomp)
% Algorithm: Total-Least-Squares registration of corresponding point sets
%
%  xobs ~ Nx3
%  yobs ~ Nx3
%  Mx ~ 3x3xN
%  My ~ 3x3xN
%

numPts = size(xobs,1);

J = zeros(3*numPts,6);
Jt_Minv = zeros(6,3*numPts);
% vTheta = zeros(maxIter,1);
% vAxis = zeros(maxIter,3);
% vTrans = zeros(maxIter,3);

if ~exist('loopCovDecomp','var') || isempty(loopCovDecomp)
  loopCovDecomp = 0;   % default to fast method
end

% Step 1: initialize
R = eye(3);
t = zeros(3,1);
J(:,4:6) = repmat(-eye(3),[numPts,1]);

if (isoInit)
  % initialize transform parameters to isotropic values
  [R,t] = ComputeRigidBodyXfm_Iso(xobs,yobs);
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
  % update covariances
  M = mtimesx(Fxi,mtimesx(Mx,Fxi,'t')) + My;
  
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
  %[dAxis, dTheta] = rot2AxisAngle(dR);
  dTheta = dTheta*180/pi;
  dt_norm = norm(dt);
  
  numIter = numIter + 1;
  % vTheta(numIter) = dTheta;
  % vAxis(numIter,:) = dAxis';
  % vTrans(numIter,:) = dt';
end 

data{1} = [R t];
data{2} = numIter;
% data{3} = vTheta(1:numIter)';
% data{4} = vAxis(1:numIter,:);
% data{5} = vTrans(1:numIter,:);
end







