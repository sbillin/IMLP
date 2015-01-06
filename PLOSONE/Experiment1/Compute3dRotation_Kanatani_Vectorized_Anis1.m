function [data] = Compute3dRotation_Kanatani_Vectorized(xobs,yobs,Mx,My, maxIter,threshAng)
%
% Vectorized Version
%
% xobs ~ n x 3
% yobs ~ n x 3
% Mx ~ 3x3xn covariances
% My ~ 3x3xn covariances

numPts = size(xobs,1);

R = eye(3);
N = zeros(4,4);
Mx_m_My = Mx-My;
Mx_p_My = Mx+My;

% Step 1: Init Xa
Xa = zeros(3,4,numPts);
Xa(:,1,:) =  reshape((yobs-xobs)',3,1,numPts);
Xa(:,2:4,:) = skewSet2(yobs + xobs);

% Step 2: Init c, Wa
c = 0;
Wa = repmat(eye(3),1,1,numPts);

numIter = 0;
lambda = 1;
dTheta = 1;

%while ( abs(lambda) > 0.1 && numIter <= maxIter )    % Kanatani's termination condition
while ( dTheta > threshAng && numIter < maxIter)

  % Step 3: M
  M = sum(mtimesx(Xa,'t',mtimesx(Wa,Xa)),3);

  % Step 4: N
  n0 = sum(sum(sum(Wa.*Mx_p_My)));
  A = sum(mtimesx(Wa,Mx_m_My),3);
  A = (A-A');
  n = -[A(3,2) A(1,3) A(2,1)]';
  Np = zeros(3,3);
  skewW1 = skewSet2(reshape(Wa(:,1,:),3,numPts,1)');
  skewW2 = skewSet2(reshape(Wa(:,2,:),3,numPts,1)');
  skewW3 = skewSet2(reshape(Wa(:,3,:),3,numPts,1)');
  Np(:,1) = sum(mtimesx(skewW2,Mx_p_My(:,3,:)) - mtimesx(skewW3,Mx_p_My(:,2,:)),3);
  Np(:,2) = sum(-mtimesx(skewW1,Mx_p_My(:,3,:)) + mtimesx(skewW3,Mx_p_My(:,1,:)),3);
  Np(:,3) = sum(mtimesx(skewW1,Mx_p_My(:,2,:)) - mtimesx(skewW2,Mx_p_My(:,1,:)),3);  
  N(1,1) = n0;
  N(2:4,1) = n;
  N(1,2:4) = n;
  N(2:4,2:4) = Np;    

  % Step 5: eigenvalue of Mh
  Mh = M - c*N;
  [V,D] = eig(Mh);
  lambda = min(diag(D));
  idx = find(diag(D) == lambda);
  q = V(:,idx);
  q0 = q(1);
  ql = q(2:4);
  
  % Step 6: update c & Wa
  c = c + lambda/(q'*N*q);
  skew_ql = skew(ql);
  qM = mtimesx(skew_ql,Mx_p_My);
  T2 = mtimesx(skew_ql,qM,'t');
  qMq = permute(T2,[2,1,3]);  
  
  B = mtimesx(skew_ql,Mx_m_My);
  Bt = permute(B, [2 1 3]); % transpose in 3rd dim
  S = (B + Bt);
  Tmp = q0*q0*Mx_p_My - q0*S + qMq;
  % compute W by inversion of W^-1
  idxI = repmat(reshape(1:3*numPts,3,1,numPts),[1 3 1]);
  idxJ = repmat(reshape(1:3*numPts,1,3,numPts),[3 1 1]);
  sparseTmp = sparse(idxI(:),idxJ(:),Tmp(:));
  Wa = reshape((sparseTmp \ repmat(eye(3),numPts,1))', 3,3,numPts);  

  Rprev = R;
  R = quat2rot(q);
  
  dR = R*Rprev';
  [~, dTheta] = rot2AxisAngle(dR);
  % [dAxis, dTheta] = rot2AxisAngle(dR);
  dTheta = dTheta*180/pi;
  numIter = numIter + 1;
  
  % vTheta(numIter) = dTheta;
  % vAxis(numIter,:) = dAxis';
end

data{1} = R;
data{2} = numIter;
% data{3} = vTheta(1:numIter)';
% data{4} = vAxis(1:numIter,:);
end


function [v] = EddingtonEps(i,j,k)
  v = (j-i)*(k-j)*(k-i)/2;
end



