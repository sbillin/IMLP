function [data] = ComputeRigidBodyXfm_Estepar(xobs,yobs,Mx,My,maxIter,maxIterR,threshAng,threshPos,isoInit,loopCovDecomp)
% Implementation of algorithm:
%
%   Estepar, et. al., "Robust Generalized Total Least Squares Iterative Closest Point
%   Registration", MICCAI 2004

numPts = size(xobs,1);

if ~exist('loopCovDecomp','var') || isempty(loopCovDecomp)
  loopCovDecomp = 0;   % default to fast method
end

%== Algorithm: Kanatani Rotation with Translation Estimation ==%

% q = zeros(4,1);
% X = zeros(3*numPts,4);
% Xt_W = zeros(4,3*numPts);
% N = zeros(4,4);

LoopData = struct(...
  'dThetas', 0,...
  'dAxis', 0,...
  'dTrans', 0,...
  'numIterR', 0 );
LoopData(maxIter).numIterR = 0; % initialize size

% Init
R = eye(3);
t = zeros(3,1);

if (isoInit)
  [R,t] = ComputeRigidBodyXfm_Iso(xobs,yobs);
end

%--- Compute t with fixed R ---%
%
%  Note: in case of a large initial offset, it is more stable
%        to estimate t first; this is a modification of
%        the method described in Estepar, et al
xrot = xobs*R';

% compute covariances
M = mtimesx(R,mtimesx(Mx,R,'t')) + My;

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

% solve t
sumWt = sum(Minv,3);
diff = (yobs - xrot)';
sumWt_xp = sum(mtimesx(Minv,reshape(diff,3,1,numPts)),3);
t = sumWt \ sumWt_xp;   % solve t = inv(sumWt) * sumWt_xp

numTotalIter = 0;
numLoopIter = 0;
dt_norm = threshPos + 1;

while ( (dt_norm > threshPos) && (numLoopIter < maxIter))
  numLoopIter = numLoopIter + 1;
  
  %--- Compute R with fixed t ---%
  
  % apply translation estimate to yi's
  ytrans = yobs - repmat(t',[numPts,1]);
  % compute new rotation estimate
  data = Compute3dRotation_Kanatani_Vectorized( xobs,ytrans,Mx,My,maxIterR,threshAng );
  R = data{1};
  numIterR = data{2};
  numTotalIter = numTotalIter + numIterR;
  LoopData(numLoopIter).numIterR = numIterR;
  % LoopData(numLoopIter).dThetas = data{3};
  % LoopData(numLoopIter).dAxis = data{4};
  
  if numIterR >= maxIterR
    % algorithm went unstable
    break;
  end  
  
  %--- Compute t with fixed R ---%
  
  % update covariances
  %  it makes most sense to do this here
  M = mtimesx(R,mtimesx(Mx,R,'t')) + My;
  
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
  
  t_prev = t;
  xrot = xobs*R';
  diff = (yobs - xrot)';
  sumWt = sum(Minv,3);
  sumWt_xp = sum(mtimesx(Minv,reshape(diff,3,1,numPts)),3);
  t = sumWt \ sumWt_xp;   % solve t = inv(sumWt) * sumWt_xp

  % update covariances
  %M = mtimesx(R,mtimesx(Mx,R,'t')) + My;
  
  dt = t - t_prev;
  dt_norm = norm(dt);
  %LoopData(numLoopIter).dTrans = dt;
end

data{1} = [R t];
data{2} = numTotalIter;
data{3} = LoopData(1:numLoopIter);
data{4} = numLoopIter;

end

