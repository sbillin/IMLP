function [data] = Compute3dRotation_TotalLeastSquares(xobs,yobs,Mxi,Myi,maxIter,threshAng)


numPts = size(xobs,1);

%== Algorithm: Total Least Squares ==%

R = eye(3);
I3 = eye(3);
J = zeros(3*numPts,3);
Jt_Minv = zeros(3,3*numPts);


%xobs0 = xobs;
%yobs0 = yobs;

% Step 1: initialize
alpha = zeros(3,1);
R = eye(3);
%xcalc = xobs;
%ycalc = yobs;

numIter = 0;
dTheta = threshAng + 1;

%maxIter = 50;
vThetas = zeros(maxIter,1);
vAxis = zeros(maxIter,3);

while (dTheta > threshAng && numIter < maxIter)
  %tic

  % Step 2: F0
  Rxobs = xobs*R';
  Rxobs_t = Rxobs';
  yobs_t = yobs';
  F0 = yobs_t(:) - Rxobs_t(:);  % stack into a single vector [f1x f1y f1z f2x...]'  
  %---------------
  %yobs_t = yobs';
  %xobs_t = xobs';
  %F0 = yobs_t(:) - xobs_t(:);  % stack into a single vector [y1x y1y y1z y2x...]'
  
  % Step 3: J & Fx
  %   Note: Fx is block diagonal with all sub-blocks being equal
  %         to each other (equal to Fxi)
  J = skewSet(Rxobs);
  Fxi = -R;
  %Fxi = -rodrigues2rot(alpha);
  %---------------
  %J = skewSet(xobs);
  %J(1:3:end,:) = skew(xobs(1:end,:));
  %Fxi = -(I3 + skew(alpha));  % identity when alpha = 0
  
  % Step 4: solve dAlpha by least squares
  Mi = Fxi*Mxi*Fxi' + Myi;
  Minv = inv(Mi);
  for i=1:numPts
      % Minv is block diagonal => only multiply the sub-blocks
      ix = (i-1)*3 + 1;
      Jt_Minv(:,ix:ix+2) = J(ix:ix+2,:)'*Minv;
  end
  Jt_Minv_J = Jt_Minv * J;
  nJt_Minv_F0 = Jt_Minv * (-F0);
  % Ax = b
  %  Note: most of the algorithm run-time is for this SVD call
  [U S V] = svd(Jt_Minv_J, 'econ');
  Sinv = diag(1./diag(S));
  dAlpha = V*Sinv*U'*nJt_Minv_F0;

  % Step 5: R
  dR = rodrigues2rot(dAlpha);
  Rprev = R;
  R = dR*R;

%   % This step is unneccessary
%   % Step 6: xcalc
%   Wi = R*Mxi*R' + Myi;
%   Winv = inv(Wi);
%   Mx_Rt_Winv = Mxi*R'*Winv;
%   res = yobs - xobs*R';
%   xcalc = xobs + res * Mx_Rt_Winv';
%   ycalc = xcalc * dR';

%   % Plotting:
%   %  find calculated xi relative to their initial location
%   %  (rotate xcalc back to starting position)
%   xcalc0 = xcalc * Rprev;
%   plot3(xact(:,1),xact(:,2),xact(:,3), 'ob');
%   hold on
%   plot3(xobs0(:,1),xobs0(:,2),xobs0(:,3), '.b');
%   plot3(yact(:,1),yact(:,2),yact(:,3), 'or');
%   plot3(yobs0(:,1),yobs0(:,2),yobs0(:,3), '.r');
%   plot3(xcalc0(:,1),xcalc0(:,2),xcalc0(:,3), '.k');
%   plot3(ycalc(:,1),ycalc(:,2),ycalc(:,3), '.g');
%   xlabel('X')
%   ylabel('Y')
%   zlabel('Z')
%   view([0 90])
%   axis equal
%   hold off
%
%   msg = ['Click for Next Iteration'];
%   h = msgbox(msg);
%   uiwait(h);

  %--------------
  % Step 7: xobs
  %xobs = xobs0*R';
  %alpha = zeros(3,1);
  
  [dAxis dTheta] = rot2AxisAngle(dR);
  dTheta = dTheta*180/pi;
  numIter = numIter + 1;
  
  vTheta(numIter) = dTheta;
  vAxis(numIter,:) = dAxis';
  %toc
end 

data{1} = R;
data{2} = numIter;
data{3} = vTheta(1:numIter)';
data{4} = vAxis(1:numIter,:);
end







