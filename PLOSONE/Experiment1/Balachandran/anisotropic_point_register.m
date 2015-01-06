function [R,t,n] = anisotropic_point_register(X,Y,M,maxIter,threshAng,threshPos,isoInit)
%
% SDB: modifications from implementation by Balachandran (marked in code by
% initials SDB)
% 
%  - FRE removed (to improve efficienty)
%  - isotropic initialization optional
%  - changed termination condition threshold on change in rotation and translation
%    plus a maximum iteration count
%  - using covariance (M) of y point set as input rather than M^-1/2
%
%function [R,t,FRE,n] = anisotropic_point_register(X,Y,W,threshold)
% X is the moving set, which is registered to the static set Y. Both are 3
% by N, where N is the number of fiducials. W is a 3-by-3-by-N array, with
% each page containing the weighting matrix for the Nth pair of points.
% THRESHOLD is the size of the change to the moving set above which the
% iteration continues.
%
% Creation:
% R. Balachandran and J. M. Fitzpatrick
% December 2008
if nargin<3,error('X and Y must be given as input');end
if size(X,1)~=3 && size(Y,1)~=3,error('X and Y must be 3 by N.');end
N = size(X,2);
if size(Y,2)~=N,error('X and Y must have the same number of points.');end

% SDB: compute W
%  I don't know any way to compute Mi^(1/2) for i=1..N without using
%  a loop to operate on each Mi individually
sqrtDinv = zeros(3,3);
for i=1:N
  [V,D] = eig(M(:,:,i));
  sqrtDinv(1:4:end) = sqrt(1./diag(D)); % set diagonal values
  W(:,:,i) = V*sqrtDinv*V'; % W = M^(-1/2)
  %W(:,:,i) = sqrtm(M(:,:,i));  % this is slower
end  

% Initial estimate of the transformation assumes anisotropy: 
%[R,t,FRE] = point_register(X,Y);   % SDB

% SDB
R = eye(3);
t = zeros(3,1);
if (isoInit)
  % Initial estimate of the transformation assumes anisotropy:
  [R,t] = point_register(X,Y);      % SDB
end

if nargin<3 % if W not given, then give the isotropic solution
  n = 0;
  return
end

%if nargin<4,threshold = 1e-6;end   % SDB
n = 0; % iteration index = 0;
%config_change = Inf;   % SDB
Xold = R*X+repmat(t,1,N);

% SDB
dTheta = 1;
dt_norm = 1;

while ( ((dTheta > threshAng) | (dt_norm > threshPos)) & (n < maxIter)) % SDB
% while (config_change>threshold) % SDB
  n = n+1;
    
  C = C_maker(Xold,W);
  e = e_maker(Xold,Y,W);
  q = C\e;
  if n > 1,q = (q + oldq)/2; end %damps oscillations
  oldq = q;
  delta_t = [q(4) q(5) q(6)]';
  delta_theta = [1 -q(3) q(2); q(3) 1 -q(1); -q(2) q(1) 1];
  [U,Lambda,V] = svd(delta_theta);
  delta_R = U*V';
  R = delta_R*R; % update rotation
  t = delta_R*t+delta_t; % update translation
  Xnew = R*X+repmat(t,1,N); % update moving points
  % SDB: no longer need this
  % config_change = sqrt(sum(sum((Xnew-Xold).^2))/...
  %   sum(sum((Xold-repmat(mean(Xold,2),1,N)).^2)));
  Xold = Xnew;
  
  % SDB
  [~, dTheta] = rot2AxisAngle(delta_R);
  dTheta = dTheta*180/pi;
  dt_norm = norm(delta_t);  
end

% SDB: do not compute FRE here since it may impact the
%      algorithm run-time
% for ii = 1:N
%   D = W(:,:,ii)*(Xnew(:,ii)-Y(:,ii));
%   FRE(ii) = D'*D;
% end
% FRE = sqrt(mean(FRE));

end
