function [R,t,FRE,n] = anisotropic_point_register(X,Y,W,threshold)
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

% Initial estimate of the transformation assumes anisotropy:
[R,t,FRE] = point_register(X,Y);

if nargin<3 % if W not given, then give the isotropic solution
  n = 0;
  return
end

if nargin<4,threshold = 1e-6;end

n = 0; % iteration index = 0;
config_change = Inf;
Xold = R*X+repmat(t,1,N);
while (config_change>threshold)
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
  config_change = sqrt(sum(sum((Xnew-Xold).^2))/...
    sum(sum((Xold-repmat(mean(Xold,2),1,N)).^2)));
  Xold = Xnew;
end

for ii = 1:N
  D = W(:,:,ii)*(Xnew(:,ii)-Y(:,ii));
  FRE(ii) = D'*D;
end

FRE = sqrt(mean(FRE));