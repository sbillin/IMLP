function [ avgTRE, TRE, residuals, avgTREnorms, TREnorms ] = ...
  ComputeTREWithScale( Freg, scalereg, Fgt, scalegt, Xgt, Xn_gt )
% input:
%   Freg     ~ rigid body component of similarity xfm (registered)
%   scalereg ~ scale component of similarity xfm (registered)
%   Fgt      ~ rigid body component of similarity xfm (ground truth)
%   scalegt  ~ scale component of similarity xfm (ground truth)
%   Xgt     ~  n x 3   (positions)
%   Xn_gt   ~  n x 3   (normals)
% output:
%   avgTRE        ~  double
%   TRE           ~  n x 1
%   residuals     ~  n x 3
%   avgTREnorms   ~  double  (degrees)
%   TREnorms      ~  n x 1   (degrees)
% compute the TRE for each point in Xgt defined in CT space
%  TRE = ||Xgt - Xreg|| = ||Xgt - Freg * Fgt^-1 * Xgt||

% positional TRE
X = applyFrm3(invFrm3(Fgt), Xgt);
Xreg = applyFrm3(Freg, X);
residuals = Xreg - Xgt;
TRE = sqrt(sum(residuals.^2,2));    
avgTRE = mean(TRE);

if exist('Xn_gt','var')
  % rotational TRE (about a normal vector)
  nPts = size(Xn_gt,1);
  TREnorms = zeros(nPts,1);
  Xn = Xn_gt * getRot(invFrm3(Fgt))';
  Xn_reg = Xn * getRot(Freg)';
  for i = 1:nPts
    TREnorms(i) = acos( min(1, Xn_reg(i,:) * Xn(i,:)') ) * 180/pi;
  end
  avgTREnorms = mean(TREnorms);
end

end

