function [ avgTRE, TRE, residuals, avgTREnorms, TREnorms ] = ...
  ComputeTRE( Freg, Fgt, Xgt, Xn_gt, bSimilarityXfm )
% input:
%   Xgt     ~  n x 3   (positions)
%   Xn_gt   ~  n x 3   (normals)
%   bSimilarityXfm  ~ boolean whether Freg & Fgt are similarity xfms
% output:
%   avgTRE        ~  double
%   TRE           ~  n x 1
%   residuals     ~  n x 3
%   avgTREnorms   ~  double  (degrees)
%   TREnorms      ~  n x 1   (degrees)
% compute the TRE for each point in Xgt defined in CT space
%  TRE = ||Xgt - Xreg|| = ||Xgt - Freg * Fgt^-1 * Xgt||

if ~exist('bSimilarityXfm','var') || isempty(bSimilarityXfm)
  bSimilarityXfm = false;
end

% positional TRE
Xreg = applyFrm3(Freg * invFrm3(Fgt, bSimilarityXfm), Xgt);
residuals = Xreg - Xgt;
TRE = sqrt(sum(residuals.^2,2));
avgTRE = mean(TRE);

if exist('Xn_gt','var') && ~isempty(Xn_gt)
  % rotational TRE (about a normal vector)
  nPts = size(Xn_gt,1);
  TREnorms = zeros(nPts,1);
  if bSimilarityXfm
    [Rgtinv,~] = getRotScale(invFrm3(Fgt));
    [Rreg,~] = getRotScale(Freg);
  else
    Rgtinv = getRot(invFrm3(Fgt));
    Rreg = getRot(Freg);
  end  
  Xn_reg = Xn_gt * (Rreg*Rgtinv)';
  for i = 1:nPts
    TREnorms(i) = acos( min(1, Xn_reg(i,:) * Xn_gt(i,:)') ) * 180/pi;
  end
  avgTREnorms = mean(TREnorms);
end

end

