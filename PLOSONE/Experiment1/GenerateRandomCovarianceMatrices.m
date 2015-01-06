function [Cov] = GenerateRandomCovarianceMatrices( minEig, maxEig, num )
%function [Cov, Weight] = GenerateRandomCovarianceMatrices( minEig, maxEig, num )
% Generate random 3x3 covariance matrices having eigenvalues between
% minEig and maxEig
%  Cov ~ V*S*V'            (cell array)   where S = diag matrix of inverse variance
%  Weight ~ V*S^(1/2)*V'   (cell array)

Cov = cell(num,1);
%Weight = cell(num,1);

% Generate random rotations
Rot = GenerateRandomRotations(0,180,num);

eigVal = minEig*ones(num,3) + maxEig*rand(num,3);
for i=1:num
  eigM = diag(eigVal(i,:));
  Cov{i} = Rot{i}*eigM*Rot{i}';
  %Weight{i} = Rot{i}*sqrt(eigM)*Rot{i}';
end

end