function [Rot Ax Ang] = GenerateRandomRotations( minTheta, maxTheta, num )
% Generate num random rotations with offsets between
% minTheta and maxTheta degrees.
%
% Note: to the results repeatable, randomization seed should be set by
%       calling rand('seed',<seed_value>) prior to calling this function

minThetaRad = minTheta*pi/180;
maxThetaRad = maxTheta*pi/180;
rng = maxThetaRad - minThetaRad;
Rot = cell(num,1);

% generate num random axis on the unit sphere
Ax = rand(num,3) - 0.5*ones(num,3);
AxNorm = sqrt(sum(Ax.^2,2));
AxNorm = repmat(AxNorm,[1,3]);
Ax = Ax./AxNorm;
% generate num theta offsets
Ang = minThetaRad*ones(num,1) + rng*rand(num,1);

for i=1:num
  % generate rotation from random axis / theta
  Rot{i} = rodrigues2rot(Ax(i,:)*Ang(i));
end
  
% plot3(Ax(:,1),Ax(:,2),Ax(:,3),'.');
% axis equal

end