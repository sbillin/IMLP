function [ Finv ] = invTrans(F, bSimilarityXfm)
% Finds the inverse of a frame transformation
%  F   ~  rigid body (without scale) or similarity (with scale) xfm
%  bSimilarityXfm  ~ boolean for similarity xfm

% check size
szF = size(F);
if (szF(1) ~= 4 || szF(2) ~= 4)
   error('Incorrect frame transform size');
end
% check for homogeneous transform
if (any(F(4,:) ~= [0 0 0 1]))
   error('Non=homogeneous transform'); 
end

if ~exist('bSimilarityXfm','var') || isempty(bSimilarityXfm)
  bSimilarityXfm = false;
end

sR = F(1:3,1:3);
t = F(1:3,4);
if ~bSimilarityXfm
  sRinv = sR';
else
  sRinv = inv(sR);
end
Finv = [sRinv, -(sRinv*t); 0 0 0 1];

%Finv = F^-1;
end