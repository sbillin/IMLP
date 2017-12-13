function [ y ] = applyFrm3( F, x )
% [y] = applyTrans(F,x) applies transform F to vectors in x
% F = 4x4 homogeneous transform (rigid body or similarity xfm)
% x = nx3 vector
% y = nx3 vector

% check size
szF = size(F);
szX = size(x);
if (szF(1) ~= 4 || szF(2) ~= 4)
   error('Incorrect frame transform size');
end
if (szX(2) ~= 3)
   error('Incorrect data array size');
end
% check for homogeneous transform
if (any(F(4,:) ~= [0 0 0 1]))
   error('Non=homogeneous transform'); 
end

sR = F(1:3,1:3);
t = F(1:3,4);

% bsxfun is about 2x as fast as repmat
y = bsxfun(@plus, x*sR', t'); 
%y = x*R' + repmat(t',szX(1),1);    % for x as nx3
%y = R*x + repmat(t,1,szX(2));    % for x as 3xn

end

