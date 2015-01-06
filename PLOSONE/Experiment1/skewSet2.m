function [ skewMatrix ] = skewSet2( V )
% Cross product operator for a set of 3D vectors
%  V ~ nx3
%  skewMatrix ~ 3 x 3 x n
n = size(V,1);
skewMatrix = zeros(3,3,n);
nZeros = zeros(n,1);

skewMatrix(1,:,1:end) = reshape([ nZeros,     -V(1:end,3),  V(1:end,2)]',1,3,n);
skewMatrix(2,:,1:end) = reshape([ V(1:end,3),  nZeros,     -V(1:end,1)]',1,3,n);
skewMatrix(3,:,1:end) = reshape([-V(1:end,2),  V(1:end,1),  nZeros]',1,3,n);

% skewMatrix(1:3:end,:) = [ nZeros,     -V(1:end,3),  V(1:end,2)];
% skewMatrix(2:3:end,:) = [ V(1:end,3),  nZeros,     -V(1:end,1)];
% skewMatrix(3:3:end,:) = [-V(1:end,2),  V(1:end,1),  nZeros];

% skewMatrix(1:3:end,:) = [...
% 0, -V(1:end,3), V(1:end,2);
% V(1:end,3), 0, -V(1:end,1);
% -V(1:end,2), V(1:end,1), 0;
% ];

% skewMatrix = [...
% 0, -V(3), V(2);
% V(3), 0, -V(1);
% -V(2), V(1), 0;
% ];

end