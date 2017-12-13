function [Tctr] = MeshGetTriangleCenters( V, T )
% Computes neighbors for each triangle
%
% inputs:
%   V - vertices                   (numV x 3)
%   T - triangles                  (numT x 3)
%  
% outputs:
%   Tctr - triangle center points  (numT x 3)
%

disp('Computing Triangle Center Points')

numT = size(T,1);
Tctr = zeros(numT,3);

for Tx = 1:numT
  Vi = V(T(Tx,:),:);
  Tctr(Tx,:) = mean(Vi);  
end

end
