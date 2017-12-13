function [ hMesh ] = PlotMesh( V,T, faceAlpha,edgeAlpha,color,color_map, color_range, edgeColor )

if ~exist('faceAlpha','var') || isempty(faceAlpha)
  faceAlpha = 1.0;
end
if ~exist('edgeAlpha','var') || isempty(edgeAlpha)
  edgeAlpha = 0.1;
end
if ~exist('color','var') || isempty(color)
  color = 0.85;
end
if ~exist('color_map','var') || isempty(color_map)
  color_map = 'Bone';
end
if ~exist('color_range','var') || isempty(color_range)
  color_range = [0,1];
end
if ~exist('edgeColor','var') || isempty(edgeColor)
  edgeColor = [0 0 0];
end

numT = size(T,1);
colormap(color_map)         % set colormap
if isscalar(color)
  C = ones(numT,1)*color;     % set mesh color
else
  C = color;
end

% plot mesh
hMesh = trisurf(T, V(:,1),V(:,2),V(:,3),C);
set(hMesh, 'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha, 'FaceAlpha', faceAlpha ) %'FaceColor', [1 0 0], );
caxis(color_range)            % set colormap range (must be done after trimesh)
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal

end