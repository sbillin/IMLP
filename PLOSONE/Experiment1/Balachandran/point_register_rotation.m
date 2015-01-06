function [R] = point_register_rotation(X, Y)

[Ncoords Npoints] = size(X);
[Ncoords_Y Npoints_Y] = size(Y);

if Ncoords ~= 3 | Ncoords_Y ~= 3
  error('Each argument must have exactly three rows.')
elseif (Ncoords ~= Ncoords_Y) | (Npoints ~= Npoints_Y)
  error('X and Y must have the same number of columns.');
elseif Npoints < 3
  error('X and Y must each have 3 or more columns.');
end

% Xbar = mean(X,2); % X centroid
% Ybar = mean(Y,2); % Y centroid
% Xtilde = X-repmat(Xbar,1,Npoints); % X relative to centroid
% Ytilde = Y-repmat(Ybar,1,Npoints); % Y relative to centroid
% H = Xtilde*Ytilde'; % cross covariance matrix

H = X*Y';
[U S V] = svd(H); % U*S*V' = H
R = V*diag([1, 1, det(V*U)])*U';
%t = Ybar - R*Xbar;
% SDB: leave this out for accurate run-time
% FREvect = R*X + repmat(t,1,Npoints) - Y;
% FRE = sqrt(mean(sum(FREvect.^2,1)));
end