function [R,dR] = dR_Rodrigues( a )
% computes the rotation matrix described by Rodrigues vector "a" 
%  and the partial derivatives of the rotation matrix wrt each 
%  element of "a"
%
% Note: each partial derivatives is represented as a seperate 
%       3x3 matrix storing the change in each element of the rotation
%       matrix with respect to one of the elements of "a"
%

theta = norm(a);
if (theta <= eps)  % don't divide by zero
  alpha = [0,0,0]; 
  theta = 0;
  %disp('==== WARNING: tiny theta ====')
else
  alpha = a/theta;
end

%% precompute
gTheta = alpha;
theta2 = theta*theta;
theta3 = theta2*theta;
aaT = a*a';

I = eye(3);
sk_alpha = skew(alpha);
AlphaAlphaT_I = alpha*alpha' - I;
sTheta = sin(theta);
cTheta = cos(theta);

%% compute current rotation matrix
R = I + sTheta*sk_alpha + (1-cTheta)*AlphaAlphaT_I;

%% compute rotation differential wrt to each Rodrigues parameter
dR = cell(1,3);
dai = {[1,0,0]',[0,1,0]',[0,0,1]'}; % wrt dax, day, daz
for i=1:3
  da = dai{i};
  if (theta == 0)
    % set alpha in direction of da
    sk_alpha = skew(da);
    dR{i} = sk_alpha;
  else    
    dTheta = gTheta'*da;
    % TODO: can speed up each line below
    daat_adat = da*a' + a*da';    
    dR{i} = cTheta*dTheta*sk_alpha ...
               + (sTheta/theta2)*skew(theta*da - a*dTheta) ...
               + sTheta*dTheta*AlphaAlphaT_I ...
               + ((1-cTheta)/theta3)*(theta*(daat_adat) - 2*dTheta*aaT);  
  end
end

end