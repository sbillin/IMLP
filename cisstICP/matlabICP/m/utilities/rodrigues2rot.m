function [ R ] = rodrigues2rot(v)
     
% if rotation is less than a small tolerance, then return identity
%  to prevent NaN's in calculation
if (norm(v) < eps)
    R = eye(3);
else
    I = eye(3);
    theta = norm(v);
    k = v/norm(v);
    kSkew = [...
        0 -k(3) k(2)
        k(3) 0 -k(1)
        -k(2) k(1) 0];

    R = I + kSkew*sin(theta) + (1-cos(theta))*(kSkew*kSkew);
    %R = I*cos(theta) + kSkew*sin(theta) + (1-cos(theta))*(k*k')
end

end