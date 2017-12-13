function R = quat2rot(Q0xyz)
% Note: this function follows the code provided by NDI for Polaris tracker

qnorm = norm(Q0xyz);

% %test for normalized quaternion
% slack = 0.0001;
% if (qnorm > 1 + slack) || (qnorm < 1 - slack)
%     warning('Quaternion has poor normalization: Q = [%f %f %f %f], norm = %f', Q0xyz(1), Q0xyz(2), Q0xyz(3), Q0xyz(4), norm(Q0xyz));    
% end

% renormalize quaternion
Q0xyz = Q0xyz / qnorm;

% find squared components of quaternion 
fQ0Q0 = Q0xyz(1) * Q0xyz(1);
fQxQx = Q0xyz(2) * Q0xyz(2);
fQyQy = Q0xyz(3) * Q0xyz(3);
fQzQz = Q0xyz(4) * Q0xyz(4);
% find other component multiplications of quaterion
fQ0Qx = Q0xyz(1) * Q0xyz(2);
fQ0Qy = Q0xyz(1) * Q0xyz(3);
fQ0Qz = Q0xyz(1) * Q0xyz(4);
fQxQy = Q0xyz(2) * Q0xyz(3);
fQxQz = Q0xyz(2) * Q0xyz(4);
fQyQz = Q0xyz(3) * Q0xyz(4); 

% determine the rotation matrix elements
R(1,1) = fQ0Q0 + fQxQx - fQyQy - fQzQz;
R(1,2) = 2.0 * (-fQ0Qz + fQxQy);
R(1,3) = 2.0 * (fQ0Qy + fQxQz);
R(2,1) = 2.0 * (fQ0Qz + fQxQy);
R(2,2) = fQ0Q0 - fQxQx + fQyQy - fQzQz;
R(2,3) = 2.0 * (-fQ0Qx + fQyQz);
R(3,1) = 2.0 * (-fQ0Qy + fQxQz);
R(3,2) = 2.0 * (fQ0Qx + fQyQz);
R(3,3) = fQ0Q0 - fQxQx - fQyQy + fQzQz;

end
