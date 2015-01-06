function R = quat2rot(Q0xyz)

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

% test for normalized quaternion
slack = 0.001;
if ((norm(Q0xyz) > 1 + slack) || (norm(Q0xyz) < 1 - slack))
    error('Quaternion is not normalized: Q = [%f %f %f %f], norm = %f', Q0xyz(1), Q0xyz(2), Q0xyz(3), Q0xyz(4), norm(Q0xyz));
end

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
