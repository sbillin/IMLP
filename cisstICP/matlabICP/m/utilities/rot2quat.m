function [ Q0xyz ] = rot2quat(Rot)

% Converts a 3D rotation matrix to quaternion form
% This code is based on C++ code from the CISST library with minor 
%  modifications, as well as corrections for errors!
%   cisstVector / vctQuaternionRotation3Base.cpp /
%   vctQuaternionRotation3BaseFromRaw()

% Note: for CISST library, the default tolerance for a float is set to 
%       1.0e-5f and for a double it is set to 1.0e-9.
tolerance = 1e-9;

trace = Rot(1,1) + Rot(2,2) + Rot(3,3);

if (trace > tolerance)
    S = sqrt(trace + 1) * 2;    % S = 4*q0
    Q0xyz(1) = 0.25 * S;
    Q0xyz(2) = (Rot(3, 2) - Rot(2, 3)) / S;
    Q0xyz(3) = (Rot(1, 3) - Rot(3, 1)) / S;
    Q0xyz(4) = (Rot(2, 1) - Rot(1, 2)) / S;
    
    disp('case 0');    
    
elseif (Rot(1, 1) > Rot(2, 2) && Rot(1, 1) > Rot(3, 3))

    disp('case 1');
    
    S = sqrt(1 + Rot(1, 1) - Rot(2, 2) - Rot(3, 3)) * 2;   % S = 4*qx
    
    %Q0xyz(1) = (Rot(2, 3) - Rot(3, 2)) / S;     % CISST BUG!!!
    Q0xyz(1) = (Rot(3, 2) - Rot(2, 3)) / S;  
    Q0xyz(2) = 0.25 * S;
    Q0xyz(3) = (Rot(1, 2) + Rot(2, 1)) / S;
    Q0xyz(4) = (Rot(3, 1) + Rot(1, 3)) / S;
    
elseif (Rot(2, 2) > Rot(3, 3))
    
    disp('case 2');    
    
    S = sqrt(1 + Rot(2, 2) - Rot(1, 1) - Rot(3, 3)) * 2;    % S = 4*qy

    %Q0xyz(1) = (Rot(3, 1) - Rot(1, 3)) / S;    % CISST BUG!!!
    Q0xyz(1) = (Rot(1, 3) - Rot(3, 1)) / S;
    Q0xyz(2) = (Rot(1, 2) + Rot(2, 1)) / S;
    Q0xyz(3) = 0.25 * S;
    Q0xyz(4) = (Rot(2, 3) + Rot(3, 2)) / S;

else

    disp('case 3');    
    
    S = sqrt(1 + Rot(3, 3) - Rot(1, 1) - Rot(2, 2)) * 2;    % S = 4*qz
    
    %Q0xyz(1) = (Rot(1, 2) - Rot(2, 1)) / S;    % CISST BUG!!!
    Q0xyz(1) = (Rot(2, 1) - Rot(1, 2)) / S;  
    Q0xyz(2) = (Rot(3, 1) + Rot(1, 3)) / S;
    Q0xyz(3) = (Rot(2, 3) + Rot(3, 2)) / S;
    Q0xyz(4) = 0.25 * S;

end

end
