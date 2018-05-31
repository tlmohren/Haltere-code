function out_matrix = euler_angle( axisLetter, angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if axisLetter == 'X'
        out_matrix  = [ 1       0                     0;...
                    0,  cos( angle ),  - sin(  angle )  ; ...
                   0   sin(  angle ) cos(  angle )  ];
    elseif axisLetter == 'Y'
        out_matrix  = [cos(  angle)      0       sin( angle) ; ...
                    0              1       0 ;...
                    -sin( angle)    0        cos(  angle)]; 
    elseif axisLetter == 'Z'
        out_matrix = [cos(  angle)   sin(  angle)    0 ; ...
                     -sin(  angle) cos(  angle)     0 ;...
                     0          0               1];
    end


end

