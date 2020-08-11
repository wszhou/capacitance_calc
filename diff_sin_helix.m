%% Function to be used by len_sin_helix.m for calculating the length of spherical solenoids

%  Created by Wenshen Zhou on 7 Aug 2020

%  Introduction:
%  The code can be used by len_sin_helix.m for calculating the length of 
%  spherical solenoids


% Input:
% radius:          radius of the spherical solenoid
% N1:              tapering factor of the spherical solenoid
% t:               position parameter

% Output:     
% y:               sqrt(x'(t)^2 + y'(t)^2 + z'(t)^2)

% An example of using the function:
% To calculate the length of each turn of a spherical solenoid coil 
% with a radius of 40 mm, a tapering factor of 10, and 4 turns which
% corresponding to t1 = -4*pi and t2 = 4*pi

% Len_loop = len_sin_helix(0.04, 10, -4*pi, 4*pi)

function y = diff_sin_helix(radius, N1, t)
    x_d = -radius.*(1/N1 * sin(t/N1).*cos(t)+cos(t/N1).*sin(t));  % derivative of x(t)
    y_d = radius.*(cos(t/N1).*cos(t)-1/N1*sin(t/N1).*sin(t));     % derivative of y(t)
    z_d = radius * 1/N1 * cos(t/N1);                              % derivative of z(t)
    y = sqrt(x_d.^2 + y_d.^2 + z_d.^2);                               
end
