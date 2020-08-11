%% Function for calculating the length of each turn of spherical solenoids

%  Created by Wenshen Zhou on 7 Aug 2020

%  Introduction:
%  The code can be used for calculating the length of spherical 
%  solenoid coils. A curve in 3D space can be described with 3 parametric
%  functions x(t), y(t) and z(t). The length of the curve can be calculated
%  by integrating sqrt(x'(t)^2 + y'(t)^2 + z'(t)^2).

% Input:
% radius:          radius of the spherical solenoid
% N1:              tapering factor of the spherical solenoid
% t1:              position parameter of the starting point of the coil
% t2:              position parameter of the ending point of the coil

% Output:     
% Len_loop:        1*N array, length of the N turns of the coil

% Functions needed:
% diff_sin_helix.m

% An example of using the function:
% To calculate the length of each turn of a spherical solenoid coil 
% with a radius of 40 mm, a tapering factor of 10, and 4 turns which
% corresponding to t1 = -4*pi and t2 = 4*pi

% Len_loop = len_sin_helix(0.04, 10, -4*pi, 4*pi)

function [Len_loop] = len_sin_helix(radius, N1, t1, t2)
    N = (t2 - t1)/(2 * pi);     % number of turns
    for n = 1:N
        Len_loop(n) = integral(@(t)diff_sin_helix(radius,N1,t),t1+(n-1)*2*pi,t1+n*2*pi);  
    end
end