%% Function for calculating the capacitance (C) of cylindrical solenoids with uniform pitches

%  Created by Wenshen Zhou on 7 Aug 2020

%  Introduction:
%  The code can be used for calculating the capacitance of cylindrical 
%  solenoid coils with uniform pitches. The methods for capacitance 
%  calculation are illustrated in the paper: W. Zhou, and S. Y. Huang, 
%  "Modeling the Self-Capacitance of an Irregular Solenoid".

% Input:
% pitch:           pitch of the solenoid
% N:               number of turns of the solenoid
% r_w:             radius of the wire
% radius:          radius of the solenoid
% t:               thickness of the insulation coating
% epsilon_r:       permittivity of the insulation coating

% Output:     
% C:               capacitance of the solenoid


% An example of using the function:
% To calculate the self-capacitance of a cylindrical solenoid coil with a
% pitch of 15 mm, 4 turns, a radius of 40 mm, wound with wire having a 
% diameter of 1.024 mm. The insulation layer of the wire is 60 um thick
% with a relative permittivity of 3.

% C = func_capacitance_cylindrical(0.015, 4, 1.024e-3/2, 0.04, 60e-6, 3)

function [C] = func_capacitance_cylindrical(pitch, N, r_w, radius, t, epsilon_r)
    mu_0 = 4*pi*10^-7;
    epsilon_0 = 8.854187817e-12;                    
    mu_r = 1;
    sigma = 5.96e7;                   % conductivity of copper
    d_w = 2 * r_w;                    % wire diameter
    D = 2 * radius;                   % solenoid diameter   
    din = d_w - 2 * t;                % wire diameter without insulation layer
    
    l_t = sqrt((pi * D)^2 + pitch^2);            % length of one turn
    l_t_2 = sqrt((pi * D)^2 + (2 * pitch)^2);    % length of two turns
    
    C_ins = pi/3 * epsilon_0 * epsilon_r / (log(r_w/(r_w-t)) * 2);   % capacitance in the insulation layer
    
    C_1st_air = epsilon_0 * pi/(acosh(pitch/d_w));                   % capacitance in the air for C_{NN}
    C_1st_unit = 1/(2/C_ins + 1/C_1st_air);                          % C_{NN} in unit length
    C_1st = C_1st_unit * l_t;                                        % total C_{NN}
    
    C_2nd_air = epsilon_0 * pi/(acosh(2 * pitch/d_w));               % capacitance in the air for C_{2nd-NN}
    C_2nd_unit = 1/(2/C_ins + 1/C_2nd_air);                          % C_{2nd-NN} in unit length
    C_2nd = C_2nd_unit * l_t_2;                                      % total C_{2nd-NN}
    
    C = C_1st/(N - 1) + C_2nd/(N-2);                                 % total C
end