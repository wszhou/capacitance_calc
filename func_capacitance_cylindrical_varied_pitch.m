%% Function for calculating the capacitance (C) of cylindrical solenoids with non-uniform pitches

%  Created by Wenshen Zhou on 7 Aug 2020

%  Introduction:
%  The code can be used for calculating the capacitance of cylindrical 
%  solenoid coils with non-uniform pitches. The methods for capacitance 
%  calculation are illustrated in the paper: W. Zhou, and S. Y. Huang, 
%  "Modeling the Self-Capacitance of an Irregular Solenoid".

% Input:
% pitch:           1*N array, pitches of the N turns of the solenoid
% N:               number of turns of the solenoid
% r_w:             radius of the wire
% radius:          radius of the solenoid
% t:               thickness of the insulation coating
% epsilon_r:       permittivity of the insulation coating

% Output:     
% C:               capacitance of the solenoid


% An example of using the function:
% To calculate the self-capacitance of a 4-turn cylindrical solenoid coil 
% with non-uniform pitches of 2, 4, 6 and 8 mm for each loop, a radius of 
% 40 mm, wound with wire having a diameter of 1.024 mm. The insulation 
% layer of the wire is 60 um thick with a relative permittivity of 3.

% C = func_capacitance_cylindrical_varied_pitch([0.002 0.004 0.006 0.008], 4, 1.024e-3/2, 0.04, 60e-6, 3)


function [C] = func_capacitance_cylindrical_varied_pitch(pitch, N, r_w, radius, t, epsilon_r)
    mu_0 = 4*pi*10^-7;
    epsilon_0 = 8.854187817e-12;                    
    mu_r = 1;
    sigma = 5.96e7;                   % conductivity of copper
    d_w = 2 * r_w;                    % wire diameter
    D = 2 * radius;                   % solenoid diameter   
    din = d_w - 2 * t;                % wire diameter without insulation layer
    
    l_coil = sqrt((pi * D)^2 + pitch.^2);    % length of the wire of each turn of the solenoid
    height = sum(pitch);                     % height of the solenoid
    
    C_ins = pi/3 * epsilon_0 * epsilon_r / (log(r_w/(r_w-t)) * 2);   % capacitance in the insulation layer
    
    C_NN_total = 0;                   % initilization of C_{NN} and C_{2nd-NN}       
    C_2nd_NN_total = 0;
    
    for i = 1:(N - 1)                                           % use a loop to calculate C_{NN} for each turn
        pitch_NN(i) = 1/2 * (pitch(i) + pitch(i+1));            % average pitch
        len_NN(i) = sqrt((pi * D)^2 + pitch_NN(i)^2);           % calculate length with average pitch
        C_NN_air(i) = epsilon_0 * pi/(acosh(pitch_NN(i)/d_w)); 
        C_NN_ins(i) = C_ins;
        
        C_NN_unit = 1/(2/C_NN_ins(i) + 1/C_NN_air(i));          % unit C_{NN}
        C_NN(i) = C_NN_unit * pi * D;                           % total C_{NN}
        
        if C_NN_total == 0                                      % decide whether this is the first turn
            C_NN_total = C_NN(i);
        else
            C_NN_total = C_NN_total * C_NN(i)/(C_NN_total + C_NN(i));  % C_{NN} in series with each other
        end
    end
    
    for j = 1:(N - 2)                                                                    % use a loop to calculate C_{2nd-NN}
        pitch_2nd_NN(j) = 1/2 * ((pitch(j) + pitch(j+1)) + (pitch(j+1) + pitch(j+2)));   % average pitch
        len_2nd_NN(j) = sqrt((pi * D)^2 + pitch_2nd_NN(j)^2);                            % calculate length with average pitch
        C_2nd_NN_air(j) = epsilon_0 * pi/(acosh(pitch_2nd_NN(j)/d_w));
        C_2nd_NN_ins(j) = C_ins;
        
        C_2nd_NN_unit = 1/(2/C_2nd_NN_ins(j) + 1/C_2nd_NN_air(j));                       % unit C_{2nd-NN}
        C_2nd_NN(j) = C_2nd_NN_unit * pi * D;                                            % total C_{2nd-NN}
        
        if C_2nd_NN_total == 0                                                           % decide whether this is the first turn   
            C_2nd_NN_total = C_2nd_NN(j);
        else
            C_2nd_NN_total = C_2nd_NN_total*C_2nd_NN(j)/(C_2nd_NN_total+C_2nd_NN(j));    % C_{2nd-NN} in series with each other
        end
    end
    
    C = C_NN_total + C_2nd_NN_total;                            % total capacitance
end