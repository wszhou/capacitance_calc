%% Function for calculating the capacitance (C) of spherical solenoids

%  Created by Wenshen Zhou on 7 Aug 2020

%  Introduction:
%  The code can be used for calculating the capacitance of spherical 
%  solenoid coils. The methods for capacitance calculation are illustrated
%  in the paper: W. Zhou, and S. Y. Huang, "Modeling the Self-Capacitance 
%  of an Irregular Solenoid".

% Input:
% N:               number of turns of the solenoid
% N1:              tapering factor of the spherical solenoid
% r_w:             radius of the wire
% radius:          radius of the solenoid
% s:               number of segments of the coil
% t:               thickness of the insulation coating
% epsilon_r:       permittivity of the insulation coating

% Output:     
% C:               capacitance of the solenoid

% Functions needed:
% len_sin_helix.m  
% diff_sin_helix.m

% An example of using the function:
% To calculate the self-capacitance of a 4-turn spherical solenoid coil 
% with a tapering factor of 10 , a radius of 40 mm, wound with wire having 
% a diameter of 1.024 mm. The insulation layer of the wire is 60 um thick
% with a relative permittivity of 3.

% C = func_capacitance_spherical(4, 10, 1.024e-3/2, 0.04, 200, 60e-6, 3)

function [C] = func_capacitance_spherical(N, N1, r_w, radius, s, t, epsilon_r)
    mu_0 = 4*pi*10^-7;
    epsilon_0 = 8.854187817e-12;                    
    mu_r = 1;
    sigma = 5.96e7;                   % conductivity of copper
    d_w = 2 * r_w;                    % wire diameter
    D = 2 * radius;                   % solenoid diameter   
    din = d_w - 2 * t;                % wire diameter without insulation layer
    
    n = 0;                          
    for ts = -N*pi:2*N*pi/s:N*pi
        n = n+1;
        coil(n,1) = radius*cos(ts/10)*cos(ts);          % X coordinate of point n
        coil(n,2) = radius*cos(ts/10)*sin(ts);          % Y coordinate of point n
        coil(n,3) = radius*sin(ts/10);                  % Z coordinate of point n
    end
    dl = coil(2:s+1, :) - coil(1:s, :);                 % vectors of the coil segments
    
    len_loop = len_sin_helix(radius, N1, -N*pi, N*pi);  % calculate the length of each loop with function len_sin_helix.m
    
    Sph_position = zeros(N,3);
    n2 = 0;
    for t2 = -(N-1)*pi:2*pi:(N-1)*pi                    % pick the middle segment of each loop
        n2 = n2 + 1;
        Rc_s(n2) = radius * cos(t2/N1);                 % radius of circle in xy-plane
    end
    
    for n3 = 1:(N-1)
        R_NN(n3) = 1/2 * (Rc_s(n3) + Rc_s(n3+1));               % radius of circle for C_{NN} calculation
        len_NN(n3) = 1/2 * (len_loop(n3) + len_loop(n3+1));     % average length for C_{NN} calculation
    end
    
    for n3 = 1:(N-2)
        R_2nd_NN(n3) = 1/2 * (Rc_s(n3) + Rc_s(n3+2));           % radius of circle for C_{2nd-NN} calculation
        len_2nd_NN(n3) = 1/2 * (len_loop(n3) + len_loop(n3+2)); % average length for C_{2nd-NN} calculation
    end
    
    for n4 = 1:N
        Sph_position(n4,:) = coil(1 + (n4-1) * s/N,:);          % position of the first element of each turn
    end
    Sph_position(N+1,:) = coil(s,:);
    
    C_NN_total = 0;                                                  % initilization of C_{NN} and C_{2nd-NN}
    C_2nd_NN_total = 0;
    
    C_ins = pi/3 * epsilon_0 * epsilon_r / (log(r_w/(r_w-t)) * 2);   % capacitance in the insulation layer
    
    for n4 = 1:(N-1)                                                 % calculation of C_{NN}
        pitch_NN(n4,:) = 1/2 * (norm(Sph_position(n4+1)-Sph_position(n4))+norm(Sph_position(n4+2)-Sph_position(n4+1)));  % average pitch
        C_sph_NN_air(n4,:) = epsilon_0 * pi/(acosh(pitch_NN(n4,:)/d_w));
        C_sph_NN_ins(n4,:) = C_ins;
        
        C_sph_NN_unit = 1/(2/C_sph_NN_ins(n4,:) + 1/C_sph_NN_air(n4,:));  % unit C_{NN}
        C_sph_NN(n4,:) = C_sph_NN_unit*2*pi*R_NN(n4);                     % total C_{NN}
        
        if C_NN_total == 0                                                % decide whether this is the first turn
            C_NN_total = C_sph_NN(n4,:);
        else
            C_NN_total = C_NN_total*C_sph_NN(n4,:)/(C_NN_total+C_sph_NN(n4,:));  % C_{NN} in series with each other
        end
    end
    
    for n4 = 1:(N-2)                                                      % calculation of C_{2nd-NN}
        pitch_2nd_NN(n4,:) = 1/2 * (norm(Sph_position(n4+2)-Sph_position(n4))+norm(Sph_position(n4+3)-Sph_position(n4+1))); % average pitch
        C_sph_2nd_NN_air(n4,:) = epsilon_0 * pi/(acosh(pitch_2nd_NN(n4,:)/d_w));
        C_sph_2nd_NN_ins(n4,:) = C_ins;
        
        C_sph_2nd_NN_unit = 1/(2/C_sph_2nd_NN_ins(n4,:) + 1/C_sph_2nd_NN_air(n4,:));  % unit C_{2nd-NN}
        C_sph_2nd_NN(n4,:) = C_sph_2nd_NN_unit*2*pi*R_2nd_NN(n4);                     % total C_{2nd-NN}
        
        if C_2nd_NN_total == 0                                            % decide whether this is the first turn
            C_2nd_NN_total = C_sph_2nd_NN(n4,:);
        else
            C_2nd_NN_total = C_2nd_NN_total*C_sph_2nd_NN(n4,:)/(C_2nd_NN_total+C_sph_2nd_NN(n4,:));  % C_{2nd-NN} in series with each other
        end
    end
    
    C = C_2nd_NN_total + C_NN_total;                                      % total capacitance
end

