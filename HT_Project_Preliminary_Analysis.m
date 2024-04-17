clear; clear all;
% Set up problem
%air_int.T = 21.1; % Initial inside air temperature (C)
h_e_inf = 50; % Heat transfer coefficient on the exterior side (W/m^2K)
concrete.k = 2.25; % Thermal conductivity of concrete (W/m-K)
concrete.t = 0.1524;% Thickness of concrete part of wall (m)
pcm.k = 0.12; % Thermal conductivity of PCM (W/m-K)
pcm.t = 0.0254; % Thickness of PCM layer (m)
A_wall = 96.39*10; % Surface area of wall (m^2)
q_gen = 100000000; % Due to servers (W)

time_steps =3600;  % Total number of time steps
dt = 1;% Time step (s)

concrete.R = concrete.t/(A_wall*concrete.k);
pcm.R = pcm.t/(A_wall*pcm.k);
air_int.h = 20; % W/m^2-K estimated through common heat transfer coefficients
air_ext.h = 10; % W/m^2-K
air_int.R = 1/(A_wall*air_int.h);
air_ext.R = 1/(A_wall*air_ext.h);
R_wall = concrete.t/(concrete.k*A_wall);


time = 0:dt:time_steps;
T_wall_int = zeros(length(time),1);
T_air_int = zeros(length(time),1);
T_air_int(1) = 21.1;
air_ext.T = linspace(26,32,length(time));     % Outside air temperature (C)
 % External wall temperature (C). This is another assumed value to keep 

 T_wall_ext = zeros(length(time),1);
 q_conv_ext = zeros(length(time),1);
 q_conv_i = zeros(length(time),1);
 q_cond_wall = zeros(length(time),1);
%iterate
for i = 2:length(time)
    T_wall_ext(i) = 28;
     %q_gen/(4*A_wall); %Internal convective heat transfer rate
    q_conv_ext(i) = (T_wall_ext(i)-air_ext.T(i))/air_ext.R; % external convective heat transfer rate
    q_cond_wall(i) = q_conv_ext(i);
    T_wall_int(i) = -q_cond_wall(i)*R_wall+T_wall_ext(i); % Heat conduction through wall
    q_conv_i(i)= air_int.R*(T_wall_int(i)-T_air_int(i-1));

    %T_wall_int(i) = air_int.R*q_conv_i(i)+T_air_int(i);
    T_air_int(i) = T_wall_int(i)-(q_conv_i(i)*dt)*air_int.R; %This term should be making the internal temperature change
    
    %T_wall_int(i+1) = T_wall_int(i)+dt*(q_conv_i(i)-q_cond_wall(i));
end

% Plot
figure;
plot(time, T_air_int, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Inside Air Temperature (°C)');
title('Inside Air Temperature vs. Time');
grid on;