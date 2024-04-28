function [rink, roof, amb] = MEGN571_Example6_13_S2024()
%
% Adapted from Bergman and Lavine 7: Ice hockey rinks must have ceilings with 
% high reflectivity, or condensation may occur on the ceiling and cause 
% water to drip on the ice. Condensation occurs when the surface of the 
% ceiling drops below the dew point of the air inside the rink.  Imagine 
% the rink is a cylindrical structure with a diameter D = 50 m and a 
% height L = 10 m.  The temperature of the ice and sidewalls are 
% T_ice = -5°C and Twall = 15°C respectively.  The indoor air is 
% equilibrated with the walls at T_air = 15°C and has a relative humidity 
% phi_air = 0.70. The heat transfer coefficient between the air and the 
% ceiling is estimated to be h_conv = 5 W/^2.K.  The ceiling has a 
% thickness of t_ceiling = 0.3 m and a very low average thermal 
% conductivity k_ceiling = 0.035 W/m.K due to insulation. Outdoor air 
% is at T_outside = -5°C. The ice and sidewalls are approximated as 
% blackbodies, and the ceiling is assumed to be a diffuse grey surface.
% .
% a)	Calculate T_ceiling for different eps_ceiling and determine 
%       if condensation will occur.
% b)	Determine what the effects of ceiling (insulation) thickness 
%       is for a given emissivity.
%
%... Universal constants
sigma = 5.67e-8;  % Stefan-Boltzmann constant in W/m^2*K^4
%
%... Set number of surfaces in enclosure calculation
rink.n_surf = 3;
roof.n_surf = 1;
%
% Geometric parameters for rink
rink.y = 10;    % rink height [m]
rink.D = 50;   % rink diameter [m]
roof.t = 0.30;    % rink ceiling thickness [m]
%
%... Calculate areas 1- ice,  2- walls, 3 - ceiling, 4 - roof
rink.A = zeros(rink.n_surf,1);
rink.A(1) = 0.25*pi*rink.D^2;  % ice surface area [m^2]
rink.A(2) = pi*rink.D*rink.y;  % cylindrical wall surface area [m^2]
rink.A(3) = rink.A(1);  % internal ceiling surface area [m^2]
roof.A(1) = rink.A(3);  % external roof surface area [m^2]
%
%...Set conductivity and internal and external heat transfer cofficients 
%   for ceiling 
roof.k_cond_c = 0.035;  % roof thermal conductivity [W/m*K]
rink.h_conv(3) = 5.0; % ceiling heat transfer coefficient [W/m^2*K]
roof.h_conv(1) = 500.0; %  roof heat transfer coefficient [W/m^2*K]
%
%... Set emissivities of enclosure surfaces:
%    1: ice, 2: side walls, 3: internal ceiling,
rink.eps = zeros(rink.n_surf,1);
rink.eps(1) = 0.8;
rink.eps(2) = 1.0;
rink.eps(3) = 0.05;   % reflective ceiling emissivity
% rink.eps(3) = 0.94;  % painted ceiling emissivity
roof.eps(1) = 0.05;    % external roof emissivity [--]
%
%... Set known temperatures and heat fluxes
rink.T = zeros(rink.n_surf,1);
rink.q = zeros(rink.n_surf,1);  
rink.T(1) = 0 + 273.15;     % approximate external ceiling temperature guess [K]
rink.T_int  = 15 + 273.15;    % internal rink air temperature [K]
rink.T(2) = rink.T_int;    % internal rink wall temperature [K]
T_3_guess = 5 + 273.15;     % approximate internal ceiling temperature guess [K]
roof.T(1) = -5 + 273.15;    % external roof temperature [K]
amb.T = -5 + 273.15;   % ambient temperature [K]
rink.phi = 0.70;    
rink.T_dew = 9.4 + 273.15;   % dew point temperature in rink [K]
%
%...Determine radiative heat transfer externally and combine resistance to
%   convective resistance and add to conductive resistance
roof.R_t_ext = roof.t/(roof.A(1)*roof.k_cond_c); % [K/W]
roof.h_rad(1) = roof.eps(1)*sigma*(roof.T(1)^2 + amb.T^2)...
                *(roof.T(1) + amb.T);
if roof.h_conv(1) + roof.h_rad(1) > 0
    roof.R_t_ext = roof.R_t_ext + ...
        1/(roof.A(1)*(roof.h_conv(1) + roof.h_rad(1)));  % [K/W]
end
%
%% Set up flag for determining calculations for each enclosure surface type
%      surf.flag = 1  -- blackbody
%      surf.flag = 2  -- known surface temperature
%      surf.flag = 3  -- known surface heat transfer
rink.flag = [2;  1;  2];   
%
%% Find View Factors for cylindrical ice rink
rink.F_ij(1,1) = 0.0;
rink.F_ij(3,3) = 0.0;
rink.F_ij(1,3) =  F_ij_ParallelDisks(0.5*rink.D,0.5*rink.D,rink.y);
rink.F_ij(1,2) = 1 - rink.F_ij(1,1) - rink.F_ij(1,3);
rink.F_ij(2,1) = rink.F_ij(1,2)*rink.A(1)/rink.A(2);
rink.F_ij(2,3) = rink.F_ij(2,1);
rink.F_ij(2,2) = 1 - rink.F_ij(2,1) - rink.F_ij(2,3);
rink.F_ij(3,1) = rink.F_ij(1,3)*rink.A(1)/rink.A(3);
rink.F_ij(3,2) = 1 - rink.F_ij(3,1) - rink.F_ij(3,3);
%
for i = 1:rink.n_surf
    rink.A_iF_ij(i,:) = rink.F_ij(i,:)*rink.A(i);
end
%
%...Set up steadyenergy balance on internal surface of ceiling
%   0 = q_rad_int + q_conv_int + q_ext
T_3 = fzero(@(T_3) ceiling_qcalc(T_3, rink, roof, amb, sigma), T_3_guess);
rink.T(3) = T_3;
end

%   
%...Calculate the energy flux balance on the ceiling by soliving for T(3).
function Res = ceiling_qcalc(T_3, rink, roof, amb, sigma)
    
T = [rink.T(1); rink.T(2); T_3];

q_ext = (T_3 - amb.T)/roof.R_t_ext;
q_conv_int = rink.h_conv(3)*rink.A(3)*(T_3 - rink.T(1));
% q_rad_int = rink.A(3)*rink.eps(3)*sigma*T_3^4 - ...
%        rink.AF_ij(1,3)*rink.eps(3)*sigma*rink.T(1)^4 - ...
%        rink.AF_ij(2,3)*rink.eps(3)*sigma*rink.T(2)^4;
% Res = - q_ext - q_conv_int - q_rad_int;

[J_rad, T, q] = Jrad_calc(rink, rink.A_iF_ij, sigma)

q_rad_int = rink.A(3)*rink.eps(3)/(1-rink.eps(3))*(sigma*T(3)^4 - J_rad(3));
Res = - q_ext - q_conv_int - q_rad_int;
end
%--------------------------------------------------------------
%%  This function solves for the radiosities, temperatures, and heat fluxes
%    for a radiation enclosure without any convective heat transfer.  It
%    solves the equation for the gray surface exchange radiosities  (to be completed during in-class discussion)
function [J_rad, T, q] = Jrad_calc(surf, A_iF_ij, sigma)
%
%  Put known surface quantitiates into ouputs.
T = surf.T;
q = surf.q;
%...Set up the matrix that solves for the different radiosities J's with the
%   equation  A_rad * J_rad = b_rad
A_rad = zeros(surf.n_surf,surf.n_surf);
b_rad = zeros(surf.n_surf,1);
%
%...Determine appropirate equation for matrix A_rad and vector b_rad
%     if flag = 1, the surface is black
%     if flag = 2, the surface has a known temperature
%     if flag >= 3, the surface has a known heat flux
for i_surf = 1:surf.n_surf
   if surf.flag(i_surf) == 1     
      A_rad(i_surf,i_surf) = 1.0;
      b_rad(i_surf) = sigma*surf.T(i_surf)^4;   
%
   elseif surf.flag(i_surf) == 2  
      A_rad(i_surf,:) = - A_iF_ij(i_surf,:);
      A_rad(i_surf,i_surf) = surf.eps(i_surf)*surf.A(i_surf)/(1 - surf.eps(i_surf)) + sum(A_iF_ij(i_surf,:)) - A_iF_ij(i_surf,i_surf);
      b_rad(i_surf) = surf.eps(i_surf)*surf.A(i_surf)/(1 - surf.eps(i_surf))*sigma*surf.T(i_surf)^4;   
%
   elseif surf.flag(i_surf) >= 3   % 
      A_rad(i_surf,:) = - A_iF_ij(i_surf,:);
      A_rad(i_surf,i_surf) = sum(A_iF_ij(i_surf,:)) - A_iF_ij(i_surf,i_surf);
      b_rad(i_surf) = surf.q(i_surf);
   end
end
%
%...Solve for J's with matrix inversion of A_rad 
J_rad = A_rad\b_rad;
%
%...Find heat fluxes leaving each surface for those for which q is not 
%   specified and temperatures for those surfaces which Temp is not
%   specified
for i_surf = 1:surf.n_surf
%
   if surf.flag(i_surf) == 1  
     q(i_surf) = A_iF_ij(i_surf,:)*(J_rad(i_surf) - J_rad);
%
   elseif surf.flag(i_surf) == 2  % 
      q(i_surf) = surf.eps(i_surf)*surf.A(i_surf)/(1 - surf.eps(i_surf))*(sigma*surf.T(i_surf)^4 - J_rad(i_surf));
%
   elseif surf.flag(i_surf) >= 3  % 
       T(i_surf) = (((1 - surf.eps(i_surf))/(surf.eps(i_surf)*surf.A(i_surf))*surf.q(i_surf) ...
                          + J_rad(i_surf))/sigma)^0.25;
   end
end  
end
%--------------------------------------------------------------
%
%% This function calculates a view factor for parallel disks based on correlation
%   in textbook where L is the separation disk and r_i and r_j are the
function F_ij= F_ij_ParallelDisks(r_i,r_j,L_y)
    R_i = r_i/L_y;
    R_j = r_j/L_y;
    S = 1 + (1 + R_j^2)/R_i^2;
    F_ij = 0.5*(S - (S^2 - 4*(r_j/r_i)^2)^(0.5));
end


