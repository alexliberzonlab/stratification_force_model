% f_particle calculates the trajectory of an inertial particle traversing
% a stratified fluid layer. The z-axis is defined downwards, and
% layer 1 (the top layer) is for z<zu and layer 2 (the bottom layer) is for
% z>zl. 
%
% SYNTAX [t, zp, V] = f_particle(z0, tend, rhop, d, g, 
%                                zu, zl, rho1, rho2, nu1, nu2, fitfun, 
%                                options)
%
% INPUT ARGUMENTS
% ---------------
%        z0 : initial particle position        [m]
%      tend : simulation time                  [s]
%      rhop : particle density                 [kg/m3]
%         d : particle diameter                [m]
%         g : gravitational acceleration       [m/s2]
%        zu : start of density interface       [m]
%        zl : end of density interface         [m]
%      rho1 : density of top layer             [kg/m3]
%      rho2 : density of bottom layer          [kg/m3]
%       nu1 : kinematic viscosity of top layer [m2/s]
%       nu2 : kinematic viscosity of bottom layer [m2/s]
%    fitfun : a function of with arguments (z, zu, zl, f1, f2) to create
%             continuous density and viscosity profiles
%   options : ODE option function
%
% OUTPUT ARGUMENTS
% ----------------
%         t : time interval     [s]
%        zp : particle position [m]
%         V : particle velocity [m/s]

function [t, zp, V] = f_particle(z0, tend, rhop, d, g, ...
                                 zu, zl, rho1, rho2, nu1, nu2, lam,options)

Vp = (pi*d^3)/6;                             % volume of the sphere [m3]
h  = zl-zu;                                  % interface thickness  [m]
N  = (2*g*(rho2-rho1)/h/(rho1+rho2))^0.5;    % buoyancy frequency   [1/s]

% Construct density and viscosity functions
rho = @(z) rho2 - 0.5*(rho2-rho1)*(1-tanh((z-0.5*(zl + zu))/(lam*h)));
nu  = @(z)  nu2 + 0.5*(nu1-nu2)*(1-tanh((z-0.5*(zl + zu))/(lam*h)));

% Auxiliary functions
Cd         = @(Re) 0.25 + (24./Re) + (6./(1+Re.^(0.5)));  % drag law according to White 1974
trec_d2nu2 = @(Re) 13./Re;                                % recovery time 
Vc0_f      = @(Fr,Re,Vp)  0.13*Fr.^0.75 * Vp + 0.*Re;     % caudal volume

% calculate settling velocity
V1 = settlingvelocity(rhop,rho1,g,d,nu1);
V2 = settlingvelocity(rhop,rho2,g,d,nu2);

% Determine Fr1, Vc0 and trec
Fr1  = abs(V1) / (N * d);
Re1  = abs(V1) * d / (nu1);
trec = trec_d2nu2(V2*d/nu2)*(d^2/nu2);
Vc0  = max(Vc0_f(Fr1,Re1,Vp),0);

% Solve equation of motion
y0 = [z0; V1];  
odefun = @(t, y) particle_ode(t, y,rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0,trec);
[t, y] = ode15s(odefun, [0 tend], y0, options);
zp  = y(:, 1);  V  = y(:, 2);

function dydt = particle_ode(t, y, rho1, rho, nu , Cd, g, rhop, d, ...
                             zu, zl, Vc0,trec)
                             
zp   = y(1);                % particle position      [m]
V    = y(2);                % particle velocity      [m/s] 
Ap=(pi*d^2)/4;              % surface of the sphere  [m2]
Vp=(pi*d^3)/6;              % volume  of the sphere  [m3]
Cam = 0.5;                  % added mass coefficient [-]
global tzl                  % time at which the particle reaches zl

% determine Vc
if (zp <= zl)
    Vc = Vc0;
    tzl = t;
else
    Vc = Vc0 * exp(-(t-tzl)/trec);
end

% stratification force
FS = (rho(zp) - rho1)*Vc*g;

% force balance on the particle 
dzpdt = V;
dVdt  = ( (rhop-rho(zp)) * Vp * g ...
        - 0.5 *Cd(V*d/nu(zp))* rho(zp) * Ap * abs(V) * V ...
        - FS ...
        ) / (rhop*Vp + Cam*rho(zp)*Vp);

dydt = [dzpdt; dVdt;];
