% f_particle calculates the trajectory of an inertial particle traversing
% through a stratified fluid layer. The z-axis is defined downwards, and
% layer 1 (the top layer) is for z<zu and layer 2 (the bottom layer) is for
% z>zl. 
%
% SYNTAX [t, zp, V] = f_particle(z0, tend, rhop, d, g, zu, zl, rho1, rho2, nu1, nu2)
%
% INPUT ARGUMENTS
% ---------------
%        z0 : initial particle position [m]
%      tend : simulation time [s]
%      rhop : particle density [kg/m3]
%         d : particle diameter [m]
%         g : gravitational acceleration [m/s2]
%        zu : start of density interface [m]
%        zl : end of density interface [m]
%      rho1 : density of top layer [kg/m3]
%      rho2 : density of bottom layer [kg/m3]
%       nu1 : kinematic viscosity of top layer [m2/s]
%       nu2 : kinematic viscosity of bottom layer [m2/s]
%
% OUTPUT ARGUMENTS
% ----------------
%         t : time interval [s]
%        zp : particle position [m]
%         V : particle velocity [m/s]

function [t, zp, V] = f_particle(z0, tend, rhop, d, g, zu, zl, rho1, rho2, nu1, nu2,lam)

Vp = (pi*d^3)/6;                             % volume of the sphere [m3]
h  = zl-zu;                                  % interface thickness  [m]
N  = (2*g*(rho2-rho1)/h/(rho1+rho2))^0.5;    % buoyancy frequency [1/s]

% Construct density and viscosity functions
% lam = 0.35;
rho = @(z) rho2 - 0.5*(rho2-rho1)*(1-tanh((z-0.5*(zl + zu))/(lam*h)));
nu  = @(z)  nu2 + 0.5*(nu1-nu2)*(1-tanh((z-0.5*(zl + zu))/(lam*h)));

% Visualise the density and viscosity functions

% z=linspace(0,0.12,500);  % z linspace vector z-axis
% figure
% subplot(1,2,1)
% plot(rho(z),z, ...
%      [rho1, rho2], zu * [1,1], '--k', [rho1, rho2], zl * [1,1], '--k'); 
% xlabel('\rho'); ylabel('z')
% subplot(1,2,2)
% plot(nu(z), z, ...
%     [nu1, nu2], zu * [1,1], '--k', [nu1, nu2], zl * [1,1], '--k'); 
% xlabel('\nu'); ylabel('z')

% Auxiliary functions
Cd         = @(Re) 0.25 + (24./Re) + (6./(1+Re.^(0.5)));  % drag law
Vc0_Vp     = @(Fr) 1-(1.85)./Fr;                               % caudal volume
trec_d2nu2 = @(Re) (41.439).*Re.^(-1.1884);                                % recovery time

% calculate settling velocity
V1 = settlingvelocity(rhop,rho1,g,d,nu1);
V2 = settlingvelocity(rhop,rho2,g,d,nu2);

% Determine Fr1, Vc0 and trec
Fr1 = abs(V1) / (N * d);
Vc0 = Vc0_Vp(Fr1) * Vp;
trec = trec_d2nu2(V2*d/nu2)*(d^2/nu2);

% Solve equation of motion
y0 = [z0; V1];  
options = odeset('reltol', 1e-9, 'abstol', 1e-6);
[t, y] = ode15s(@particle_ode, [0 tend], y0, options, ...
                  rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0,trec);
zp  = y(:, 1);  V  = y(:, 2);

function dydt = particle_ode(t,y,rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0,trec)
zp   = y(1);                % particle position
wp   = y(2);                % particle velocity

Ap=(pi*d^2)/4;              % surface of the sphere [m2]
Vp=(pi*d^3)/6;              % volume  of the sphere [m3]
Cam = 0.5;                  % added mass coefficient [-]
global tzl                  % time at which the particle reaches zl

% determine Vc
if (zp <= zl)
    Vc = Vc0;
    tzl = t;
else
    Vc = Vc0 * exp(-3*(t-tzl)/trec);
end

% stratification force
FS = (rho(zp) - rho1)*Vc*g;

% force balance on the particle 
dzpdt = wp;
dwpdt = ( (rhop-rho(zp)) * Vp * g ...
        - 0.5 *Cd(wp*d/nu(zp))* rho(zp) * Ap * abs(wp) * wp ...
        - FS ...
        ) / (rhop*Vp + Cam*rho(zp)*Vp);

dydt = [dzpdt; dwpdt;];