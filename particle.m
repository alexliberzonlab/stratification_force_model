%% Particle trajectory visualisation

%% 
clc 
clear all
close all

%% parameters
g    = 9.81;         % gravitational acceleration [m/s2]
zu   = 0.037;        % upper limit of the interface (entrance) [m]
zl   = 0.05;         % bottom limit of the interface (exit) [m]
rho1 = 976;          % density top layer [kg/m3]
rho2 = 1025;         % density bottom layer [kg/m3]
nu1  = 1.43e-6;      % viscosity top layer [m2/s]
nu2  = 1.012e-6;     % viscosity bottom layer [m2/s]
zl   = 0.05;         % bottom limit of the interface (exit) [m]
lam=0.35;
%% Load the experimental data
%Select the particle family (flag=1 =>P1,flag=2 =>P2, flag=3 =>P3,flag=4 =>P4, flag=5 => Fernando's particles) 
flag=4;

% Falling particles
if flag ==1;load('Data/TrajP1_complete.mat'); TRAJP =  eval(['TrajP',num2str(flag)]);
end
if flag ==2;load('Data/TrajP2_complete.mat'); TRAJP =  eval(['TrajP',num2str(flag)]);
end
if flag ==3;load('Data/TrajP3_complete.mat'); TRAJP =  eval(['TrajP',num2str(flag)]);
end


% Rising particles
if flag ==4;load('Data/TrajP4_complete.mat'); TRAJP =  eval(['TrajP',num2str(flag)]);
    g=-g;
    zu   = -0.003;       % upper limit of the interface (entrance) [m]
    zl   = 0.0075;       % bottom limit of the interface (exit) [m]
    rho2 = 976;          % density top layer [kg/m3]
    rho1 = 1025;         % density bottom layer [kg/m3]
    nu2  = 1.43e-6;      % viscosity top layer [m2/s]
    nu1  = 1.012e-6;     % viscosity bottom layer [m2/s]
end

% Fernando's particles
if flag ==5;load('Data/TrajPF_complete.mat');
    zu   = 0.125;         % upper limit of the interface (entrance) [m]
    zl   = 0.138;          % bottom limit of the interface (exit) [m]
    rho1 = 1015;          % density top layer [kg/m3]
    rho2 = 1045;          % density bottom layer [kg/m3]
    nu1  = 1.17e-06;      % viscosity top layer [m2/s]
    nu2  = 1.19e-06;      % viscosity bottom layer [m2/s]
    TRAJP =  eval(['TrajPF']); lam=0.25;
end

h=zl-zu;

Ind_vec=[1:length(TRAJP)];

  options = odeset('reltol', 1e-6, 'abstol', 1e-6);

figure
for i=1:length(Ind_vec)
     
    name_to_disp=['trajP',num2str(i)];
    disp(name_to_disp)
    
    rhop=TRAJP(i).rhop;
    d=TRAJP(i).d;
    z_exp  = TRAJP(i).z;
    vz_exp = TRAJP(i).vz;
    t_exp  = TRAJP(i).time;
    t0= t_exp(1);                     % initial time         [s]
    z0= z_exp(1);                     % initial position     [m]
    tend= t_exp(end);                 % final time           [s]
    

    
    % calculate the trajectory and plot the results
    [t, zp, V] = f_particle(z0, tend, rhop, d, g, zu, zl, rho1, rho2, nu1, nu2, lam, options);

    subplot(311)
    hold on
    plot(t_exp-t0,z_exp,'s','markersize',5)
    plot(t,zp,'r','linewidth',1.5)
    xlabel('t'), ylabel('z_p')
    
    subplot(312)
    hold on
    plot(t_exp-t0,vz_exp,'s','markersize',5)
    plot(t,V,'r','linewidth',1.5)
    xlabel('t'), ylabel('V')
    
    subplot(313)
    hold on
    plot(z_exp,vz_exp,'s','markersize',5)
    plot(zp,V,'r','linewidth',1.5)
    xlabel('z'), ylabel('V')
end
%% function f_particle

type f_particle

%% function settlingvelocity

type settlingvelocity