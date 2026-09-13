% EFFICIENCIES_RAY Trapping efficiencies for a single ray on a sphere
%
% Calculation of the trapping efficiencies corresponding to the scattering
% of a ray on a spherical particle as a function of the incidence angle.
%
% See also Vector, Ray, ParticleSpherical.
%
% This example is part of the OTGO - Optical Tweezers in Geometrical Optics
% software package, which complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

% Author: Agnese Callegari
% Date: 2014/01/01
% Version: 1.0.0

%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Medium refractive index

% Particle
R = 1e-6; % Particle radius [m]
np = 1.5; % Particle refractive index

%% Simulation

% Particle
c = Point(0,0,0); % Particle center [m]
bead = ParticleSpherical(c,R,nm,np);

% Rays
theta = [0:1:89.9]/180*pi; % Incidence angles [rad]
v = Vector(-2*R*ones(size(theta)),R*sin(theta),zeros(size(theta)),...
    ones(size(theta)),zeros(size(theta)),zeros(size(theta))); % Direction
P = ones(size(theta)); % Power [W]
pol = Vector(zeros(size(theta)),zeros(size(theta)),zeros(size(theta)),...
    zeros(size(theta)),ones(size(theta)),ones(size(theta)));
pol = v*pol;
pol = pol.versor(); % Polarization % Makes it a unit vector
r = Ray(v,P,pol);

% Scattering coefficients
f = bead.force(r,1e-12,6);
cm = PhysConst.c0/bead.nm; % Speed of light in medium [m/s]
Q = f./(r.P/cm); % Trapping efficiency

%% Figure 
figure

hold on

plot(theta/pi*180, sqrt(Q.Vx.^2+Q.Vy.^2), 'k', 'linewidth', 2.5)
plot(theta/pi*180,Q.Vy,'b','linewidth',2.5) %Gradient trapping efficiency
plot(theta/pi*180,Q.Vx,'m','linewidth',2.5) %Scattering trapping efficiency

title('Trapping efficiencies')
xlabel('\theta [deg]')
ylabel('Q, Q_g, Q_s')
legend('Q','Q_g','Q_s','Location','NorthWest')
box on
grid on
xlim([0 90])