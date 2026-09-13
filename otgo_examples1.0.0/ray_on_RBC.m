% RAY_ON_RBC Scattering of a ray on a RBC particle and relative optical forces
%
% Calculation of the scattering of a ray on a rbc particle and the
% relative optical forces as a function of its incidence angle. 
% The red arrows plot the transmitted and reflected rays.
% The black arrow represents the total optical force.
% The magenta arrow represents the scattering force.
% The blue arrow represents the gradient force.
%
% See also Vector, Ray, ParticleRBC.
%
% This example is part of the OTGO - Optical Tweezers in Geometrical Optics
% software package, which complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Agnese Callegari
%   Date: 2014/01/01
%   Version: 1.0.0

%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Medium refractive index

% Particle
R = 3.91; % Particle radius [m]
t_min = 0.81;
t_max = 2.52;
d = 2.76;
np = 1.50; % Particle refractive index
c = Point(0,0,0); % Particle center [m]
or = Vector(0,0,0,0,0,1); % Orientation of RBC

%% Initialization

% Particle
particle = ParticleRBC(c,or,R,t_min,t_max,d,nm,np);

% Figure
figure
box on
axis equal
grid on
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
title('Scattering and forces on a RBC')
hold on

%% Simulation
for theta = (0:1:89.9)/180*pi % Incidence angles [rad]
    disp("-------------------------------------")
    disp(theta)
    disp("-------------------------------------")
    % Ray
    v = Vector(-2*R,R*sin(theta),0,1,0,0); % Direction
    P = 1; % Power [W]
    pol = Vector(0,0,0,0,1,1);
    pol = v*pol;
    pol = pol.versor(); % Polarization
    r = Ray(v,P,pol);
    
    % Scattered rays
    s = particle.scattering(r,1e-12,6);
    
    % Scattering coefficients
    f = particle.force(s,r);
    
    % Clear plot
    cla;
    
    % Plot bead
    particle.plot('scale', 1e+6, ...
        'range', 20, ...
        'FaceColor', [0 0.75 0], ...
        'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.2, ...
        'EdgeAlpha', 0.2 ...
        );
    
    % Plot incident ray
    r.v.plot('scale', [1e+6 2-cos(theta)], ...
        'color', [1 0 0], ...
        'LineWidth', 2 ...
        );
    
    % Plot perpendicular line at the incidence point
    ip = Point(s(1).r.v.X,s(1).r.v.Y,s(1).r.v.Z);
    lin = particle.rbc.perpline(ip);
    lin.plot('range', [0 2e+6], ...
        'color', [0 0 0], ...
        'LineStyle', '-.', ...
        'LineWidth', 1 ...
        );
    
    % Plot first reflected ray
    s(1).r.v.plot('scale', [1e+6 1], ...
        'color', [1 0 0], ...
        'LineWidth', 2 ...
        );
    
    % Plot first transmitted ray
    s(1).t.v.plot('scale', [1e+6 2*sqrt(1-(nm/np*sin(theta)).^2)], ...
        'color', [1 0 0], ...
        'LineWidth', 2 ...
        );
    
    % Plot subsequent transmitted and reflected rays
    for j = 2:1:6

        % j-th transmitted ray
        s(j).t.v.plot('scale', [1e+6 1], ...
            'color', [1 0 0], ...
            'LineWidth', 2 ...
            );
        
        % j-th reflected ray
        s(j).r.v.plot('scale', [1e+6 2*sqrt(1-(nm/np*sin(theta)).^2)],...
            'color', [1 0 0], ...
            'LineWidth', 2 ...
            );
        
    end
    
    % Plot total force
    f.plot('scale', [1e+6 0.75e+9], ...
        'color', [0 0 0], ...
        'LineWidth', 2 ...
        );
    
    % Plot scattering force
    w = Vector(0,0,0,f.Vx,0,0);
    w.plot('scale', [1e+6 0.75e+9], ...
        'color', [1 0 1], ...
        'LineWidth', 2 ...
        );
    
    % Plot gradient force
    w = Vector(0,0,0,0,f.Vy,0);
    w.plot('scale', [1e+6 0.75e+9], ...
        'color', [0 0 1], ...
        'LineWidth', 2 ...
        );
       
    xlabel(['\theta = ' num2str(theta/pi*180) ' deg']);   

    drawnow(); %draw the output on the figure before passing to the next value of theta in the for-cycle
         
end