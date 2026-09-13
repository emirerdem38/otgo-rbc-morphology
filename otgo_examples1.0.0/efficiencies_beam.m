% EFFICIENCIES_BEAM Scattering of a beam on a spherical particle
%
% Calculation of the trapping efficiencies corresponding to the scattering
% of a focused beam on a spherical particle as a function of the position
% of the particle with respect to the focal point.
%
% See also Ray, BeamGauss, ParticleSpherical.
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

% Spherical particle
R = 1e-6; % Particle radius [m]
np = 1.50; % Particle refractive index

% Focusing
f = 10e-6; % Focal length [m]
NA = 1.30; % numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 100e-6; % Beam waist [m]
Nphi = 40; % Azimuthal divisions
Nr = 40; % Radial divisions
power = 5e-3; % Power [W]

% Positions where to calculate the trapping efficiencies
z = [-1.6:0.02:1.6]*R; % Longitudinal (z) direction [m]
x = [0:0.02:1.6]*R; % Transversal (x) direction [m]

%% Initialization

% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to optical beam
r = Ray.beam2focused(bg,f);

%% Simulation - Longitudinal

% Preallocation variables

% Force components
fz_z = zeros(1,numel(z));
fz_x = zeros(1,numel(z));

% Trapping efficiencies
qz_x = zeros(1,numel(z));   
qz_z = zeros(1,numel(z));   

% Scattering and gradient forces
fz_s = zeros(3,numel(z)); % Scattering force
fz_g = zeros(3,numel(z)); % Gradient force

% Scattering and gradient trapping efficiencies
qz_s = zeros(3,numel(z)); % Scattering trapping efficiency
qz_g = zeros(3,numel(z)); % Gradient trapping efficiency

for i = 1:1:length(z)
    
    % Display update message
    disp(['Longitudinal force calculation ' int2str(i) '/' int2str(length(z))])
    
    % Spherical particle
    bead = ParticleSpherical(Point(0,0,z(i)),R,nm,np);
    
    % Calculate force
    forces = bead.force(r);
    force = Vector(0,0,z(i), ...
        sum(forces.Vx(isfinite(forces.Vx))), ...
        sum(forces.Vy(isfinite(forces.Vy))), ...
        sum(forces.Vz(isfinite(forces.Vz))) ...
        );
    
    fz_z(i) = force.Vz;
    fz_x(i) = force.Vx;
    
    qz_z(i) = fz_z(i)/power*(PhysConst.c0/nm);
    qz_x(i) = fz_x(i)/power*(PhysConst.c0/nm);
    
    % Scattering force
    fs = r.v .* (forces.*r.v);
    fz_s(1,i) = sum(fs.Vx(isfinite(fs.Vx)));
    fz_s(2,i) = sum(fs.Vy(isfinite(fs.Vy)));
    fz_s(3,i) = sum(fs.Vz(isfinite(fs.Vz)));
    
    % Scattering trapping efficiencies
    qz_s(1,i) = fz_s(1,i)/power*PhysConst.c0/nm;
    qz_s(2,i) = fz_s(2,i)/power*PhysConst.c0/nm;
    qz_s(3,i) = fz_s(3,i)/power*PhysConst.c0/nm;
    
    % Gradient force
    fg = forces - fs;
    fz_g(1,i) = sum(fg.Vx(isfinite(fg.Vx)));
    fz_g(2,i) = sum(fg.Vy(isfinite(fg.Vy)));
    fz_g(3,i) = sum(fg.Vz(isfinite(fg.Vz)));
    
    % Gradient trapping efficiencies
    qz_g(1,i) = fz_g(1,i)/power*PhysConst.c0/nm;
    qz_g(2,i) = fz_g(2,i)/power*PhysConst.c0/nm;
    qz_g(3,i) = fz_g(3,i)/power*PhysConst.c0/nm;
    
end

%% Figure - Longitudinal

figure

hold on

plot(z*1e+6,abs(qz_z),'k','LineWidth',1.5)
plot(z*1e+6,qz_g(3,:),'b','LineWidth',1.5)
plot(z*1e+6,qz_s(3,:),'m','LineWidth',1.5)

title('Trapping efficiencies - Longitudinal')
xlabel('z   [\mum]')
ylabel('Q, Q_g, Q_s [adimensional]')
legend('Q','Q_g','Q_s','Location','SouthWest')
box on
grid on
xlim([z(1) z(end)]*1e+6)

drawnow()

%% Simulation - Transversal

% Preallocation variables

% Force components
fx_z = zeros(1,numel(x));
fx_x = zeros(1,numel(x));

% Trapping efficiencies
qx_z = zeros(1,numel(x));
qx_x = zeros(1,numel(x));

% Scattering and gradient force
fx_s = zeros(3,numel(x)); % Scattering force
fx_g = zeros(3,numel(x)); % Gradient force

% Scattering and gradient trapping efficiencies
qx_s = zeros(3,numel(x)); % Scattering trapping efficiency
qx_g = zeros(3,numel(x)); % Gradient trapping efficiency

for i = 1:1:length(x)

    % Display update message
    disp(['Transversal force calculation ' int2str(i) '/' int2str(length(x))])

    % Spherical particle
    bead = ParticleSpherical(Point(x(i),0,0),R,nm,np);
    
    % Calculate force
    forces = bead.force(r);
    force = Vector(x(i),0,0, ...
        sum(forces.Vx(isfinite(forces.Vx))), ...
        sum(forces.Vy(isfinite(forces.Vy))), ...
        sum(forces.Vz(isfinite(forces.Vz))) ...
        );
    
    fx_z(i) = force.Vz;
    fx_x(i) = force.Vx;
    
    qx_x(i) = fx_x(i)/power*PhysConst.c0/nm;
    qx_z(i) = fx_z(i)/power*PhysConst.c0/nm;
    
    % Scattering force
    fs = r.v .* (forces.*r.v);
    fx_s(1,i) = sum(fs.Vx(isfinite(fs.Vx)));
    fx_s(2,i) = sum(fs.Vy(isfinite(fs.Vy)));
    fx_s(3,i) = sum(fs.Vz(isfinite(fs.Vz)));
    
    % Scattering trapping efficiencies    
    qx_s(1,i) = fx_s(1,i)/power*PhysConst.c0/nm;
    qx_s(2,i) = fx_s(2,i)/power*PhysConst.c0/nm;
    qx_s(3,i) = fx_s(3,i)/power*PhysConst.c0/nm;
    
    % Gradient force
    fg = forces - fs;
    fx_g(1,i) = sum(fg.Vx(isfinite(fg.Vx)));
    fx_g(2,i) = sum(fg.Vy(isfinite(fg.Vy)));
    fx_g(3,i) = sum(fg.Vz(isfinite(fg.Vz)));
    
    % Gradient trapping efficiencies    
    qx_g(1,i) = fx_g(1,i)/power*PhysConst.c0/nm;
    qx_g(2,i) = fx_g(2,i)/power*PhysConst.c0/nm;
    qx_g(3,i) = fx_g(3,i)/power*PhysConst.c0/nm;
    
end

%% Figure - Transversal

figure

hold on

plot([-x(end:-1:1) x]*1e+6,[abs(qx_x(end:-1:1)) abs(qx_x)],'k','LineWidth',1.5)
plot([-x(end:-1:1) x]*1e+6,[-qx_g(1,end:-1:1) qx_g(1,:)],'b','LineWidth',1.5)
plot([-x(end:-1:1) x]*1e+6,[qx_s(3,end:-1:1) qx_s(3,:)],'m','LineWidth',1.5)

title('Trapping efficiencies - Transversal')
xlabel('x [\mum]')
ylabel('Q, Q_g, Q_s [adimensional]')
legend('Q','Q_g','Q_s','Location','SouthWest')
box on
grid on
xlim([-x(end) x(end)]*1e+6)