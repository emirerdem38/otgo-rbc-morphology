% ELLIPSOID Brownian motion of an ellipsoidal particle in a optical trap
%
% Brownian motion of an optically trapped ellipsoidal particle. It can be
% seen how the particle gets alligned along the longitudinal direction.
%
% See also Point, Vector, SLine, Ray, BeamGauss, ParticleEllipsoidal.
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
nm = 1.33;  % Medium refractive index

% Ellipsoid
a = 2.00000e-06; % semiaxis [m]
b = 2.00000e-06; % semiaxis [m]
c = 3.33333e-06; % semiaxis [m]
% Initial position 
x0 = 0.9e-6; % [m]
y0 = 0.0e-6; % [m]
z0 = 1.4e-6; % [m]
% Initial orientation of the semiaxes
roty0 = 80*pi/180; % [rad]
rotx0 = 0; % [rad]

np = 1.5; % Particle refractive index

% Diffusion Matrix (from Hydro++)
% already symmetrized and in SI units
D = [ ...
    [ 1.74400e-12   1.17855e-18   1.58755e-18   4.48200e-14  -6.19300e-14   1.45600e-13 ]; ...
    [ 1.17855e-18   1.74400e-12   1.94700e-19  -6.19300e-14  -8.48100e-14  -2.62100e-12 ]; ...
    [ 1.58755e-18   1.94700e-19   1.93400e-12   1.45600e-13  -2.62100e-12  -4.07800e-14 ]; ...
    [ 4.48200e-14  -6.19300e-14   1.45600e-13   7.28700e-02   1.05490e-09  -5.85000e-10 ]; ...
    [-6.19300e-14  -8.48100e-14  -2.62100e-12   1.05490e-09   7.28600e-02   8.34700e-09 ]; ...
    [ 1.45600e-13  -2.62100e-12  -4.07800e-14  -5.85000e-10   8.34700e-09   1.08100e-01 ] ...
    ];

% Focusing
f = 100e-6;	% Focal length [m]
NA = 1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 50e-6; % Beam waist [m]
Nphi = 20; % Azimuthal divisions
Nr = 20; % Radial divisions
power = 5e-3; % power [W]

% Brownian motion
T = 293; % Temperature [K]
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
eta = 0.001; % viscosity [Pa*s] = [10 P] (water)
dt = 1e-3; % timestep [s]
N = 4e+2; % number of steps


%% Initialization

% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);

% Incident rays - starting coordinates
z0ray = 5e-6;
stpoin = Point((r.v.Vx ./ r.v.Vz)*(-z0ray), (r.v.Vy ./ r.v.Vz)*(-z0ray), ones(r.v.size())*(-z0ray));

% Particle
sa = Vector(0,0,0,1,0,0);
sa = sa.yrotation(roty0);
sa = sa.xrotation(rotx0);

sb = Vector(0,0,0,0,1,0);
sb = sb.yrotation(roty0);
sb = sb.xrotation(rotx0);

sc = Vector(0,0,0,0,0,1);
sc = sc.yrotation(roty0);
sc = sc.xrotation(rotx0);

bead = ParticleEllipsoidal(Point(x0,y0,z0),sa*a,sb*b,sc*c,nm,np);

%% Simulation - Inizialization

x = x0; % x-coordinate center of mass [m]
y = y0; % y-coordinate center of mass [m]
z = z0; % z-coordinate center of mass [m]

u1 = [sa.Vx; sa.Vy; sa.Vz]; % 1st direction cosine
u2 = [sb.Vx; sb.Vy; sb.Vz]; % 2nd direction cosine
u3 = [sc.Vx; sc.Vy; sc.Vz]; % 3rd direction cosine


% Changing frame reference: Lab->Particle   and   Particle ->Lab
% Transformation matrix from PARTICLE TO LAB reference frame
Mp2l = [u1 u2 u3];
% Transformation matrix from LAB TO PARTICLE reference frame
Ml2p = Mp2l';

% Rotation matrix - instantaneous rotation of the principal axes of the
% particle due to the Brownian noise
% initially R is set as unitary matrix
R = [1 0 0; 0 1 0; 0 0 1];

%% Figure - Inizialization

figure

box on 
grid off
axis equal
xlim([-1 1]*5);
ylim([-1 1]*5);
zlim([-1 1]*5);
view(0,0);

%% Simulation
for n = 1:1:N
    
    % Display update message
    disp(['time = ' num2str(n*dt) ' s / ' num2str(N*dt) ' s'])

    % Update the position of the bead
    bead.elli.c.X = x;
    bead.elli.c.Y = y;
    bead.elli.c.Z = z;
    
    % Update the orientation of the bead
    bead.elli.sa.Vx = a*u1(1);
    bead.elli.sa.Vy = a*u1(2);
    bead.elli.sa.Vz = a*u1(3);
    bead.elli.sb.Vx = b*u2(1);
    bead.elli.sb.Vy = b*u2(2);
    bead.elli.sb.Vz = b*u2(3);
    bead.elli.sc.Vx = c*u3(1);
    bead.elli.sc.Vy = c*u3(2);
    bead.elli.sc.Vz = c*u3(3);
    
    % Force due to the ray on the particle - LAB reference frame
    forces = bead.force(r,1e-12);
    force = Vector(x,y,z, ...
        sum(forces.Vx(isfinite(forces.Vx))), ...
        sum(forces.Vy(isfinite(forces.Vy))), ...
        sum(forces.Vz(isfinite(forces.Vz))) ...
        );
    
    fx = force.Vx;
    fy = force.Vy;
    fz = force.Vz;
    
    % Torque due to the ray on the particle - LAB reference frame
    torques = bead.torque(r,1e-12);
    torque = Vector(x,y,z, ...
        sum(torques.Vx(isfinite(torques.Vx))), ...
        sum(torques.Vy(isfinite(torques.Vy))), ...
        sum(torques.Vz(isfinite(torques.Vz))) ...
        );
    
    tx = torque.Vx;
    ty = torque.Vy;
    tz = torque.Vz;
    
    % Force due to the ray on the particle - PART reference frame
    f_lab = [force.Vx; force.Vy; force.Vz];
    f_part = Ml2p * f_lab;
    
    % Torque due to the ray on the particle - PART reference frame
    t_lab = [torque.Vx; torque.Vy; torque.Vz];
    t_part = Ml2p * t_lab;
    
    FTpart = [f_part(1) f_part(2) f_part(3) t_part(1) t_part(2) t_part(3)];
    
    % Now FTpart (generalized force and torque) is set in the proper way
    dQ_optical = dt*D*FTpart'/(kB * T);
    
    % Thermal noise - PART reference frame
    dQ_thermal = sqrt(2*dt)*mvnrnd(zeros(1,6), D)'; %stochastic (brownian) increments
    
    % Find the new increments for the generalized coordinates
    dq = dQ_thermal + dQ_optical;
    
    % Spatial
    dx = dq(1);       dy = dq(2);       dz = dq(3);
    % Angular
    dphi1 = dq(4);    dphi2 = dq(5);    dphi3 = dq(6);
    
    % Displacement of the center of diffusion in the LAB reference frame
    x = x + Mp2l(1,:)*[dx;dy;dz];
    y = y + Mp2l(2,:)*[dx;dy;dz];
    z = z + Mp2l(3,:)*[dx;dy;dz];
       
    % Find the new orientation of the axes on the particle
    R1 = [1 0 0; 0 cos(dphi1) -sin(dphi1); 0 sin(dphi1) cos(dphi1)];
    R2 = [cos(dphi2) 0 sin(dphi2); 0 1 0; -sin(dphi2) 0 cos(dphi2)];
    R3 = [cos(dphi3) -sin(dphi3) 0; sin(dphi3) cos(dphi3) 0; 0 0 1];
    R = R1*R2*R3;
    % the columns of R are the new director cosines in the PARTICLE reference frame
    
    % New orientation of the particle's axes in the LAB reference frame
    u1 = Mp2l*R(:,1);
    u2 = Mp2l*R(:,2);
    u3 = Mp2l*R(:,3);
    
    % Calculate the new rotation matrices for the next step
    % Transformation matrix from PARTICLE TO LAB reference frame
    Mp2l = [u1 u2 u3];
    % Transformation matrix from LAB TO PARTICLE reference frame
    Ml2p = Mp2l';
    
    
    % Plot the configuration at timestep n
    cla
    hold on
    
    bead.plot('scale', 1e+6, ...
        'FaceColor', [0 0.75 0], ...
        'EdgeColor', [0.0 0.0 0.0], ...
        'FaceAlpha', 0.2, ...
        'EdgeAlpha', 0.2 ...
        );
    
    scatt = bead.scattering(r,1e-12,2);
    
    % First scattering event
    intpoin1 = Point(scatt(1).t.v.X,scatt(1).t.v.Y,scatt(1).t.v.Z);
    
    % Second scattering event
    intpoin2 = Point(scatt(2).t.v.X,scatt(2).t.v.Y,scatt(2).t.v.Z);
    
    % Transmitted rays: final coordinates
    fipoin = Point( scatt(2).t.v.X + (z0ray-scatt(2).t.v.Z) .* (scatt(2).t.v.Vx ./ scatt(2).t.v.Vz),...
        scatt(2).t.v.Y + (z0ray-scatt(2).t.v.Z) .* (scatt(2).t.v.Vy ./ scatt(2).t.v.Vz),...
        z0ray*ones(scatt(2).t.v.size()));
    
    % Build the line respresenting the rays
    lin1 = SLine(stpoin,intpoin1);
    lin2 = SLine(intpoin1,intpoin2);
    lin3 = SLine(intpoin2,fipoin);
    
    % Plot
    lin1.plot('scale', 1e+6, 'color', 'r', 'LineWidth', 0.5);
    lin2.plot('scale', 1e+6, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    lin3.plot('scale', 1e+6, 'color', [1 0.5 0.3], 'LineWidth', 0.5);
    
    title(['t = ' num2str(n*dt) 's']);
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    zlabel('z [\mum]')

    drawnow();
    
end