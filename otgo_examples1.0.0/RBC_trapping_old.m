%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Water based medium

% RBC Particle
rad = 3.91*1e-6; % [m]
t_min = 0.81*1e-6;
t_max = 2.52*1e-6;
d = 2.76*1e-6;
np = 1.38; % RBC Particle refractive index
m = 1e-11; % RBC mass

% Focusing
f = 5*1e-6; % Focal length [m] % previously 5 % 200
NA = 1.3; %1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 5e-6; % Beam waist [m]
Nphi = 20; % Azimuthal divisions
Nr = 40; % Radial divisions
power = 1e-3; % power [W]

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K]
eta = 0.001; % Water viscosity [Pa*s]
dt = 1e-3; % timestep [s]
N = 1e+2; % number of steps

%% Initialization

% Trapping beam
r = CreateBeam(nm, f, NA, Ex0, Ey0, w0, Nphi, Nr, power);

% rotate the beam: propagation direction along the x axis
% r.v = r.v.yrotation(pi/p2);
% r.pol = r.pol.yrotation(pi/2);s

% Since the RBC equation and parameters are working in a larger scale, rays should be adjusted 
for i = 1:size(r,1)
    for j = 1:size(r,2)

        ordX = floor(log10(abs(r.v.X(i,j)))); % Order of Magnitude
        ordY = floor(log10(abs(r.v.Y(i,j))));
        ordZ = floor(log10(abs(r.v.Z(i,j))));
        if ordX < 0 && ordY < 0 && ordX < 0
            A = 10^(max([ordX,ordY,ordZ]));
        else
            A = 1;
        end

        % lncX = r.v.Vx(i,j)*A;
        % lncY = r.v.Vy(i,j)*A;
        % lncZ = r.v.Vz(i,j)*A;

        r.v.Vx(i,j) = r.v.Vx(i,j)*A;
        r.v.Vy(i,j) = r.v.Vy(i,j)*A;
        r.v.Vz(i,j) = r.v.Vz(i,j)*A;

    end
end

% Define the Difision Tensor (Results From Article)
Dtt = [ 7.43*1e-14 -4.83*1e-20 6.23*1e-21;
        5.93*1e-21  7.43*1e-14  5.99*1e-21;
       -8.74*1e-20 -5.42*1e-19 6.28*1e-14];

Drt = [-6.18*1e-15 -2.52*1e-15 -1.75*1e-15;
       -2.52*1e-15  8.85*1e-16 -2.72*1e-15;
       -1.75*1e-15 -2.72*1e-15 -2.20*1e-16];

Dtr = transpose(Drt);

Drr = [4.04*1e-3  3.63*1e-11  1.06*1e-10;
       1.04*1e-9  4.04*1e-3  -7.85*1e-10;
       1.02*1e-10 3.17*1e-10  3.36*1e-3];

D = [Dtt Dtr;
     Drt Drr];

% Ensure D is symmetric
D = (D + D.') / 2;

% D = ones(6,6);
% D = D*2.146098814590747e-13;

Dxx = D(1,1);
Dyy = D(2,2);
Dzz = D(3,3);

gamma_x = (kB*T)/Dxx; % Drag coefficient in x direction
gamma_y = (kB*T)/Dyy; % Drag coefficient in y direction
gamma_z = (kB*T)/Dzz; % Drag coefficient in z direction

gamma = min([gamma_x, gamma_y, gamma_z]);
   
tau_m = m / gamma; % Momentum relaxation time

k_x = 1.6*1e-6; % Trap stiffness in x direction

tau_OT = gamma/k_x;% Time scale on which the restoring force acts

if dt > tau_OT || dt < tau_m
    fprintf('Time scale limits: %.4f s > Δt > %.4f s\n', tau_OT, tau_m);
    error('Please Adjust the time step!')
end

Mx = [1 0 0;
      0 0 1;
      0 -1 0];

My = [0 0 1;
      0 1 0;
      -1 0 0];

% Initial RBC particle
particle = ParticleRBC(Point(0,0,0),Vector(0,0,0,0,1,0),rad,t_min,t_max,d,nm,np);

angles = orientation(particle.rbc,false);
x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 

Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];

R = Rz * Ry * Rx;

X = R*[1;0;0];
Y = R*[0;1;0];
Z = [particle.rbc.or.Vx;particle.rbc.or.Vy;particle.rbc.or.Vz];

Mlp = [X Y Z]; % Lab to Particle
Mpl = transpose(Mlp); % Particle to Lab

forces = [];
torques = [];
centers = [];
orientations = [];
particles = [];

%% Dynamic Analysis

figure;
ax = axes;

for n = 1:1:N
    
    % Display update message
    disp(['time = ' num2str(n*dt) ' s / ' num2str(N*dt) ' s'])

    % Display particle information
    disp(particle.rbc.c)
    disp(particle.rbc.or)

    cla(ax);
    plot3(ax, particle.rbc.c.X, particle.rbc.c.Y, particle.rbc.c.Z, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    hold on;
    xlim(ax, [-5 5]);
    ylim(ax, [-5 5]);
    zlim(ax, [-5 5]);
    xlabel(ax,'x')
    ylabel(ax,'y')
    zlabel(ax,'z')
    
    % Plot the particle
    particle.plot('scale', 1e+6, ...
        'facecolor', [0 0.75 0], ...
        'edgecolor', [0 0 0], ...
        'facealpha', .2, ...
        'edgealpha', 0.1 ...
        );
    drawnow;
      
    r_vec = scattering(particle,r,1e-12,4);
    f = particle.force(r_vec,r);
    t = particle.torque(r_vec,r);

    force_lab = [sum(f.Vx(isfinite(f.Vx)));
                 sum(f.Vy(isfinite(f.Vy)));
                 sum(f.Vz(isfinite(f.Vz)));
                ];

    torque_lab = [sum(t.Vx(isfinite(t.Vx)));
                  sum(t.Vy(isfinite(t.Vy)));
                  sum(t.Vz(isfinite(t.Vz)));
                 ];

    force_particle = Mlp * force_lab;
    torque_particle = Mlp * torque_lab;

    plot(Vector(particle.rbc.c.X,particle.rbc.c.Y,particle.rbc.c.Z,force_particle(1),force_particle(2),force_particle(3)),'Scale',1e+7)

    forces = [forces, [force_particle(1); force_particle(2); force_particle(3)]];
    torques = [torques, [torque_particle(1); torque_particle(2); torque_particle(3)]];
    centers = [centers, [particle.rbc.c.X; particle.rbc.c.Y; particle.rbc.c.Z]];
    orientations = [orientations, [particle.rbc.or.Vx; particle.rbc.or.Vy; particle.rbc.or.Vz]];
    particles = [particles particle.rbc];

    % White noise associated with translation and rotation with respect to x,y,z axes
    whiteNoise = mvnrnd(zeros(1,6), D)';

    d_particle = ((D*dt)/(kB*T)) * [force_particle; torque_particle] + sqrt(2*dt) * whiteNoise;
    
    % Update center coordinates
    new_center = [particle.rbc.c.X; particle.rbc.c.Y; particle.rbc.c.Z] + Mpl * d_particle(1:3);
    particle.rbc.c = Point(new_center(1),new_center(2),new_center(3));
    
    deltaAlpha = d_particle(4); deltaBeta = d_particle(5); deltaGamma = d_particle(6);

    % Calculate the rotation axis vector
    d_axis = [deltaAlpha; deltaBeta; deltaGamma];

    % Normalize the rotation axis vector
    d_axis = d_axis / norm(d_axis);

    % Calculate the rotation angle
    angle = norm(axis);

    % Skew-symmetric matrix of the rotation axis vector
    skewMatrix = [0, -d_axis(3), d_axis(2);
                  d_axis(3), 0, -d_axis(1);
                  -d_axis(2), d_axis(1), 0];

    % Calculate the rotation matrix using the Rodriguez rotation formula
    rotationMatrix = eye(3) + sin(angle) * skewMatrix + (1 - cos(angle)) * skewMatrix^2;
    
    % Update RBC orientation, Lab to Particle and Particle to Lab matrices
    Mlp = Mlp * rotationMatrix;
    particle.rbc.or = normalize(Vector(0,0,0,Mlp(1,3),Mlp(2,3),Mlp(3,3)));
    Mpl = transpose(Mlp); % Particle to Lab
    
end

figure;

subplot(2,3,1)
title('Forces in x axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(forces(1,:))

subplot(2,3,2)
title('Forces in y axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(forces(2,:))

subplot(2,3,3)
title('Forces in z axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(forces(3,:))

subplot(2,3,4)
title('Torques in x axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(torques(1,:))

subplot(2,3,5)
title('Torques in y axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(torques(2,:))

subplot(2,3,6)
title('Torques in z axis')
xlabel('x')
ylabel('y')
zlabel('z')
plot(torques(3,:))

figure;
plot3(centers(1,:), centers(2,:), centers(3,:), 'rx', 'MarkerSize', 4);
hold on
plot(particle)

figure;
hold on
plot(particle.rbc)
for i = 1:size(orientations,2)
    plot(Vector(0,0,0,orientations(1,i),orientations(2,i),orientations(3,i)),'Scale',1e-6)
end

% Visualization of the simulation
figure;
ax = axes;
hold on
for i=1:size(centers,2)
    cla(ax);
    rbc = RBC(Point(centers(1,i),centers(2,i),centers(3,i)),Vector(0,0,0,orientations(1,i),orientations(2,i),orientations(3,i)),rad,t_min,t_max,d);
    plot(rbc)
    drawnow()
    pause(1)
end
