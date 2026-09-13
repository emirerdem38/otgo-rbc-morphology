%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Water based medium

% Particle
mag = -6;
p = 32;
shp = 'fourBump';
vol = 'vol';

np = 1.38; % CELL Particle refractive index
m = 1e-11; % CELL mass

scat = 10;

% Focusing
f = 5*1e-6; % Focal length [m]
NA = 1.3; %1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 5e-6; % Beam waist [m]
Nphi = 20; % Radial divisions
Nr = 40; % Azimuthal divisions
power = 5e-3; % power [W]

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K]
eta = 0.001; % Water viscosity [Pa*s]
dt = 2e-3; % timestep [s]
N = 5000; % number of steps

%% Initialization

% Trapping beam
r = CreateBeam(nm, f, NA, Ex0, Ey0, w0, Nphi, Nr, power);

% rotate the beam: propagation direction along the x axis
% r.v = r.v.yrotation(pi/p2);
% r.pol = r.pol.yrotation(pi/2);s

load('diffusion_tensors.mat'); %  Load the diffusion tensors
Dname= "D_" + shp;
D = eval(Dname);

D = D_dumbbell*1e-7; % scale back to unit sizes

% Ensure a proper timescale
Dxx = D(1,1);
Dyy = D(2,2);
Dzz = D(3,3);

gamma_x = (kB*T)/Dxx; % Drag coefficient in x direction
gamma_y = (kB*T)/Dyy; % Drag coefficient in y direction
gamma_z = (kB*T)/Dzz; % Drag coefficient in z direction

gamma = min([gamma_x, gamma_y, gamma_z]);
   
tau_m = m / gamma; % Momentum relaxation time

k_x = 1.6*1e-6; % Trap stiffness in x direction (NOT CORRECT BUT KEPT FOR THE TIME BEING!)

tau_OT = gamma/k_x; % Time scale on which the restoring force acts

if dt > tau_OT || dt < tau_m
    fprintf('Time scale limits: %.4f s > Δt > %.4f s\n', tau_OT, tau_m);
    error('Please Adjust the time step!')
end

% Random displacement and orientation generation
seed = 42; % You can change this to any integer
rng(seed);

% Generate 3 rotation angles between 0 and pi/2
% rotation_angles = (pi/2) * rand(1, 3);
%rotation_angles = [0 90 0];
rotation_angles = [0.9670172108614332 0 0];

% Initial particle
shape_reference = ParticleCELL(p,shp,Point(0,0,0),Vector(0,0,0,1,0,0),Vector(0,0,0,0,1,0),Vector(0,0,0,0,0,1),vol,-6,nm,np);
or1 = shape_reference.cell.or1;
or2 = shape_reference.cell.or2;
or3 = shape_reference.cell.or3;

center = Point(0,0,0);

particle = ParticleCELL(p,shp,center,or1,or2,or3,vol,-6,nm,np);
particle.cell = particle.cell.xrotation(rotation_angles(1));
particle.cell = particle.cell.yrotation(rotation_angles(2));
particle.cell = particle.cell.zrotation(rotation_angles(3));

R = [particle.cell.or1.Vx particle.cell.or2.Vx particle.cell.or3.Vx;
    particle.cell.or1.Vy particle.cell.or2.Vy particle.cell.or3.Vy;
    particle.cell.or1.Vz particle.cell.or2.Vz particle.cell.or3.Vz];

% Transformation matrices
Mlp = R.';  % Lab to Particle (inverse of rotation matrix)
Mpl = R;    % Particle to Lab (rotation matrix)

L = chol(D, 'lower'); % 'lower' gives the lower triangular matrix

%forces = zeros(3,N);
%torques = zeros(3,N);
%centers = zeros(3,N);
%orientations_1 = zeros(3,N);
%orientations_2 = zeros(3,N);
%orientations_3 = zeros(3,N);

%% Dynamic Analysis

tic;

for n = 1:1:N
    
    r_vec = scattering(particle,r,0,scat-1);
    f = particle.force(r_vec,r);
    t = particle.torque(r_vec,r);

    force_lab = [sum(f.Vx(isfinite(f.Vx)));
             sum(f.Vy(isfinite(f.Vy)));
             sum(f.Vz(isfinite(f.Vz)))];
    torque_lab = [sum(t.Vx(isfinite(t.Vx)));
             sum(t.Vy(isfinite(t.Vy)));
             sum(t.Vz(isfinite(t.Vz)))];
    
    force_particle = Mlp * force_lab;
    torque_particle = Mlp * torque_lab;

    forces(:,n) = force_particle;
    torques(:,n) = torque_particle;

    X = particle.cell.c.X;
    Y = particle.cell.c.Y;
    Z = particle.cell.c.Z;

    centers(1,n) = X;
    centers(2,n) = Y; 
    centers(3,n) = Z;

    orientations_1(:,n) = Mlp(1,:)';
    orientations_2(:,n) = Mlp(2,:)';
    orientations_3(:,n) = Mlp(3,:)';

    z = randn(6, 1); % 6x1 vector of standard normal samples
    whiteNoise = L * z; % Transform to obtain the white noise vector

    d_particle = ((D*dt)/(kB*T)) * [force_particle; torque_particle] + sqrt(2*dt) * whiteNoise;

    new_center = [X; Y; Z] + Mlp' * d_particle(1:3);
    
    deltaAlpha = d_particle(4); deltaBeta = d_particle(5); deltaGamma = d_particle(6);

    particle.cell.c = Point(new_center(1),new_center(2),new_center(3)).zrotation(deltaAlpha);
    particle.cell.or1 = particle.cell.or1.zrotation(deltaAlpha);
    particle.cell.or2 = particle.cell.or2.zrotation(deltaAlpha);
    particle.cell.or3 = particle.cell.or3.zrotation(deltaAlpha);

    particle.cell.c = particle.cell.c.yrotation(deltaBeta);
    particle.cell.or1 = particle.cell.or1.yrotation(deltaBeta);
    particle.cell.or2 = particle.cell.or2.yrotation(deltaBeta);
    particle.cell.or3 = particle.cell.or3.yrotation(deltaBeta);

    particle.cell.c = particle.cell.c.xrotation(deltaGamma);
    particle.cell.or1 = particle.cell.or1.xrotation(deltaGamma);
    particle.cell.or2 = particle.cell.or2.xrotation(deltaGamma);
    particle.cell.or3 = particle.cell.or3.xrotation(deltaGamma);    

    Mlp = [particle.cell.or1.Vx  particle.cell.or1.Vy particle.cell.or1.Vz;
    particle.cell.or2.Vx particle.cell.or2.Vy particle.cell.or2.Vz;
    particle.cell.or3.Vx particle.cell.or3.Vy particle.cell.or3.Vz];
end

toc;

%%

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

% figure;
% hold on
% plot(particle.cell)
% for i = 1:size(orientations,2)
%     plot(Vector(0,0,0,orientations(1,i),orientations(2,i),orientations(3,i)),'Scale',1e-6)
% end

% Visualization of the simulation
figure;
ax = axes;
hold on

% Set axis limits (adjust based on your expected motion range)
xmin = min(centers(1,:)) - 1; xmax = max(centers(1,:)) + 1;
ymin = min(centers(2,:)) - 1; ymax = max(centers(2,:)) + 1;
zmin = min(centers(3,:)) - 1; zmax = max(centers(3,:)) + 1;

axis([xmin xmax ymin ymax zmin zmax]) % Fix axis limits
axis equal % Keep aspect ratio consistent
view(3) % Set a 3D view
p_vis = 32; % visualiztion precision

for i = 1:100:size(centers,2)
    cla(ax); % Clear only the axes, not the entire figure
    
    % Create cell object at the given position and orientation
    cell = Cell(p_vis, shp, Point(centers(1,i), centers(2,i), centers(3,i)), ...
                Vector(0,0,0, orientations_x(1,i), orientations_x(2,i), orientations_x(3,i)), ...
                Vector(0,0,0, orientations_y(1,i), orientations_y(2,i), orientations_y(3,i)), ...
                Vector(0,0,0, orientations_z(1,i), orientations_z(2,i), orientations_z(3,i)), vol, mag);
    
    plot(cell) % Plot the cell
    plot(r.toline)
    plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 

    % Hold the view stable
    drawnow()
    pause(0.1)
end

R1 = [1 0 0 0; 0 1/sqrt(2) -1/sqrt(2) 0; 0 1/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];
S = [1 0 0 0; 0 1/3 0 0; 0 0 1 0; 0 0 0 1];
T = [1 0 0 3; 0 1 0 -sqrt(2); 0 0 1 0; 0 0 0 1];
R2 = [1 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0; 0 -1/sqrt(2) 1/sqrt(2) 0; 0 0 0 1];

M = R2*T*S*R1;