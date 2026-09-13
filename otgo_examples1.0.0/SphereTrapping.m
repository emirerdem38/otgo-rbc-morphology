% OPTICALTRAP Optical trap and forward scattering
%
% Simulates the Brownian motion of a numerically constructed spehrical particle 
% in an optical trap under the action of both optical forces and thermal forces.

%% Parameters

% Medium
nm = 1.33; % Medium refractive index

% Spherical Particle
p = 32;
shp = 'sphere';
R = 1.0e-06; % [m]
vol = 'vol';
np = 1.38; % Particle refractive index
or1 = Vector(0,0,0,1,0,0);
or2 = Vector(0,0,0,0,1,0);
or3 = Vector(0,0,0,0,0,1);

% Focusing
f = 5e-6; % Focal length [m]
NA = 1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 5e-6; % Beam waist [m]
Nphi = 20; % Azimuthal divisions
Nr = 40; % Radial divisions
power = 5e-3; % power [W]

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K]
eta = 0.001; % Water viscosity [Pa*s]
dt = 5e-3; % timestep [s]
N = 1e+3; % number of steps

% Initial position
x = -0.2*R;
y = 0.0;
z = -0.2*R;

%% Initialization

% Trapping beam
r = CreateBeam(nm, f, NA, Ex0, Ey0, w0, Nphi, Nr, power);

% rotate the beam: propagation direction along the x axis
% r.v = r.v.yrotation(pi/2);
% r.pol = r.pol.yrotation(pi/2);

% Diffusion Constant of the Spherical Particle
D = kB*T / (6*pi*eta*R);

% Initialize matrices
center_x = NaN(1,N);
center_y = NaN(1,N);
center_z = NaN(1,N);

% Brownian motion simulation
for n = 1:1:N
    
    % Display update message
    disp(['time = ' num2str(n*dt) ' s / ' num2str(N*dt) ' s'])
    sprintf('X:%g\n',x)
    sprintf('Y:%g\n',y)
    sprintf('Z:%g\n',z)

    center_x(1,n) = x;
    center_y(1,n) = y;
    center_z(1,n) = z;
    
    % Spherical particle
    bead = ParticleCELL(p,shp,Point(x,y,z),or1,or2,or3,vol,-6,nm,np);

    r_vec = scattering(bead,r);
    
    % Optical force
    forces = bead.force(r_vec,r);
    force = Vector(x,y,z, ...
        sum(forces.Vx(isfinite(forces.Vx))), ...
        sum(forces.Vy(isfinite(forces.Vy))), ...
        sum(forces.Vz(isfinite(forces.Vz))) ...
        );
    
    % Particle position update - optical force contribution
    x = x + dt*D*force.Vx/(kB*T);
    y = y + dt*D*force.Vy/(kB*T);
    z = z + dt*D*force.Vz/(kB*T);
    
    
    % Particle position update - thermal contribution
    x = x + sqrt(2*dt*D)*randn();
    y = y + sqrt(2*dt*D)*randn();
    z = z + sqrt(2*dt*D)*randn();

end

%% Visualization of the simulation
figure;
ax = axes;
ax = plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'DisplayName', 'Origin','MarkerFaceColor','r');

hold on

for i = 1:N
    cla(ax);
    % Update the cell center dynamically
    bead.cell.c = Point(center_x(1, i), center_y(1, i), center_z(1, i));
    
    % Plot the updated bead
    plot(bead);
    plot(r.toline);
    
    % Retain the axes and viewpoint
    grid on
    view(3);
    axis equal;
    axis([-2e-6 2e-6 -2e-6 2e-6 -2e-6 2e-6]); % Keep limits fixed
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    
    drawnow()
    pause(1)
end