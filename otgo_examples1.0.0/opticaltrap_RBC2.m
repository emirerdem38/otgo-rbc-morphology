%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Medium refractive index

% RBC Particle
R = 3.91*1e-6; % [m]
t_min = 0.81*1e-6;
t_max = 2.52*1e-6;
d = 2.76*1e-6;
np = 1.50; % Particle refractive index

% Focusing
f = 5e-6; % Focal length [m]
NA = 1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 5e-6; % Beam waist [m]
Nphi = 2; % Azimuthal divisions
Nr = 2; % Radial divisions
power = 0.25e-3; % power [W]

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K]
eta = 0.001; % Water viscosity [Pa*s]
dt = 5e-3; % timestep [s]
N = 1e+3; % number of steps

% Initial position
x = 0;
y = 0;
z = 0;

%% Initialization

% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);
% rotate the beam: propagation direction along the x axis
r.v = r.v.yrotation(pi/2);
r.pol = r.pol.yrotation(pi/2);

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

        lncX = r.v.Vx(i,j)*A;
        lncY = r.v.Vy(i,j)*A;
        lncZ = r.v.Vz(i,j)*A;

        r.v.Vx(i,j) = r.v.Vx(i,j)*A;
        r.v.Vy(i,j) = r.v.Vy(i,j)*A;
        r.v.Vz(i,j) = r.v.Vz(i,j)*A;

    end
end

% Screen (RBC with radius f centered in the focal point)
screen = RBC(Point(0,0,0),Vector(0,0,0,0,0,1),R,t_min,t_max,d);

% Diffusion Constant of the Spherical Particle
D = kB*T / (6*pi*eta*R);

% %% Plot incoming beam
% 
% figure
% % Levels of the color map
% levels = 128;
% 
% hold on
% 
% % Intensity profile of the incoming beam in cartesian coordinates
% [xc,yc] = Transform.Pol2Car(bg.phi,bg.r);
% I = bg.intensity();
% 
% contourf([zeros(size(xc,1)+1,1),[xc;xc(1,:)]],...
%     [zeros(size(yc,1)+1,1),[yc;yc(1,:)]],...
%     [[I(:,1);I(1,1)],[I;I(1,:)]], ...
%     levels)
% colormap(ones(levels,3)-(0:1/(levels-1):1)'*[0 1 1])
% shading flat
% 
% % Sketch reference axis
% plot([max(max(xc)) min(min(xc))], [0 0], 'k')
% plot([0 0], [max(max(yc)) min(min(yc))], 'k')
% 
% % Sketch iris
% theta = (-1:0.01:1)*pi;
% xcirc = max(max(xc))*cos(theta);
% ycirc = max(max(xc))*sin(theta);
% plot(xcirc, ycirc, 'k', 'LineWidth', 1);
% 
% 
% axis equal tight off
% 
% drawnow()

%% Brownian motion simulation

figure
for n = 1:1:N
    
    % Display update message
    disp(['time = ' num2str(n*dt) ' s / ' num2str(N*dt) ' s'])
    
    % RBC particle
    particle = ParticleRBC(Point(1.193843550653809e-08,-4.394895341071864e-09,4.725171924600484e-08),Vector(0,0,0,0,0,1),R,t_min,t_max,d,nm,np);
    
    %% Optical force
    forces = particle.force(r);
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
    
    disp('----------------------------------------------')
    disp(force)
    disp(x)
    disp(y)
    disp(z)
    disp('----------------------------------------------')

    % Plot the particle
    particle.plot('scale', 1e+6, ...
        'facecolor', [0 0.75 0], ...
        'edgecolor', [0 0 0], ...
        'facealpha', .2, ...
        'edgealpha', 0.1 ...
        );

    axis equal
    view(3)
    plot(Point(particle.rbc.c.X,particle.rbc.c.Y,particle.rbc.c.Z)) % Plot the center
     
end