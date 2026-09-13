% KRAMERS Kramers transition 
%
% Simualtions of thermally driven Kramers transitions between two optical traps.
%
% See also Point, Vector, SLine, Ray, BeamGauss, Spherical, ParticleSpherical.
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

% Spherical Particle
np = 1.50; % Particle refractive index
R = 0.8e-6; % Radius of the particle [m]

% Focusing
f = 100e-6;	% Focal length [m]
NA1 = 1.3; % numerical aperture first trapping beam
NA2 = 1.3; % numerical aperture second trapping beam
L1 = f*NA1/nm; % Iris aperture [m] first objective
L2 = f*NA2/nm; % Iris aperture [m] second objective

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 100e-6; % Beam waist [m]
Nphi = 20; % Azimuthal divisions
Nr = 20; % Radial divisions
power1 = 0.25e-3;  % [W] power of the first  trapping beam
power2 = 0.25e-3;  % [W] power of the second trapping beam

% Brownian motion
T = 293; % Temperature [K]
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
eta = 0.001; % Water viscosity [Pa*s]
dt = 5e-3; % timestep [s]
N = 5e+3; % number of steps

%% Initialization

% Particle
bead = ParticleSpherical(Point(0,0,0),R,nm,np);

% Trapping beams
bg1 = BeamGauss(Ex0,Ey0,w0,L1,Nphi,Nr);
bg1 = bg1.normalize(power1); % Set the power
bg2 = BeamGauss(Ex0,Ey0,w0,L2,Nphi,Nr);
bg2 = bg2.normalize(power2); % Set the power

% Calculates set of rays corresponding to each optical beam
r1 = Ray.beam2focused(bg1,f);
r2 = Ray.beam2focused(bg2,f);

% Distance between the two trapping beams
dist12 = 1.5 * R;

% focal points
o_x1 = -dist12/2;  % x coordinate of the focal point of the first  beam
o_x2 =  dist12/2;  % x coordinate of the focal point of the second beam

% translation along the x axis
r1.v.X = r1.v.X + o_x1;
r2.v.X = r2.v.X + o_x2;

% combining the two bunch of rays 
% for higher computational speed
r = Ray(Vector([r1.v.X; r2.v.X],[r1.v.Y; r2.v.Y],[r1.v.Z; r2.v.Z],...
    [r1.v.Vx; r2.v.Vx],[r1.v.Vy; r2.v.Vy],[r1.v.Vz; r2.v.Vz]),...
    [r1.P; r2.P],...
    Vector([r1.pol.X; r2.pol.X],[r1.pol.Y; r2.pol.Y],[r1.pol.Z; r2.pol.Z],...
    [r1.pol.Vx; r2.pol.Vx],[r1.pol.Vy; r2.pol.Vy],[r1.pol.Vz; r2.pol.Vz]));

% combining the real part ONLY of the two bunch of rays in a single bunch
% for higher computational speed while PLOTTING
rp = Ray(Vector([r1.v.X(1:Nphi,:); r2.v.X(1:Nphi,:)],...
    [r1.v.Y(1:Nphi,:); r2.v.Y(1:Nphi,:)],...
    [r1.v.Z(1:Nphi,:); r2.v.Z(1:Nphi,:)],...
    [r1.v.Vx(1:Nphi,:); r2.v.Vx(1:Nphi,:)],...
    [r1.v.Vy(1:Nphi,:); r2.v.Vy(1:Nphi,:)],...
    [r1.v.Vz(1:Nphi,:); r2.v.Vz(1:Nphi,:)]),...
    [r1.P(1:Nphi,:); r2.P(1:Nphi,:)],...
    Vector([r1.pol.X(1:Nphi,:); r2.pol.X(1:Nphi,:)],...
    [r1.pol.Y(1:Nphi,:); r2.pol.Y(1:Nphi,:)],...
    [r1.pol.Z(1:Nphi,:); r2.pol.Z(1:Nphi,:)],...
    [r1.pol.Vx(1:Nphi,:); r2.pol.Vx(1:Nphi,:)],...
    [r1.pol.Vy(1:Nphi,:); r2.pol.Vy(1:Nphi,:)],...
    [r1.pol.Vz(1:Nphi,:); r2.pol.Vz(1:Nphi,:)]));

% Diffusion Constant of the Spherical Particle
D = kB*T / (6*pi*eta*R);

% initial position of center of the bead with respect to the origin
x0 = o_x1; % the bead is in the focus of the first beam
y0 = 0.0;
z0 = 0.0;

%% Simulation - Preallocation and initialization

x = zeros(1,N); % x-coordinate center of the bead [m]
y = zeros(1,N); % y-coordinate center of the bead [m]
z = zeros(1,N); % z-coordinate center of the bead [m]

% Initial position
x(1)=x0;  y(1)=y0;  z(1)=z0;

%% Figure - Initialization

figure('Units', 'Pixels', 'Position', [100 100 1000 660])

focalpoint1 = Point(o_x1,0,0);
focalpoint2 = Point(o_x2,0,0);

z0ray = 2.5*R*1e+6;

% Incident rays: starting coordinates (beam 1)
stpoin1 = Point((r1.v.Vx ./ r1.v.Vz)*(-z0ray)+o_x1,...
    (r1.v.Vy ./ r1.v.Vz)*(-z0ray), ones(r1.v.size())*(-z0ray));

% Incident rays: starting coordinates (beam 2)
stpoin2 = Point((r2.v.Vx ./ r2.v.Vz)*(-z0ray)+o_x2,...
    (r2.v.Vy ./ r2.v.Vz)*(-z0ray), ones(r2.v.size())*(-z0ray));

% Incident rays: starting coordinates (beams 1 2)
stpoin = Point([stpoin1.X(1:Nphi,:); stpoin2.X(1:Nphi,:)],...
    [stpoin1.Y(1:Nphi,:); stpoin2.Y(1:Nphi,:)],...
    [stpoin1.Z(1:Nphi,:); stpoin2.Z(1:Nphi,:)]);

% Screen (to determine the points where scattered rays strike on)
screen = Spherical(Point(0,0,0),10*R);

% Define the plot axes for the 3D plot
h_3d = axes ('Units','Pixels','Position',[10 10 600 600],'Color','w');

box on
grid off
axis on
xlim([-1 1]*2.5*R*1e+6)
ylim([-1 1]*2.5*R*1e+6)
zlim([-1 1]*2.5*R*1e+6)
view(0,0)

% Define the plot axes for x(t)
h_x = axes ('Units','Pixels','Position',[660 430 330 180],'Color','w');

box on
grid off
axis on
ylim([-1 1]*1.5*R*1e+6)
ylabel('x(t) [\mum]')

% Define the plot axes for y(t)
h_y = axes ('Units','Pixels','Position',[660 225 330 180],'Color','w');

box on
grid off
axis on
ylim([-1 1]*1.5*R*1e+6)
ylabel('y(t) [\mum]')

% Define the plot axes for z(t)
h_z = axes ('Units','Pixels','Position',[660 20 330 180],'Color','w');

box on
grid off
axis on
ylim([-1 1]*1.5*R*1e+6)
ylabel('z(t) [\mum]')

%% Simulation
for n = 1:1:N
    
    % Update the position of the bead
    bead.sp.c.X = x(n);
    bead.sp.c.Y = y(n);
    bead.sp.c.Z = z(n);
    
    % Force due to the first beam on the particle
    forces = bead.force(r);
    force = Vector(x(n),y(n),z(n), ...
        sum(forces.Vx(isfinite(forces.Vx))), ...
        sum(forces.Vy(isfinite(forces.Vy))), ...
        sum(forces.Vz(isfinite(forces.Vz))) ...
        );
            
    % Find the new increments for the coordinates
    dq = dt*D*[force.Vx force.Vy force.Vz]/(kB * T) + ...
        sqrt(2*dt*D)*randn(1,3);
    
    % Spatial
    dx = dq(1); dy = dq(2); dz = dq(3);
    
    % Displacement of the center of the bead 
    x(n+1) = x(n) + dx;
    y(n+1) = y(n) + dy;
    z(n+1) = z(n) + dz;
    
    % Plot every fifth iteration
    if mod(n,5)==0
        
        % Plot in 3D: incoming beams, bead, outgoing rays
        set (gcf,'CurrentAxes',h_3d);
        cla
        hold on
        title(['t = ' num2str(n*dt) ' s'])
               
        % Trajectory (blue line)
        plot3(x(1:n)*1e+6,y(1:n)*1e+6-1.5,z(1:n)*1e+6,'b');
        
        % Bead (green)
        bead.plot('scale', 1e+6, ...
            'FaceColor', [0 0.75 0], ...
            'EdgeColor', [0.0 0.0 0.0], ...
            'FaceAlpha', 0.2, ...
            'EdgeAlpha', 0.2 ...
            );
        
        
        % Beams: incoming rays, transmitted rays, outgoing rays
        
        % The rays are plotted on the basis of the "reduced" bunch rp
        scatt = bead.scattering(rp,1e-12,2);
        
        % First scattering event
        intpoin1 = Point(scatt(1).t.v.X,scatt(1).t.v.Y,scatt(1).t.v.Z);
        
        % Second scattering event
        intpoin2 = Point(scatt(2).t.v.X,scatt(2).t.v.Y,scatt(2).t.v.Z);
        
        % Transmitted rays: final coordinates
        fipoin = screen.intersectionpoint(scatt(2).t.v,2);
        
        % Nir (non intersecting rays) are determined on the basis of "starting points"
        nir = stpoin;
        nir.X(isfinite(intpoin2.X))=NaN;
        nir.Y(isfinite(intpoin2.Y))=NaN;
        nir.Z(isfinite(intpoin2.Z))=NaN;
        nire = nir;
        for iphi=(1:Nphi)
            nire.X(iphi,:)=o_x1-(nir.X(iphi,:)-o_x1);
            nire.X(iphi+Nphi,:)=o_x2-(nir.X(iphi+Nphi,:)-o_x2);
        end
        nire.Y=-nir.Y;
        nire.Z=-nir.Z;

        lin0 = SLine(nir,nire);
        lin1 = SLine(stpoin,intpoin1);
        lin2 = SLine(intpoin1,intpoin2);
        lin3 = SLine(intpoin2,fipoin);

        % Plot rays
        lin0.plot('scale',1e+6,'color',[1 0.5 0.5],'LineWidth',2);
        lin1.plot('scale',1e+6,'color','r','LineWidth',2);
        lin2.plot('scale',1e+6,'color',[1 1 1]*0.7,'LineWidth',2);
        lin3.plot('scale',1e+6,'color',[1 0.7 0.5],'LineWidth',2);
                
        % Plot x(t)
        set (gcf,'CurrentAxes',h_x);
        cla
        hold on
        % Plot lines indicating the position of the traps
        plot([1 n]*dt,[1 1]*o_x1*1e+6,'r');
        plot([1 n]*dt,[1 1]*o_x2*1e+6,'r');
        % Plot zero
        plot([1 n]*dt,[0 0],'color',[0.5 0.5 0.5]);
        % Plot trajectory
        plot((1:n)*dt,x(1:n)*1e+6,'k');
        
        % Plot y(t)
        set (gcf,'CurrentAxes',h_y);
        cla
        hold on
        % Plot zero
        plot([1 n]*dt,[0 0],'color',[0.5 0.5 0.5]);
        % Plot trajectory
        plot((1:n)*dt,y(1:n)*1e+6,'k');        
        
        % Plot z(t)
        set (gcf,'CurrentAxes',h_z);
        cla
        hold on
        % Plot zero
        plot([1 n]*dt,[0 0],'color',[0.5 0.5 0.5]);
        % Plot trajectory
        plot((1:n)*dt,z(1:n)*1e+6,'k');        
        
        drawnow()

    end
    
end