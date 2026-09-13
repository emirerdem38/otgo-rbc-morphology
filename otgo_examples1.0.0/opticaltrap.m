% OPTICALTRAP Optical trap and forward scattering
%
% Simulates the Brownian motion of a spehrical particle in an optical trap
% under the action of both optical forces and thermal forces. It also shows
% the forward scattered light that can be used to measure the particle
% position in three dimensions.
%
% Three subplots are shown:
% 1) on the left, the intensity of the incoming beam
% 2) at the center, the optically trapped particle with the incoming and
% first scattered rays
% 3) on the right, the intensity distribution of the forward scattered
% light on the back-focal-plane of the condenser
%
% See also Vector, Ray, BeamGauss, Spherical, ParticleSpherical.
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

% Spherical Particle
R = 1.0e-06; % [m]
np = 1.50; % Particle refractive index

% Focusing
f = 5e-6; % Focal length [m]
NA = 1.30; % Numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 5e-6; % Beam waist [m]
Nphi = 100; % Azimuthal divisions
Nr = 100; % Radial divisions
power = 0.25e-3; % power [W]

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K]
eta = 0.001; % Water viscosity [Pa*s]
dt = 5e-3; % timestep [s]
N = 1e+3; % number of steps

% Initial position
x = -0.6*R;
y = 0.0;
z = -0.8*R;

%% Initialization

% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);
% rotate the beam: propagation direction along the x axis
r.v = r.v.yrotation(pi/2);
r.pol = r.pol.yrotation(pi/2);

% Screen (sphere with radius f centered in the focal point)
rscreen = f;
screen = Spherical(Point(0,0,0),rscreen);

% Diffusion Constant of the Spherical Particle
D = kB*T / (6*pi*eta*R);

%% Figure Initialization

figure('Units', 'Pixels', 'Position', [100 100 940 320], 'Color', 'w')

%% Plot incoming beam

axes ('Units','Pixels','Position',[10  10  300   300],'Color','w');

% Levels of the color map
levels = 128;

hold on

% Intensity profile of the incoming beam in cartesian coordinates
[xc,yc] = Transform.Pol2Car(bg.phi,bg.r);
I = bg.intensity();

contourf([zeros(size(xc,1)+1,1),[xc;xc(1,:)]],...
    [zeros(size(yc,1)+1,1),[yc;yc(1,:)]],...
    [[I(:,1);I(1,1)],[I;I(1,:)]], ...
    levels)
colormap(ones(levels,3)-(0:1/(levels-1):1)'*[0 1 1])
shading flat

% Sketch reference axis
plot([max(max(xc)) min(min(xc))], [0 0], 'k')
plot([0 0], [max(max(yc)) min(min(yc))], 'k')

% Sketch iris
theta = (-1:0.01:1)*pi;
xcirc = max(max(xc))*cos(theta);
ycirc = max(max(xc))*sin(theta);
plot(xcirc, ycirc, 'k', 'LineWidth', 1);


axis equal tight off

drawnow()

%% Axes for central subplot
h_c = axes('Units', 'Pixels', 'Position', [320  10  300   300], 'Color', 'k');

hold on

xlim([-1 1]*1.1*L*1e+6)
ylim([-1 1]*1.1*L*1e+6)
zlim([-1 1]*1.1*L*1e+6)

axis equal off
view(0,0)

%% Axes for right subplot
h_r = axes('Units', 'Pixels', 'Position', [630  10  300   300], 'Color', 'w');

hold on

xlim([min(min(xc)) max(max(xc))])
ylim([min(min(yc)) max(max(yc))])

axis equal off

%% Brownian motion simulation

for n = 1:1:N
    
    % Display update message
    disp(['time = ' num2str(n*dt) ' s / ' num2str(N*dt) ' s'])
    
    % Spherical particle
    bead = ParticleSpherical(Point(x,y,z),R,nm,np);

    r_vec = scattering(bead,r);
    
    %% Optical force
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

    %% Figure update - Central subplot
    
    % Calculate scattered rays
    s = r_vec;
    
    % Outgoing rays
    s_out(1)=s(1).r;
    for i=2:size(s,2)
        s_out(i)=s(i).t;
    end
    
    % Selection rays propagating towards the screen, i.e., in the positive x direction
    s_out(1).P(find(s(1).r.v.Vx<=0))=NaN;
    for i=1:size(s,2)
        s_out(i).P(find(s(i).t.v.Vx<=0))=NaN;
    end
    
    % Intersection between with the screen
    for i=1:size(s,2)
        ipoint(i) = screen.intersectionpoint(s_out(i),2);
    end
    
    % Plot of the configuration: incoming focused beam - bead - scattered rays
    set(gcf,'CurrentAxes',h_c)
    
    cla
    
    % Plot the bead
    bead.plot('scale', 1e+6, ...
        'facecolor', [0 0.75 0], ...
        'edgecolor', [0 0 0], ...
        'facealpha', .2, ...
        'edgealpha', 0.1 ...
        );
    
    % Plot bead center
    cbpoint = Point(bead.sp.c.X,-0.9*L,bead.sp.c.Z);
    cbpoint.plot('scale', 1e+6);
    
    % Plot subsample of incoming rays
    for iphi = 10:10:Nphi
        for ir = 10:10:Nr
            
            % Incoming ray
            plot3([r.v.X(iphi,ir) s(1).r.v.X(iphi,ir)]*1e+6,...
                [r.v.Y(iphi,ir) s(1).r.v.Y(iphi,ir)]*1e+6,...
                [r.v.Z(iphi,ir) s(1).r.v.Z(iphi,ir)]*1e+6,...
                'color', [1 0 0], ...
                'LineWidth', 1 ...
                );
            
            % First transmitted ray (internal ray)
            plot3([s(1).r.v.X(iphi,ir) s(2).t.v.X(iphi,ir)]*1e+6,...
                [s(1).r.v.Y(iphi,ir) s(2).t.v.Y(iphi,ir)]*1e+6,...
                [s(1).r.v.Z(iphi,ir) s(2).t.v.Z(iphi,ir)]*1e+6,...
                'color', [0.7 0.7 0.7], ...
                'LineWidth', 1 ...
                );
            
            % Second transmitted ray (outgoing ray)
            plot3([s(2).t.v.X(iphi,ir) ipoint(2).X(iphi,ir)]*1e+6,...
                [s(2).t.v.Y(iphi,ir) ipoint(2).Y(iphi,ir)]*1e+6,...
                [s(2).t.v.Z(iphi,ir) ipoint(2).Z(iphi,ir)]*1e+6,...
                'color', [1 0.3 0.5], ...
                'LineWidth', 1 ...
                );
            
        end
    end
    
    % Plot objective and screen
    thmax = pi/180*80;
    theta = (-1:0.1:1)*thmax;
    xcirc = f*cos(theta)*1e+6;
    ycirc = f*sin(theta)*1e+6;
    plot3(xcirc, 0*ycirc-0.9*L*1e+6, ycirc, 'b', 'LineWidth', 3)
    plot3(-xcirc, 0*ycirc-0.9*L*1e+6, ycirc, 'b', 'LineWidth', 3)
    
    % Plot reference axis
    plot3([-1.05 1.05]*L*1e+6, [-0.9 -0.9]*L*1e+6, [0 0], 'k--')
    plot3([0 0], [-0.9 -0.9]*L*1e+6, [-1 1]*L*1e+6, 'k--')
    
    drawnow()
    
    %% Figure update - Right subplot
    
    % Forward scattering (rays scattered on the screen)
    
    % The screen is divided into Nr concentrical equidistant circular rings
    % and each ring is dividen in Nphi equivalent segments.
    % Each ray is assigned to the proper sector of the corresponding
    % intersection point ipoint(i).
    
    % Find the index for the circular ring.
    for i = 1:1:size(s,2)
        % rho = transversal distance from the x axis
        rho(i,:,:) = real(sqrt(ipoint(i).Y.^2 + ipoint(i).Z.^2));
    end
    ic_rho = ceil(rho./(L/Nr)); % interval in radial divisions to which rho belongs
    
    % Find index for the angular segment in the ring.
    for i = 1:1:size(s,2)
        for ir = 1:1:Nr
            for iphi = 1:1:2*Nphi
                % phi = angle on the transversal plane of the intersection point with respect to the x axis
                phi(i,iphi,ir) = atan2( ipoint(i).Z(iphi,ir), ipoint(i).Y(iphi,ir) ) ...
                    + abs( sign( ipoint(i).Z(iphi,ir) ) ).*( 1-sign( ipoint(i).Z(iphi,ir) ) ) ...
                    * pi;
            end
        end
    end
    ic_phi = ceil(phi./(2*pi/Nphi)); % Interval in angular divisions to which phi belongs
    
    % Power on the screen - polar grid
    Ptot=zeros(Nphi,Nr); % Total power for each division (radial and angular)
    for i = 1:1:size(s,2)
        for ir = 1:1:Nr
            for iphi = 1:1:2*Nphi
                % Select only the rays falling in the correct region of the screen
                if (ic_phi(i,iphi,ir)<=Nphi) && (ic_rho(i,iphi,ir)<=Nr) && isfinite(s_out(i).P(iphi,ir))
                    Ptot(ic_phi(i,iphi,ir),ic_rho(i,iphi,ir)) = ...
                        Ptot(ic_phi(i,iphi,ir),ic_rho(i,iphi,ir)) + s_out(i).P(iphi,ir);
                end
            end
        end
    end
    
    % Power on the screen - Cartesian grid
    nupix = 10; % Number of divisions (semiaxis)
    dx = (max(max(xc))-min(min(xc)))/2/nupix; % Half length of the pixel (horizontal)
    dy = (max(max(yc))-min(min(yc)))/2/nupix; % Half length of the pixel (vertical)
    
    xp_c = (-(nupix):1:(nupix))*2*dx;
    yp_c = (-(nupix):1:(nupix))*2*dy;
    P2_p = zeros(2*nupix+1,2*nupix+1); % Power falling in each division
    
    for iphi=1:Nphi
        for ir=1:Nr
            iy=min(find(xc(iphi,ir)<(xp_c+dx)));
            ix=min(find(yc(iphi,ir)<(yp_c+dy)));
            P2_p(ix,iy) = P2_p(ix,iy) + Ptot(iphi,ir);
        end
    end
    
    % Plot of the image on the screen
    P2_p_max = max(max(P2_p));
    
    set(gcf,'CurrentAxes',h_r)
    
    cla
    
    contourf(xp_c, yp_c, P2_p, levels);
    colormap(ones(levels,3)-(0:1/(levels-1):1)'*[0 1 1])
    shading flat
    
    % Sketch reference axis
    plot([max(max(xc)) min(min(xc))],[0 0],'k');
    plot([0 0],[max(max(yc)) min(min(yc))],'k');
    
    % Sketch iris
    theta = (-1:0.01:1)*pi;
    xcirc = max(max(xc))*cos(theta);
    ycirc = max(max(xc))*sin(theta);
    plot(xcirc, ycirc, 'k', 'LineWidth', 1);
    
    drawnow()
    
end