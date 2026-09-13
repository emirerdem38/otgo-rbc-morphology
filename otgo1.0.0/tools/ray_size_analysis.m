function [forces_torques, comp_times] = ray_size_analysis(shp, p, ray_sizes)

    % Initialize
    nm = 1.33;
    np = 1.38;
    
    focal = 5*1e-6; % Focal length [m] % previously 5 % 200
    NA = 1.3; %1.30; % Numerical aperture
    
    Ex0 = 1e+4; % x electric field [V/m]
    Ey0 = 1i*1e+4; % y electric field [V/m]
    w0 = 5e-6; % Beam waist [m]
    power = 5e-3; % power [W]
    
    scat = 10;
    
    % Random displacement and orientation generation
    seed = 42;
    rng(seed);
    
    c = Point(randn(1,1),randn(1,1),randn(1,1))*1e-6*0.25;
    
    particle = ParticleCELL(p,shp,c,Vector(0,0,0,1,0,0),Vector(0,0,0,0,1,0),Vector(0,0,0,0,0,1),'vol',-6,nm,np);
    particle = particle.xrotation(randn(1,1));
    particle = particle.yrotation(randn(1,1));
    particle = particle.zrotation(randn(1,1));
    
    %ray_sizes = [10 10; 20 20; 30 30; 40 40; 50 50; 60 60; 70 70; 80 80; 90 90; 100 100; 110 110; 120 120; 130 130; 140 140; 150 150; 160 160];
    %ray_sizes = [(15:5:160)' (15:5:160)'];
    
    forces_torques = zeros(6,size(ray_sizes,1));
    comp_times = NaN(size(ray_sizes,1),2);
    
    for i = 1:size(ray_sizes,1)
            
        Nphi = ray_sizes(i,1); % Azimuthal divisions
        Nr = ray_sizes(i,2); % Radial divisions
        %Nphi = 20;
        %Nr = 40;
        r = CreateBeam(nm, focal, NA, Ex0, Ey0, w0, Nphi, Nr, power);
        disp(numel(r))
    
        tic;
        r_vec = scattering(particle,r,0,scat-1);
        f = particle.force(r_vec,r);
        t = particle.torque(r_vec,r);
            
        force = [sum(f.Vx(isfinite(f.Vx)));
                     sum(f.Vy(isfinite(f.Vy)));
                     sum(f.Vz(isfinite(f.Vz)));
                    ];
        torque = [sum(t.Vx(isfinite(t.Vx)));
                      sum(t.Vy(isfinite(t.Vy)));
                      sum(t.Vz(isfinite(t.Vz)));
                     ];
    
        elapsedTime = toc;
        fprintf('Total time: %.2f seconds\n', elapsedTime)
    
        forces_torques(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
    
        comp_times(i,1) = numel(r);
        comp_times(i,2) = elapsedTime;
    end
end

