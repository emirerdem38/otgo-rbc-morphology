function [] = Trap(nm,f,NA,Ex0,Ey0,w0,Nphi,Nr,power,par,data,limits,beam_direction,saveFilePath,varargin)

if ~isempty(varargin)
end

try
    Point(0,0,0);
catch
    run("D:\Bilkent\CFDProje\OTGO\otgo1.0.0\OTGO.m")
end

r = CreateBeam(nm,f,NA,Ex0,Ey0,w0,Nphi,Nr,power);
scat = 10;

mkdir(saveFilePath)

or1 = Vector(0,0,0,1,0,0);
or2 = Vector(0,0,0,0,1,0);
or3 = Vector(0,0,0,0,0,1);
particle = ParticleCSVCell(Point(0,0,0),or1,or2,or3,"vol",-6,data,par.nm,par.np);

ParLoop(12);

try 
    x_bottom = limits(1,1); x_top = limits(1,2); Nx = limits(1,3);
    y_bottom = limits(2,1); y_top = limits(2,2); Ny = limits(2,3);
    z_bottom = limits(3,1); z_top = limits(3,2); Nz = limits(3,3);
    x_bottom_rot = limits(4,1); x_top_rot = limits(4,2); Nx_rot = limits(4,3);
    y_bottom_rot = limits(5,1); y_top_rot = limits(5,2); Ny_rot = limits(5,3);
    z_bottom_rot = limits(6,1); z_top_rot = limits(6,2); Nz_rot = limits(6,3);
catch
end

%% X Rotational Analysis

% Nx_rot= 100;
% x_bottom_rot = -pi/2;
% x_top_rot = pi/2;
divx_rot = linspace(x_bottom_rot,x_top_rot,Nx_rot);        

[forces_torquesx_rot, timex_rot] = Solve(r,scat,power,particle,data,Nx_rot,x_bottom_rot,x_top_rot,4,saveFilePath);
        
image_titleX_rot = sprintf(['Deformed Cell X Discretization (ROT) \n N = %d R = %dx%d P = %gW scat = %d \n' ...
'Yvec:(%.3f,%.3f,%.3f), Zvec:(%.3f,%.3f,%.3f) \n Center: (%.3f,%.3f,%.3f) \n Step 1'], Nx_rot, size(r,1), size(r,2), power, ...
scat, particle.cell.or2.Vx,particle.cell.or2.Vy,particle.cell.or2.Vz, ...
particle.cell.or3.Vx,particle.cell.or3.Vy, particle.cell.or3.Vz, particle.cell.c.X, particle.cell.c.Y, particle.cell.c.Z);

modifiedFilePath = sprintf('X%d(%dx%d)(scat%d)(P%g)',Nx_rot,size(r,1),size(r,2),scat,power);
FilePathX_rot = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleX_rot,divx_rot,forces_torquesx_rot,'X_rot',FilePathX_rot)

if beam_direction == "x"
    intervalsx_rot = 0;
else
    intervalsx_rot = SlopeIntervals(divx_rot,forces_torquesx_rot(4,:));
end
if ~isempty(intervalsx_rot)
    for i = 1:size(intervalsx_rot,1)
        image_titleX_rot_slope = sprintf('%s Slope (%d)',image_titleX_rot,i);
        plot_graph(image_titleX_rot_slope,divx_rot,forces_torquesx_rot,'X_rot',FilePathX_rot,intervalsx_rot(i,1:2),4,i)
    end
end

if ~isempty(intervalsx_rot)
    slopes_xrot = zeros(1,size(intervalsx_rot,1));
    for i = 1:size(intervalsx_rot,1)
        slopes_xrot(i) = calculate_slope(forces_torquesx_rot,divx_rot,4,intervalsx_rot(i,1),intervalsx_rot(i,2));
    end

    [~, max_slope_xrot] = max(abs(slopes_xrot));
    cross_xrot = intervalsx_rot(max_slope_xrot,3);
end

if exist('cross_xrot','var')
    if divx_rot(cross_xrot) < max(divx_rot)/10 && divx_rot(cross_xrot) > -max(divx_rot)/10
        cross_xrot_val = 0;
    else
        cross_xrot_val = divx_rot(cross_xrot);
    end
else
    cross_xrot_val = 0;
end

or1 = or1.xrotation(cross_xrot_val); or2 = or2.xrotation(cross_xrot_val); or3 = or3.xrotation(cross_xrot_val); 
particle.cell.or1 = or1; particle.cell.or2 = or2; particle.cell.or3 = or3;

%% Y Rotational Analysis

% Ny_rot= 100;
% y_bottom_rot = -pi/2;
% y_top_rot = pi/2;
divy_rot = linspace(y_bottom_rot,y_top_rot,Ny_rot);        

[forces_torquesy_rot, timey_rot] = Solve(r,scat,power,particle,data,Ny_rot,y_bottom_rot,y_top_rot,5,saveFilePath);
        
image_titleY_rot = sprintf(['Deformed Cell Y Discretization (ROT)\n N = %d R = %dx%d P = %gW scat = %d \n' ...
            'Xvec:(%.3f,%.3f,%.3f), Zvec:(%.3f,%.3f,%.3f) \n Center: (%.3f,%.3f,%.3f) \n Step 2'], Ny_rot, size(r,1), size(r,2), power, ...
            scat, particle.cell.or1.Vx,particle.cell.or1.Vy,particle.cell.or1.Vz, ...
            particle.cell.or3.Vx,particle.cell.or3.Vy,particle.cell.or3.Vz,particle.cell.c.X,particle.cell.c.Y,particle.cell.c.Z);
        

modifiedFilePath = sprintf('Y%d(%dx%d)(scat%d)(P%g)',Ny_rot,size(r,1),size(r,2),scat,power);
FilePathY_rot = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleY_rot,divy_rot,forces_torquesy_rot,'Y_rot',FilePathY_rot)

if beam_direction == "y"
    intervalsy_rot = 0;
else
    intervalsy_rot = SlopeIntervals(divy_rot,forces_torquesy_rot(5,:));
end
if ~isempty(intervalsy_rot)
    for i = 1:size(intervalsy_rot,1)
        image_titleY_rot_slope = sprintf('%s Slope (%d)',image_titleY_rot,i);
        plot_graph(image_titleY_rot_slope,divy_rot,forces_torquesy_rot,'Y_rot',FilePathY_rot,intervalsy_rot(i,1:2),5,i)
    end
end

if ~isempty(intervalsy_rot)
    slopes_yrot = zeros(1,size(intervalsy_rot,1));
    for i = 1:size(intervalsy_rot,1)
        slopes_yrot(i) = calculate_slope(forces_torquesy_rot,divy_rot,5,intervalsy_rot(i,1),intervalsy_rot(i,2));
    end

    [~, max_slope_yrot] = max(abs(slopes_yrot));
    cross_yrot = intervalsy_rot(max_slope_yrot,3);
end

if exist('cross_yrot','var')
    if divy_rot(cross_yrot) < max(divy_rot)/10 && divy_rot(cross_yrot) > -max(divy_rot)/10
        cross_yrot_val = 0;
    else
        cross_yrot_val = divy_rot(cross_yrot);
    end
else
    cross_yrot_val = 0;
end

or1 = or1.xrotation(cross_yrot_val); or2 = or2.xrotation(cross_yrot_val); or3 = or3.xrotation(cross_yrot_val); 
particle.cell.or1 = or1; particle.cell.or2 = or2; particle.cell.or3 = or3;

%% Z Rotational Analysis

% Nz_rot= 100;
% z_bottom_rot = -pi/2;
% z_top_rot = pi/2;
divz_rot = linspace(z_bottom_rot,z_top_rot,Nz_rot);        

[forces_torquesz_rot, timez_rot] = Solve(r,scat,power,particle,data,Nz_rot,z_bottom_rot,z_top_rot,6,saveFilePath);
        
image_titleZ_rot = sprintf(['Deformed Cell Z Discretization (ROT) \n N = %d R = %dx%d P = %gW scat = %d \n' ...
            'Xvec:(%.3f,%.3f,%.3f), Yvec:(%.3f,%.3f,%.3f) \n Center: (%.3f,%.3f,%.3f) \n Step 3'], Nz_rot, size(r,1), size(r,2), power, ...
            scat, particle.cell.or1.Vx,particle.cell.or1.Vy,particle.cell.or1.Vz, ...
            particle.cell.or2.Vx,particle.cell.or2.Vy, particle.cell.or2.Vz, particle.cell.c.X, particle.cell.c.Y, particle.cell.c.Z);
       
modifiedFilePath = sprintf('Z%d(%dx%d)(scat%d)(P%g)',Nz_rot,size(r,1),size(r,2),scat,power);
FilePathZ_rot = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleZ_rot,divz_rot,forces_torquesz_rot,'Z_rot',FilePathZ_rot)

if beam_direction == "z"
    intervalsz_rot = [];
else
    intervalsz_rot = SlopeIntervals(divz_rot,forces_torquesz_rot(6,:));
end

if ~isempty(intervalsz_rot)
    for i = 1:size(intervalsz_rot,1)
        image_titleZ_rot_slope = sprintf('%s Slope (%d)',image_titleZ_rot,i);
        plot_graph(image_titleZ_rot_slope,divz_rot,forces_torquesz_rot,'Z_rot',FilePathZ_rot,intervalsz_rot(i,1:2),6,i)
    end
end

if ~isempty(intervalsz_rot)
    slopes_zrot = zeros(1,size(intervalsz_rot,1));
    for i = 1:size(intervalsz_rot,1)
        slopes_zrot(i) = calculate_slope(forces_torquesz_rot,divz_rot,6,intervalsz_rot(i,1),intervalsz_rot(i,2));
    end

    [~, max_slope_zrot] = max(abs(slopes_zrot));
    cross_zrot = intervalsz_rot(max_slope_zrot,3);
end

if exist('cross_zrot','var')
    if divz_rot(cross_zrot) < max(divz_rot)/10 && divz_rot(cross_zrot) > -max(divz_rot)/10
        cross_zrot_val = 0;
    else
        cross_zrot_val = divz_rot(cross_zrot);
    end
else
    cross_zrot_val = 0;
end

or1 = or1.xrotation(cross_zrot_val); or2 = or2.xrotation(cross_zrot_val); or3 = or3.xrotation(cross_zrot_val); 
particle.cell.or1 = or1; particle.cell.or2 = or2; particle.cell.or3 = or3;

%% X Displacement Analysis

% Nx = 90;
% x_bottom = -0.3e-6;
% x_top = 0.3e-6;
divx = linspace(x_bottom,x_top,Nx);

[forces_torquesx, time_x] = Solve(r,scat,power,particle,data,Nx,x_bottom,x_top,1,saveFilePath);

image_titleX = sprintf(['Deformed Cell X Discretization \n N = %d R = %dx%d P = %gW scat = %d \n' ...
            'Xvec:(%.3f,%.3f,%.3f), Yvec:(%.3f,%.3f,%.3f), Zvec:(%.3f,%.3f,%.3f) \n Y:%.3f, Z:%.3f \n Step 4'], Nx, size(r,1), size(r,2), power, ...
            scat, particle.cell.or1.Vx,particle.cell.or1.Vy,particle.cell.or1.Vz,particle.cell.or2.Vx,particle.cell.or2.Vy,particle.cell.or2.Vz, ...
            particle.cell.or3.Vx,particle.cell.or3.Vy,particle.cell.or3.Vz,particle.cell.c.Y,particle.cell.c.Z);

modifiedFilePath = sprintf('X%d(%dx%d)(scat%d)(P%g)',Nx,size(r,1),size(r,2),scat,power);
FilePathX = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleX,divx,forces_torquesx,'X_disp',FilePathX)

intervalsx = SlopeIntervals(divx,forces_torquesx(1,:));
if ~isempty(intervalsx)
    for i = 1:size(intervalsx,1)
        image_titleX_slope = sprintf('%s Slope (%d)',image_titleX,i);
        plot_graph(image_titleX_slope,divx,forces_torquesx,'X_disp',FilePathX,intervalsx(i,1:2),1,i)
    end
end

if ~isempty(intervalsx)
    slopes_x = zeros(1,size(intervalsx,1));
    for i = 1:size(intervalsx,1)
        slopes_x(i) = calculate_slope(forces_torquesx,divx,1,intervalsx(i,1),intervalsx(i,2));
    end

    [~, max_slope_x] = max(abs(slopes_x));
    cross_x = intervalsx(max_slope_x,3);
end

if exist('cross_x','var')
    if divx(cross_x) < max(divx)/10 && divx(cross_x) > -max(divx)/10
        cross_xval = 0;
    else
        cross_xval = divx(cross_x);
    end
else
    cross_xval = 0;
end

particle.cell.c.X = cross_xval;

%% Y Displacement Analysis

% Ny = 300;
% y_bottom = -1e-6;
% y_top = 1e-6;
divy = linspace(y_bottom,y_top,Ny);

[forces_torquesy, time_y] = Solve(r,scat,power,particle,data,Ny,y_bottom,y_top,2,saveFilePath);

image_titleY = sprintf(['Deformed Cell Y Discretization \n N = %d R = %dx%d P = %gW scat = %d \n' ...
            'Xvec:(%.3f,%.3f,%.3f), Yvec:(%.3f,%.3f,%.3f), Zvec:(%.3f,%.3f,%.3f) \n X:%.3f, Z:%.3f \n Step 5'], Ny, size(r,1), size(r,2), power, ...
            scat, particle.cell.or1.Vx,particle.cell.or1.Vy,particle.cell.or1.Vz,particle.cell.or2.Vx,particle.cell.or2.Vy,particle.cell.or2.Vz, ...
            particle.cell.or3.Vx,particle.cell.or3.Vy,particle.cell.or3.Vz,particle.cell.c.X,particle.cell.c.Z);

modifiedFilePath = sprintf('Y%d(%dx%d)(scat%d)(P%g)',Ny,size(r,1),size(r,2),scat,power);
FilePathY = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleY,divy,forces_torquesy,'Y_disp',FilePathY)

intervalsy = SlopeIntervals(divy,forces_torquesy(2,:));
if ~isempty(intervalsy)
    for i = 1:size(intervalsy,1)
        image_titleY_slope = sprintf('%s Slope (%d)',image_titleY,i);
        plot_graph(image_titleY_slope,divy,forces_torquesy,'Y_disp',FilePathY,intervalsy(i,1:2),1,i)
    end
end

if ~isempty(intervalsy)
    slopes_y = zeros(1,size(intervalsy,1));
    for i = 1:size(intervalsy,1)
        slopes_y(i) = calculate_slope(forces_torquesy,divy,2,intervalsy(i,1),intervalsy(i,2));
    end

    [~, max_slope_y] = max(abs(slopes_y));
    cross_y = intervalsy(max_slope_y,3);
end

if exist('cross_y','var')
    if divy(cross_y) < max(divy)/10 && divy(cross_y) > -max(divy)/10
        cross_yval = 0;
    else
        cross_yval = divy(cross_y);
    end
else
    cross_yval = 0;
end

particle.cell.c.Y = cross_yval;

%% Z Displacement Analysis

% Nz = 300;
% z_bottom = -1e-6;
% z_top = 1e-6;
divz = linspace(z_bottom,z_top,Nz);

[forces_torquesz, time_z] = Solve(r,scat,power,particle,data,Nz,z_bottom,z_top,3,saveFilePath);

image_titleZ = sprintf(['Deformed Cell Z Discretization \n N = %d R = %dx%d P = %gW scat = %d \n' ...
            'Xvec:(%.3f,%.3f,%.3f), Yvec:(%.3f,%.3f,%.3f), Zvec:(%.3f,%.3f,%.3f) \n X:%.3f, Y:%.3f \n Step 6'], Ny, size(r,1), size(r,2), power, ...
            scat, particle.cell.or1.Vx,particle.cell.or1.Vy,particle.cell.or1.Vz,particle.cell.or2.Vx,particle.cell.or2.Vy,particle.cell.or2.Vz, ...
            particle.cell.or3.Vx,particle.cell.or3.Vy,particle.cell.or3.Vz,particle.cell.c.X,particle.cell.c.Y);

modifiedFilePath = sprintf('Z%d(%dx%d)(scat%d)(P%g)',Nz,size(r,1),size(r,2),scat,power);
FilePathZ = fullfile(saveFilePath,modifiedFilePath);

plot_graph(image_titleZ,divz,forces_torquesz,'Z_disp',FilePathZ)

intervalsz = SlopeIntervals(divz,forces_torquesz(3,:));
if ~isempty(intervalsz)
    for i = 1:size(intervalsz,1)
        image_titleZ_slope = sprintf('%s Slope (%d)',image_titleZ,i);
        plot_graph(image_titleZ_slope,divz,forces_torquesz,'Z_disp',FilePathZ,intervalsz(i,1:2),1,i)
    end
end

if ~isempty(intervalsz)
    slopes_z = zeros(1,size(intervalsz,1));
    for i = 1:size(intervalsz,1)
        slopes_z(i) = calculate_slope(forces_torquesz,divz,3,intervalsz(i,1),intervalsz(i,2));
    end

    [~, max_slope_z] = max(abs(slopes_z));
    cross_z = intervalsz(max_slope_z,3);
end

if exist('cross_z','var')
    if divz(cross_z) < max(divz)/10 && divz(cross_z) > -max(divz)/10
        cross_zval = 0;
    else
        cross_zval = divz(cross_z);
    end
else
    cross_zval = 0;
end

particle.cell.c.Z = cross_zval;

parameterTable = table(Nx, Ny, Nz, Nx_rot, Ny_rot, Nz_rot, particle.cell.or1.X, particle.cell.or1.Y, ...
    particle.cell.or1.Z, particle.cell.or1.Vx, particle.cell.or1.Vy, particle.cell.or1.Vz, particle.cell.or2.X, ...
    particle.cell.or2.Y, particle.cell.or2.Z, particle.cell.or2.Vx, particle.cell.or2.Vy, particle.cell.or2.Vz, ...
    particle.cell.or3.X, particle.cell.or3.Y, particle.cell.or3.Z, particle.cell.or3.Vx, particle.cell.or3.Vy, ...
    particle.cell.or3.Vz, x_bottom, x_top, y_bottom, y_top, z_bottom, z_top, x_bottom_rot, x_top_rot, y_bottom_rot, y_top_rot, ...
    z_bottom_rot, z_top_rot, particle.cell.c.X, particle.cell.c.Y, particle.cell.c.Z, nm, particle.np, f, NA, f*NA/nm, Ex0, Ey0, w0, Nphi, Nr, power, scat, ...
    'VariableNames', {'Nx', 'Ny', 'Nz', 'Nx_rot', 'Ny_rot', 'Nz_rot', 'or1.X', 'or1.Y', 'or1.Z', 'or1.Vx', ...
    'or1.Vy', 'or1.Vz', 'or2.X', 'or2.Y', 'or2.Z', 'or2.Vx', 'or2.Vy', 'or2.Vz', 'or3.X', 'or3.Y', 'or3.Z', ...
    'or3.Vx', 'or3.Vy', 'or3.Vz', 'x_bottom', 'x_top', 'y_bottom', 'y_top', 'z_bottom', 'z_top', 'x_bottom_rot', ...
    'x_top_rot', 'y_bottom_rot', 'y_top_rot', 'z_bottom_rot', 'z_top_rot', 'center.X', 'center.Y', 'center.Z', 'nm' ...
    'np', 'f', 'NA', 'L', 'Ex0', 'Ey0', 'w0', 'Nphi', 'Nr', 'power','scatterings'});

writetable(parameterTable,fullfile(saveFilePath,'parameters.csv'))

%% Report Analysis

fileID = fopen(fullfile(saveFilePath,'Report.txt'), 'w');
if fileID == -1; error('Error creating file: %s', saveFilePath); end

fprintf(fileID, 'Analysis Report on CSV Loaded Cell\n');
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Optical Parametes: \n \n');
fprintf(fileID, 'Medium Refractive Index: %d\n', nm);
fprintf(fileID, 'Focal Length: %d\n', f);
fprintf(fileID, 'Numerical Aperture (NA): %d\n', NA);
fprintf(fileID, 'Focal Length: %d\n', f);
fprintf(fileID, 'Beam Waist: %g\n', w0);
fprintf(fileID, 'X Electric Field Amplitude: %g\n', Ex0);
fprintf(fileID, 'Y Electric Field Amplitude: %g\n', Ey0);
fprintf(fileID, 'Azimuthal Divisions: %d \n', Nphi);
fprintf(fileID, 'Radial Divisions: %d \n', Nr);
fprintf(fileID, 'Beam Power: %.3f W\n', power);
fprintf(fileID, 'Number of Considered Scattering Events: %d\n', scat);
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Particle Parametes: \n \n');
if particle.cell.vol == "vol"
    fprintf(fileID, 'Particle Volume: Default Volume is Used \n');
else
    fprintf(fileID, 'Particle Volume: %.3fe-18 m^3 \n', particle.cell.vol);
end
fprintf(fileID, 'Particle Refractive Index: %.3f\n', particle.np);
fprintf(fileID, 'Loaded Data Path: %s \n', data);
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Analysis Parametes: \n \n');
fprintf(fileID, 'X Discretization: From %g To %g meters (%d Points Sampled) \n', x_bottom, x_top, Nx);
fprintf(fileID, 'Y Discretization: From %g To %g meters (%d Points Sampled) \n', y_bottom, y_top, Ny);
fprintf(fileID, 'Z Discretization: From %g To %g meters (%d Points Sampled) \n', z_bottom, z_top, Nz);
fprintf(fileID, 'X Rotational Discretization: From %g To %g radians (%d Points Sampled) \n', x_bottom_rot, x_top_rot, Nx_rot);
fprintf(fileID, 'Y Rotational Discretization: From %g To %g radians (%d Points Sampled) \n', y_bottom_rot, y_top_rot, Ny_rot);
fprintf(fileID, 'Z Rotational Discretization: From %g To %g radians (%d Points Sampled) \n', z_bottom_rot, z_top_rot, Nz_rot);
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Analysis Results (Intervals For Line Fitting is Program Controlled, Manually Determine Intervals For More Accurate Slope Values): \n');

fprintf(fileID, '\nDisplacement in X Direction - Intervals of Negative Slope: \n');
if ~isempty(intervalsx)
    for i = 1:size(intervalsx,1)
        if i == max_slope_x
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsx(i,1), intervalsx(i,2), cross_xval, intervalsx(i,3), slopes_x(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsx(i,1), intervalsx(i,2), cross_xval, intervalsx(i,3), slopes_x(i));
        end
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end

fprintf(fileID, '\nDisplacement in Y Direction - Intervals of Negative Slope: \n');
if ~isempty(intervalsy)
    for i = 1:size(intervalsy,1)
        if i == max_slope_y
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsy(i,1), intervalsy(i,2), cross_yval, intervalsy(i,3), slopes_y(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsy(i,1), intervalsy(i,2), cross_yval, intervalsy(i,3), slopes_y(i));
        end       
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end

fprintf(fileID, '\nDisplacement in Z Direction - Intervals of Negative Slope: \n');
if ~isempty(intervalsz)
    for i = 1:size(intervalsz,1)
        if i == max_slope_z
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsz(i,1), intervalsz(i,2), cross_zval, intervalsz(i,3), slopes_z(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsz(i,1), intervalsz(i,2), cross_zval, intervalsz(i,3), slopes_z(i));
        end
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end

if beam_direction == "x"
    fprintf(fileID, '\nRotation in X Direction (Beam Direction) - Intervals of Negative Slope: \n');
else
    fprintf(fileID, '\nRotation in X Direction - Intervals of Negative Slope: \n');
end
if ~isempty(intervalsx_rot)
    for i = 1:size(intervalsx_rot,1)
        if i == max_slope_xrot
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsx_rot(i,1), intervalsx_rot(i,2), cross_xrot_val, intervalsx_rot(i,3), slopes_xrot(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsx_rot(i,1), intervalsx_rot(i,2), cross_xrot_val, intervalsx_rot(i,3), slopes_xrot(i));
        end
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end

if beam_direction == "y"
    fprintf(fileID, '\nRotation in Y Direction (Beam Direction) - Intervals of Negative Slope: \n');
else
    fprintf(fileID, '\nRotation in Y Direction - Intervals of Negative Slope: \n');
end
if ~isempty(intervalsy_rot)
    for i = 1:size(intervalsy_rot,1)
        if i == max_slope_yrot
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsy_rot(i,1), intervalsy_rot(i,2), cross_yrot_val, intervalsy_rot(i,3), slopes_yrot(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsy_rot(i,1), intervalsy_rot(i,2), cross_yrot_val, intervalsy_rot(i,3), slopes_yrot(i));
        end
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end

if beam_direction == "z"
    fprintf(fileID, '\nRotation in Z Direction (Beam Direction) - Intervals of Negative Slope: \n');
else
    fprintf(fileID, '\nRotation in Z Direction - Intervals of Negative Slope: \n');
end
if ~isempty(intervalsz_rot)
    for i = 1:size(intervalsz_rot,1)
        if i == max_slope_zrot
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g (MAX SLOPE) \n', i, intervalsz_rot(i,1), intervalsz_rot(i,2), cross_zrot_val, intervalsz_rot(i,3), slopes_zrot(i));
        else
            fprintf(fileID, 'Interval %d: From %g To %g, Zero Cross at %g (Index %d) with slope %g \n', i, intervalsz_rot(i,1), intervalsz_rot(i,2), cross_zrot_val, intervalsz_rot(i,3), slopes_zrot(i));
        end
    end
else
    fprintf(fileID, 'No Negative Slopes Can Be Obtained For This Case! \n');
end
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Effective Trapping Position (May Not Be Effectively Trapped In One or More Directions): \n \n');
fprintf(fileID, 'Particle Center: (%g,%g,%g) \n', particle.cell.c.X, particle.cell.c.Y, particle.cell.c.Z);
fprintf(fileID, 'Particle X Orientation: (%.3f,%.3f,%.3f), Effective When Oriented %.3f Radians \n', particle.cell.or1.Vx, particle.cell.or1.Vy, particle.cell.or1.Vz, cross_xrot_val);
fprintf(fileID, 'Particle Y Orientation: (%.3f,%.3f,%.3f), Effective When Oriented %.3f Radians \n', particle.cell.or2.Vx, particle.cell.or2.Vy, particle.cell.or2.Vz, cross_yrot_val);
fprintf(fileID, 'Particle Z Orientation: (%.3f,%.3f,%.3f), Effective When Oriented %.3f Radians \n', particle.cell.or3.Vx, particle.cell.or3.Vy, particle.cell.or3.Vz, cross_zrot_val);
fprintf(fileID, '------------------------------------------------- \n');
fprintf(fileID, 'Time Analysis: \n \n');
fprintf(fileID, 'X Displacement Analysis With %d Discretization Points Took %.2f Seconds To Complete. \n', Nx, time_x);
fprintf(fileID, 'Y Displacement Analysis With %d Discretization Points Took %.2f Seconds To Complete. \n', Ny, time_y);
fprintf(fileID, 'Z Displacement Analysis With %d Discretization Points Took %.2f Seconds To Complete. \n', Nz, time_z);
fprintf(fileID, 'X Rotational Analysis With %d Discretization Points Took %.2f To Seconds Complete. \n', Nx_rot, timex_rot);
fprintf(fileID, 'Y Rotational Analysis With %d Discretization Points Took %.2f To Seconds Complete. \n', Ny_rot, timey_rot);
fprintf(fileID, 'Z Rotational Analysis With %d Discretization Points Took %.2f To Seconds Complete. \n\n', Nz_rot, timez_rot);
fprintf(fileID, 'Total Time: %.2f Seconds', (time_x+time_y+time_z+timex_rot+timey_rot+timez_rot))
fclose(fileID);

fprintf('Analysis is completed, report is successfully saved to %s\n', fullfile(saveFilePath,'Report.txt'));

end

