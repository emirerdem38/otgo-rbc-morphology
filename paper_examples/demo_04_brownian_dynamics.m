%% demo_04_brownian_dynamics.m
% Short overdamped Brownian-dynamics trajectory in a single-beam trap.
%
% Requires otgo1.0.0/diffusion_tensors.mat (keys: D_dumbbell, D_fourBump,
% D_eightBump, D_neck, D_oblate75). Diffusion tensors are scaled by 1e-7 as
% in the paper workflow.
%
% This is a short demo (Nsteps = 200). Paper production runs used longer
% trajectories (thousands of steps, dt = 2 ms).

clear; clc; close all;
setup_paths;
P = paper_optical_params();

%% User settings
shp = 'dumbbell';
dt = 2e-3;          % s
Nsteps = 30;        % short demo (paper runs used thousands of steps)
seed = 42;
fast_demo = true;   % true -> p=32; false -> paper production p

p = P.p_of(shp);
if fast_demo
    p = min(p, 32);
end
kB = PhysConst.kB;
T = 293;

beam = CreateBeam(P.nm, P.f, P.NA, P.Ex0, P.Ey0, P.w0, P.Nphi, P.Nr, P.power);

% Diffusion tensor
Dfile = fullfile(fileparts(mfilename('fullpath')), '..', 'otgo1.0.0', 'diffusion_tensors.mat');
S = load(Dfile);
Dname = ['D_' shp];
if ~isfield(S, Dname)
    error('diffusion_tensors.mat has no field %s. Available: %s', ...
        Dname, strjoin(fieldnames(S), ', '));
end
D = S.(Dname) * 1e-7;
Lchol = chol(D, 'lower');

rng(seed);
rotation_angles = (pi / 2) * rand(1, 3);

particle = make_particle(shp, p, Point(0, 0, 0), [], [], [], P);
particle.cell = particle.cell.xrotation(rotation_angles(1));
particle.cell = particle.cell.yrotation(rotation_angles(2));
particle.cell = particle.cell.zrotation(rotation_angles(3));

Mpl = [particle.cell.or1.Vx, particle.cell.or2.Vx, particle.cell.or3.Vx; ...
       particle.cell.or1.Vy, particle.cell.or2.Vy, particle.cell.or3.Vy; ...
       particle.cell.or1.Vz, particle.cell.or2.Vz, particle.cell.or3.Vz];
Mlp = Mpl.';

centers = zeros(3, Nsteps);
fprintf('Brownian dynamics: %s, N = %d, dt = %.3f s\n', shp, Nsteps, dt);
tic;
for n = 1:Nsteps
    [F_lab, T_lab] = sum_force_torque(particle, beam, P.scat);
    f_p = Mlp * F_lab;
    t_p = Mlp * T_lab;

    whiteNoise = Lchol * randn(6, 1);
    dq = ((D * dt) / (kB * T)) * [f_p; t_p] + sqrt(2 * dt) * whiteNoise;

    new_c = [particle.cell.c.X; particle.cell.c.Y; particle.cell.c.Z] + Mpl * dq(1:3);
    particle.cell.c = Point(new_c(1), new_c(2), new_c(3));

    particle = particle.zrotation(dq(4));
    particle = particle.yrotation(dq(5));
    particle = particle.xrotation(dq(6));

    Mpl = [particle.cell.or1.Vx, particle.cell.or2.Vx, particle.cell.or3.Vx; ...
           particle.cell.or1.Vy, particle.cell.or2.Vy, particle.cell.or3.Vy; ...
           particle.cell.or1.Vz, particle.cell.or2.Vz, particle.cell.or3.Vz];
    Mlp = Mpl.';

    centers(:, n) = [particle.cell.c.X; particle.cell.c.Y; particle.cell.c.Z];
    if mod(n, 10) == 0 || n == 1 || n == Nsteps
        r3 = norm(centers(:, n));
        fprintf('  step %d/%d  r = %.3f um\n', n, Nsteps, r3 * 1e6);
    end
end
fprintf('  done in %.1f s\n', toc);

t = (1:Nsteps) * dt;
r = sqrt(sum(centers.^2, 1));

figure('Color', 'w', 'Name', sprintf('%s BD demo', shp));
subplot(1, 2, 1);
plot(t, r * 1e6, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('r = |c| [\mum]');
title(sprintf('%s: center--focus distance', shp));
subplot(1, 2, 2);
plot3(centers(1, :) * 1e6, centers(2, :) * 1e6, centers(3, :) * 1e6, '-');
grid on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
title('Trajectory [\mum]');

out = struct('shp', shp, 'p', p, 'dt', dt, 'centers', centers, 't', t, 'r', r);
outFile = fullfile(fileparts(mfilename('fullpath')), sprintf('out_bd_%s.mat', shp));
save(outFile, '-struct', 'out');
fprintf('  saved %s\n', outFile);
