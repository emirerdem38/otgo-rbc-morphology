%% demo_02_static_displacement.m
% Coarse force--displacement sweep along laboratory x (and optionally y).
%
% Matches the paper optical setup. Uses a coarse grid so the demo finishes
% quickly; production paper figures used denser sampling (N ~ 100).
%
% For four-/eight-bumped cells, set theta_eq to the stable tilt from the
% torque--rotation curve (paper: ~0.97 and ~1.29 rad). Use 0 for reference
% orientation (dumbbell, neck, oblate75).

clear; clc; close all;
setup_paths;
P = paper_optical_params();

%% User settings
shp = 'dumbbell';
theta_eq = 0;          % rad; x-rotation of analysis pose (0 = reference)
axis_list = {'x'};     % {'x'}, {'y'}, or {'x','y'}
N = 11;                % demo grid (paper used ~100)
x_range = [-1e-6, 1e-6];
fast_demo = true;      % true -> p=32; false -> paper production p

p = P.p_of(shp);
if fast_demo
    p = min(p, 32);
end
beam = CreateBeam(P.nm, P.f, P.NA, P.Ex0, P.Ey0, P.w0, P.Nphi, P.Nr, P.power);

% Fixed orientation: single lab-x rotation to the analysis pose
particle0 = make_particle(shp, p, Point(0, 0, 0), [], [], [], P);
if abs(theta_eq) > 0
    particle0 = particle0.xrotation(theta_eq);
end
or1 = particle0.cell.or1;
or2 = particle0.cell.or2;
or3 = particle0.cell.or3;

div = linspace(x_range(1), x_range(2), N);

for ia = 1:numel(axis_list)
    ax = lower(axis_list{ia});
    F = zeros(3, N);
    T = zeros(3, N);
    fprintf('Displacement along %s (%d points)...\n', ax, N);
    tic;
    for i = 1:N
        switch ax
            case 'x', c = Point(div(i), 0, 0);
            case 'y', c = Point(0, div(i), 0);
            case 'z', c = Point(0, 0, div(i));
            otherwise, error('Unknown axis %s', ax);
        end
        particle = make_particle(shp, p, c, or1, or2, or3, P);
        [F(:, i), T(:, i)] = sum_force_torque(particle, beam, P.scat);
        if mod(i, 5) == 0 || i == 1 || i == N
            fprintf('  %d/%d\n', i, N);
        end
    end
    fprintf('  done in %.1f s\n', toc);

    % Near-origin linear stiffness estimate (pN/um)
    mid = ceil(N / 2);
    i1 = max(1, mid - 2); i2 = min(N, mid + 2);
    u_um = div * 1e6;
    F_pN = F * 1e12;
    switch ax
        case 'x', Fi = F_pN(1, :);
        case 'y', Fi = F_pN(2, :);
        case 'z', Fi = F_pN(3, :);
    end
    pf = polyfit(u_um(i1:i2), Fi(i1:i2), 1);
    k_est = -pf(1);
    fprintf('  rough k_%s ≈ %.3f pN/um (local linear fit)\n', ax, k_est);

    figure('Color', 'w', 'Name', sprintf('%s F vs %s', shp, ax));
    plot(u_um, F_pN(1, :), '-o', u_um, F_pN(2, :), '-s', u_um, F_pN(3, :), '-^');
    grid on; xlabel(sprintf('%s [\\mum]', ax)); ylabel('Force [pN]');
    legend('F_x', 'F_y', 'F_z', 'Location', 'best');
    title(sprintf('%s: force vs %s-displacement (\\theta_x = %.3f rad)', shp, ax, theta_eq));

    out = struct('shp', shp, 'p', p, 'axis', ax, 'div', div, ...
        'F', F, 'T', T, 'theta_eq', theta_eq, 'k_est_pN_per_um', k_est);
    outFile = fullfile(fileparts(mfilename('fullpath')), ...
        sprintf('out_disp_%s_%s.mat', shp, ax));
    save(outFile, '-struct', 'out');
    fprintf('  saved %s\n', outFile);
end
