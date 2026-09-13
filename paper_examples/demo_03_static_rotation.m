%% demo_03_static_rotation.m
% Coarse torque--rotation sweep about laboratory x.
%
% Fresh particle at each angle, then particle.xrotation(theta), matching the
% paper static protocol. Coarse N for a quick demo; paper used denser grids.

clear; clc; close all;
setup_paths;
P = paper_optical_params();

%% User settings
shp = 'fourBump';
N = 11;                       % demo grid (paper used ~100)
theta_range = [-pi / 2, pi / 2];  % use [-pi, pi] for eightBump full axis
fast_demo = true;             % true -> p=32; false -> paper production p

p = P.p_of(shp);
if fast_demo
    p = min(p, 32);
end
beam = CreateBeam(P.nm, P.f, P.NA, P.Ex0, P.Ey0, P.w0, P.Nphi, P.Nr, P.power);
angles = linspace(theta_range(1), theta_range(2), N);

T = zeros(3, N);
fprintf('Rotation about x for %s (%d angles)...\n', shp, N);
tic;
for i = 1:N
    particle = make_particle(shp, p, Point(0, 0, 0), [], [], [], P);
    particle = particle.xrotation(angles(i));
    [~, T(:, i)] = sum_force_torque(particle, beam, P.scat);
    if mod(i, 5) == 0 || i == 1 || i == N
        fprintf('  %d/%d  theta = %.3f rad\n', i, N, angles(i));
    end
end
fprintf('  done in %.1f s\n', toc);

Tx_pNum = T(1, :) * 1e18;   % pN·µm
Ty_pNum = T(2, :) * 1e18;
Tz_pNum = T(3, :) * 1e18;

% Rough kappa near theta = 0 (or nearest sample)
[~, i0] = min(abs(angles));
i1 = max(1, i0 - 2); i2 = min(N, i0 + 2);
pf = polyfit(angles(i1:i2), Tx_pNum(i1:i2), 1);
kappa_est = -pf(1);
fprintf('  rough kappa_x near 0 ≈ %.3f pN·um/rad\n', kappa_est);
fprintf('  (For bumped cells, refit near the stable tilt, not near 0.)\n');

figure('Color', 'w', 'Name', sprintf('%s T vs theta_x', shp));
plot(angles, Tx_pNum, '-o', angles, Ty_pNum, '-s', angles, Tz_pNum, '-^');
grid on; xlabel('\theta_x [rad]'); ylabel('Torque [pN·\mum]');
legend('T_x', 'T_y', 'T_z', 'Location', 'best');
title(sprintf('%s: torque vs rotation about lab x', shp));

out = struct('shp', shp, 'p', p, 'angles', angles, 'T', T, ...
    'kappa_est_near0', kappa_est);
outFile = fullfile(fileparts(mfilename('fullpath')), ...
    sprintf('out_rot_%s_x.mat', shp));
save(outFile, '-struct', 'out');
fprintf('  saved %s\n', outFile);
