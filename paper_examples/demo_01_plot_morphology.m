%% demo_01_plot_morphology.m
% Plot one parametric morphology used in the paper.
%
% Shape names: 'dumbbell' | 'fourBump' | 'eightBump' | 'neck' | 'oblate75'
% (Paper "oblate" = 'oblate75' in the Cell constructor.)

clear; clc; close all;
setup_paths;
P = paper_optical_params();

shp = 'dumbbell';                 % <-- change morphology here
p = P.p_of(shp);

particle = make_particle(shp, p, Point(0, 0, 0), [], [], [], P);

figure('Color', 'w', 'Name', sprintf('Morphology: %s', shp));
hold on; axis equal; grid on; view(3);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('%s (p = %d)', shp, p));
plot(particle);
particle.cell.c.plot('marker', '.', 'color', 'k');
% Cell.vol / surArea are already in um^3 / um^2 for the parametric mesh
% (same units as Table 1 in the paper).
fprintf('Volume = %.4f um^3, surface area = %.4f um^2\n', ...
    particle.cell.vol, particle.cell.surArea);
