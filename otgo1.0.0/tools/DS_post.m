% Define simulation parameters
total_time = 5;  % seconds

% Compute time vectors for each timestep
dt_1ms = 0.001; t1 = 0:dt_1ms:(size(centers_1ms, 2)-1)*dt_1ms;
dt_2ms = 0.002; t2 = 0:dt_2ms:(size(centers_2ms, 2)-1)*dt_2ms;
dt_3ms = 0.003; t3 = 0:dt_3ms:(size(centers_3ms, 2)-1)*dt_3ms;
dt_5ms = 0.005; t5 = 0:dt_5ms:(size(centers_5ms, 2)-1)*dt_5ms;
dt_10ms = 0.010; t10 = 0:dt_10ms:(size(centers_10ms, 2)-1)*dt_10ms;

% Calculate distances from trap center (assuming trap center at origin)
dist_1ms  = sqrt(sum(centers_1ms.^2, 1));
dist_2ms  = sqrt(sum(centers_2ms.^2, 1));
dist_3ms  = sqrt(sum(centers_3ms.^2, 1));
dist_5ms  = sqrt(sum(centers_5ms.^2, 1));
dist_10ms = sqrt(sum(centers_10ms.^2, 1));

% Plot
figure;
plot(t1,  dist_1ms,  'b-', 'LineWidth', 1.2); hold on;
plot(t2,  dist_2ms,  'c-', 'LineWidth', 1.2);
plot(t3,  dist_3ms,  'm-', 'LineWidth', 1.2);
plot(t5,  dist_5ms,  'r-', 'LineWidth', 1.2);
plot(t10, dist_10ms, 'g-', 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Distance from trap center');
title('Distance vs. Time for Different Time Steps');
legend({'\Delta t = 1 ms', '\Delta t = 2 ms', '\Delta t = 3 ms', ...
    '\Delta t = 5 ms', '\Delta t = 10 ms'}, 'Location', 'best');
grid on;
