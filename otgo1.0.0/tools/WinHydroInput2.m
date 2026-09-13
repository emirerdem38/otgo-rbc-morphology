function WinHydroInput2(cell, FileName, plot_flag, varargin)
    % WinHydroInput generates a text file with nodal coordinates and optimized sphere radii.
    % It filters out points that are closer than a specified minimum distance to achieve even spacing.
    %
    % Parameters:
    %   cell - Structure with Cartesian coordinates (cell.cart.x, y, z)
    %   FileName - Output text file name
    %   plot_flag - true/false for plotting
    %   varargin - Optional [min_radius, max_radius, min_distance]

    % Parse input arguments
    if length(varargin) >= 1
        radius_limits = varargin{1}; 
    else
        radius_limits = [0, 1];
    end

    if length(varargin) >= 2
        min_dist = varargin{2}; 
    else
        min_dist = 0.5;  % Default: spheres at least 0.5 units apart
    end

    SaveDir = '/Users/emirerdem/Library/Mobile Documents/com~apple~CloudDocs/Optics Project';
    SavePath = fullfile(SaveDir, FileName);

    x_all = cell.cart.x(:);
    y_all = cell.cart.y(:);
    z_all = cell.cart.z(:);

    coords_all = [x_all, y_all, z_all];
    accepted = [];

    % Iteratively accept points that are not too close to already accepted points
    for i = 1:size(coords_all, 1)
        pt = coords_all(i, :);
        if isempty(accepted)
            accepted = pt;
        else
            dists = sqrt(sum((accepted - pt).^2, 2));
            if all(dists > min_dist)
                accepted = [accepted; pt]; %#ok<AGROW>
            end
        end
    end

    % Final coordinates
    x = accepted(:,1);
    y = accepted(:,2);
    z = accepted(:,3);
    num_nodes = length(x);

    % Use average radius within bounds
    radius = mean(radius_limits);
    adjusted_radii = radius * ones(num_nodes, 1);

    % Write file
    fileID = fopen(SavePath, 'w');
    fprintf(fileID, '1.E-07, !Unit of length for coordinates and radii, cm (10 A)\n');
    fprintf(fileID, '%d, !Number of beads\n', num_nodes);
    for i = 1:num_nodes
        %fprintf(fileID, '%10.4f %10.4f %10.4f %10.4f\n', x(i)*1e1, y(i)*1e1, z(i)*1e1, adjusted_radii(i)*1e1);
        fprintf(fileID, '%10.4f %10.4f %10.4f %10.4f\n', x(i)*1e1, y(i)*1e1, z(i)*1e1, 1);
    end
    fclose(fileID);

    fprintf('Filtered file saved at: %s\n', SavePath);

    % Optional plot
    if plot_flag
        figure;
        [X, Y, Z] = sphere(20);
        hold on;
        for i = 1:num_nodes
            surf(X*adjusted_radii(i) + x(i), Y*adjusted_radii(i) + y(i), Z*adjusted_radii(i) + z(i), ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        axis equal;
        title('Evenly Spaced Spheres');
        xlabel('x'); ylabel('y'); zlabel('z');
        view(3); grid on;
        hold off;
    end
end