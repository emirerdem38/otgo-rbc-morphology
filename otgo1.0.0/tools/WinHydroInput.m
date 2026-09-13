function WinHydroInput(cell, FileName, plot_flag, varargin)
   % WinHydroInput generates a text file with nodal coordinates and optimized sphere radii.
    % It optionally plots the spheres if plot_flag is set to true.
    %
    % Parameters:
    %   cell - Structure containing nodal Cartesian coordinates (cell.cart.x, cell.cart.y, cell.cart.z)
    %   FileName - Name of the output text file (without directory)
    %   plot_flag - Boolean flag (true to plot spheres, false to skip plotting)
    %   varargin - Optional argument specifying radius limits [min_radius, max_radius]

    % Check for optional radius limits
    if ~isempty(varargin)
        radius_limits = varargin{1}; 
    else
        radius_limits = [0, 1]; % Default: No strict limit
    end
    
    % Save directory
    SaveDir = '/Users/emirerdem/Library/Mobile Documents/com~apple~CloudDocs/Optics Project';
    
    % Construct full file path
    SavePath = fullfile(SaveDir, FileName);

    % Extract nodal coordinates
    x = cell.cart.x;
    y = cell.cart.y;
    z = cell.cart.z;

    % Remove polar-related points (every 1st and 13th from each 13-point row)
    idx = true(length(x), 1);
    for i = 1:floor(length(x)/13)
        base = (i-1)*13;
        idx(base + 1) = false;      % 1st in group
        idx(base + 13) = false;     % 13th in group
    end

    % Apply the mask
    x = x(idx);
    y = y(idx);
    z = z(idx);
    
    % Number of nodes
    num_nodes = length(x);
    
    % Initialize radii (starting with a reasonable guess)
    radius = mean(radius_limits);
    
    % Iteratively adjust the radius to ensure no overlaps
    adjusted_radii = radius * ones(num_nodes, 1); % Initial uniform radius
    
    % Initialize minimum distance to a large number
    min_dist = inf;
    
    % Loop to find the smallest distance between any two points
    for i = 1:num_nodes
        for j = i+1:num_nodes
            dist = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2 + (z(i) - z(j))^2);
            if dist < min_dist
                min_dist = dist;
            end
        end
    end
    
    % Set the radius to half of the minimum distance (just enough to avoid overlap)
    radius = min_dist / 2;

    % Ensure it's within the radius_limits
    radius = max(radius_limits(1), min(radius, radius_limits(2)));
    adjusted_radii(:) = radius;
    
    % Open file for writing
    fileID = fopen(SavePath, 'w');
    
    % Write header
    fprintf(fileID, '1.E-07, !Unit of length for coordinates and radii, cm (10 A)\n');
    fprintf(fileID, '%d, !Number of beads\n', num_nodes);
    
    % Write nodal coordinates and radius with proper alignment
    for i = 1:num_nodes
        fprintf(fileID, '%10.4f %10.4f %10.4f %10.4f\n', x(i)*1e+1, y(i)*1e+1, z(i)*1e+1, adjusted_radii(i)*1e+1);
        %fprintf(fileID, '%10.4f %10.4f %10.4f %10.4f\n', x(i)*1e+1, y(i)*1e+1, z(i)*1e+1, 0.134);
    end
    
    % Close file
    fclose(fileID);
    
    fprintf('File saved successfully at: %s\n', SavePath);

    % Plot spheres if plot_flag is true
    if plot_flag
        figure;
        hold on;
        axis equal;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('3D Visualization of Spheres');

        % Plot spheres
        [sx, sy, sz] = sphere(20); % Generate sphere surface coordinates

        for i = 1:num_nodes
            % Scale and translate the sphere
            surf(radius * sx + x(i), radius * sy + y(i), radius * sz + z(i), ...
                 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Transparent spheres
        end

        % Plot nodal points
        scatter3(x, y, z, 50, 'r', 'filled'); % Red points for node centers

        hold off;
    end
end