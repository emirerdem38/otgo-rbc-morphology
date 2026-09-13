function d = ray_cell_distance(t, ln_p1, lnc, cart_size, cart_vec, n_div)

    % Coordinates corresponding to different t values
    line_x = ln_p1(1) + t' .* lnc(1);
    line_y = ln_p1(2) + t' .* lnc(2);
    line_z = ln_p1(3) + t' .* lnc(3);

    % Duplicate vectors to the appropriate size
    line_x = repmat(line_x, 1, cart_size);
    line_y = repmat(line_y, 1, cart_size);
    line_z = repmat(line_z, 1, cart_size);

    % Bring cartesian coordiantes of the cell to the
    % appropriate size
    cell_x = repmat(cart_vec(1,:), n_div, 1);
    cell_y = repmat(cart_vec(2,:), n_div, 1);
    cell_z = repmat(cart_vec(3,:), n_div, 1);

    dx = line_x - cell_x;
    dy = line_y - cell_y;
    dz = line_z - cell_z;

    d = (dx.^2 + dy.^2 + dz.^2).^(1/2);
end