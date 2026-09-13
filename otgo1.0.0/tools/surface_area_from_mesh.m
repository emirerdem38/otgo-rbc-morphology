function A = surface_area_from_mesh(x, y, z)

    N = numel(x);
    p = (sqrt(2*N + 1) - 1)/2;
    if abs(p - round(p)) > 1e-10
        error('Point count is not consistent with the Gauss-uniform grid.');
    end
    p = round(p);

    % reshape to grid: (p+1) x (2p)
    X = reshape(x, p+1, 2*p);
    Y = reshape(y, p+1, 2*p);
    Z = reshape(z, p+1, 2*p);

    A = 0;

    % loop over quads in the parametric grid
    for i = 1:p
        for j = 1:2*p
            j2 = j + 1;
            if j == 2*p
                j2 = 1; % periodic wrap in azimuth
            end

            % four corners of one quad
            P1 = [X(i,j),   Y(i,j),   Z(i,j)];
            P2 = [X(i+1,j), Y(i+1,j), Z(i+1,j)];
            P3 = [X(i,j2),  Y(i,j2),  Z(i,j2)];
            P4 = [X(i+1,j2),Y(i+1,j2),Z(i+1,j2)];

            % split quad into two triangles: (P1,P2,P3) and (P2,P4,P3)
            A1 = 0.5 * norm(cross(P2-P1, P3-P1));
            A2 = 0.5 * norm(cross(P4-P2, P3-P2));

            A = A + A1 + A2;
        end
    end
end