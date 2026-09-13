function V = volume_from_mesh(x,y,z)

    N = numel(x);
    p = (sqrt(2*N + 1) - 1)/2;
    if abs(p - round(p)) > 1e-10
        error('Point count is not consistent with the Gauss-uniform grid.');
    end
    p = round(p);

    X = reshape(x, p+1, 2*p);
    Y = reshape(y, p+1, 2*p);
    Z = reshape(z, p+1, 2*p);

    V = 0;

    for i = 1:p
        for j = 1:2*p
            j2 = j + 1;
            if j == 2*p
                j2 = 1;
            end

            P1 = [X(i,j),   Y(i,j),   Z(i,j)];
            P2 = [X(i+1,j), Y(i+1,j), Z(i+1,j)];
            P3 = [X(i,j2),  Y(i,j2),  Z(i,j2)];
            P4 = [X(i+1,j2),Y(i+1,j2),Z(i+1,j2)];

            % triangle 1
            V = V + dot(P1, cross(P2, P3))/6;

            % triangle 2
            V = V + dot(P2, cross(P4, P3))/6;
        end
    end

    V = abs(V);
end