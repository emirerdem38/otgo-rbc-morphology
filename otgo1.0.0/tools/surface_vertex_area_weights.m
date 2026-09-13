function w = surface_vertex_area_weights(x, y, z, sz)
% SURFACE_VERTEX_AREA_WEIGHTS  Per-vertex patch areas on a structured surface mesh.
%
%   w = surface_vertex_area_weights(x, y, z)
%   w = surface_vertex_area_weights(x, y, z, sz)
%
% x,y,z : surface coordinates (vectors or matrices), same numel.
% sz    : optional [n1 n2] mesh size. If omitted, uses the Cell layout
%         (p+1) x (2p) expected by surface_area_from_mesh.m.
%
% Layout conventions in this project:
%   Cell  (Gauss-Legendre x uniform): sz = [p+1, 2*p], periodic in columns
%   RBC   (genMeshRBC):               sz = [2*p, p+1], periodic in rows
%
% Each quad area is split equally over its four vertices. Returned w sums
% to the total triangulated surface area (not normalized to 1).

    x = x(:); y = y(:); z = z(:);
    N = numel(x);
    if N ~= numel(y) || N ~= numel(z)
        error('x, y, z must have the same number of elements.');
    end

    if nargin < 4 || isempty(sz)
        p = (sqrt(2*N + 1) - 1) / 2;
        if abs(p - round(p)) > 1e-8
            error('N=%d is not a structured OTGO grid with N=(p+1)*2p.', N);
        end
        p = round(p);
        sz = [p+1, 2*p];          % Cell default
        periodDim = 2;
    else
        sz = sz(:).';
        if numel(sz) ~= 2 || prod(sz) ~= N
            error('sz must be [n1 n2] with n1*n2 = numel(x).');
        end
        % Longer index is the periodic azimuthal-like direction
        if sz(2) >= sz(1)
            periodDim = 2;        % Cell-like
        else
            periodDim = 1;        % RBC-like
        end
    end

    X = reshape(x, sz);
    Y = reshape(y, sz);
    Z = reshape(z, sz);
    w = zeros(sz);
    n1 = sz(1); n2 = sz(2);

    if periodDim == 2
        for i = 1:n1-1
            for j = 1:n2
                j2 = j + 1; if j == n2, j2 = 1; end
                P1 = [X(i,j),   Y(i,j),   Z(i,j)];
                P2 = [X(i+1,j), Y(i+1,j), Z(i+1,j)];
                P3 = [X(i,j2),  Y(i,j2),  Z(i,j2)];
                P4 = [X(i+1,j2),Y(i+1,j2),Z(i+1,j2)];
                aquad = triangle_pair_area(P1, P2, P3, P4);
                w(i,j)     = w(i,j)     + aquad/4;
                w(i+1,j)   = w(i+1,j)   + aquad/4;
                w(i,j2)    = w(i,j2)    + aquad/4;
                w(i+1,j2)  = w(i+1,j2)  + aquad/4;
            end
        end
    else
        for i = 1:n1
            i2 = i + 1; if i == n1, i2 = 1; end
            for j = 1:n2-1
                P1 = [X(i,j),   Y(i,j),   Z(i,j)];
                P2 = [X(i2,j),  Y(i2,j),  Z(i2,j)];
                P3 = [X(i,j+1), Y(i,j+1), Z(i,j+1)];
                P4 = [X(i2,j+1),Y(i2,j+1),Z(i2,j+1)];
                aquad = triangle_pair_area(P1, P2, P3, P4);
                w(i,j)     = w(i,j)     + aquad/4;
                w(i2,j)    = w(i2,j)    + aquad/4;
                w(i,j+1)   = w(i,j+1)   + aquad/4;
                w(i2,j+1)  = w(i2,j+1)  + aquad/4;
            end
        end
    end

    w = w(:);
    if any(~isfinite(w)) || sum(w) <= 0
        error('Invalid area weights (non-finite or zero total area).');
    end
end

function a = triangle_pair_area(P1, P2, P3, P4)
    a = 0.5*norm(cross(P2-P1, P3-P1)) + 0.5*norm(cross(P4-P2, P3-P2));
end
