function geom = compute_shape_anisotropy(x, y, z, w)
% COMPUTE_SHAPE_ANISOTROPY
% Equivalent-ellipsoid axes, aspect ratio, and asphericity from surface points.
%
%   geom = compute_shape_anisotropy(x, y, z)
%   geom = compute_shape_anisotropy(x, y, z, w)
%
% INPUT:
%   x, y, z : surface coordinates (vectors or matrices)
%   w       : optional per-point surface-area weights (same numel as x).
%             If omitted, weights are built from the structured OTGO mesh
%             via surface_vertex_area_weights (Cell layout by default).
%             Pass w = [] to force equal (unweighted) point averaging —
%             not recommended for theta/phi-type grids.
%
% The centroid and covariance are the area-weighted surface measures
%   bar{x} = sum_i w_i x_i / sum_i w_i
%   C      = sum_i w_i (x_i-bar{x})(x_i-bar{x})^T / sum_i w_i
% with a = sqrt(lambda_1) >= b >= c = sqrt(lambda_3), aspect ratio a/c,
% and asphericity
%   A = ((a-b)^2+(b-c)^2+(c-a)^2) / (2(a+b+c)^2).
%
% OUTPUT fields: a,b,c, aspect_ratio, asphericity, eigenvalues, eigenvectors,
%                weights_used ('area'|'uniform'|'custom')

    x = x(:); y = y(:); z = z(:);
    valid = ~(isnan(x) | isnan(y) | isnan(z));
    x = x(valid); y = y(valid); z = z(valid);

    if nargin < 4
        % Default: area-weighted structured mesh (Agnese / ESI clarification)
        w = surface_vertex_area_weights(x, y, z);
        weights_used = 'area';
    elseif isempty(w)
        w = ones(size(x));
        weights_used = 'uniform';
    else
        w = w(:);
        w = w(valid);
        if numel(w) ~= numel(x)
            error('Weight vector must match the number of valid points.');
        end
        weights_used = 'custom';
    end

    if any(w < 0) || ~all(isfinite(w)) || sum(w) <= 0
        error('Weights must be finite, non-negative, and sum to a positive value.');
    end

    w = w / sum(w);   % normalize for centroid / covariance

    xc = sum(w .* x);
    yc = sum(w .* y);
    zc = sum(w .* z);

    X = [x - xc, y - yc, z - zc];
    C = X' * (w .* X);   % sum_i w_i X_i X_i^T

    [V, D] = eig(C);
    eigvals = diag(D);
    [eigvals_sorted, idx] = sort(eigvals, 'descend');
    V = V(:, idx);

    a = sqrt(max(eigvals_sorted(1), 0));
    b = sqrt(max(eigvals_sorted(2), 0));
    c = sqrt(max(eigvals_sorted(3), 0));

    aspect_ratio = a / max(c, eps);
    asphericity = ((a - b)^2 + (b - c)^2 + (c - a)^2) / ...
                  (2 * (a + b + c)^2);

    geom.a = a;
    geom.b = b;
    geom.c = c;
    geom.aspect_ratio = aspect_ratio;
    geom.asphericity = asphericity;
    geom.eigenvalues = eigvals_sorted;
    geom.eigenvectors = V;
    geom.weights_used = weights_used;
    geom.centroid = [xc, yc, zc];
end
