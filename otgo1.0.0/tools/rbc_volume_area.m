function [V, A] = rbc_volume_area(r, t_min, t_max, d, varargin)
%RBC_VOLUME_AREA Enclosed volume and surface area of an Evans–Fung RBC.
%
%   [V,A] = RBC_VOLUME_AREA(r,t_min,t_max,d) integrates the axisymmetric
%   full-thickness profile
%
%       Z(rho) = sqrt(1-(rho/r)^2) * (C0 + C2*(rho/r)^2 + C4*(rho/r)^4)
%
%   with coefficients from RBC_EVANS_COEFFS. Volume and area use the same
%   length unit as the inputs (e.g. µm → µm^3 / µm^2).
%
%   [V,A] = RBC_VOLUME_AREA(...,'N',N) sets the radial quadrature size
%   (default 200001).
%
%   Do not use volume_from_mesh on genMeshRBC output: that helper assumes a
%   Gauss-uniform vesicle grid, not the (2N)x(N+1) RBC mesh.

    N = 200001;
    for k = 1:2:numel(varargin)
        if strcmpi(varargin{k}, 'n')
            N = varargin{k+1};
        end
    end

    [C0, C2, C4] = rbc_evans_coeffs(r, t_min, t_max, d);
    rho = linspace(0, r, N);
    s = (rho ./ r).^2;
    Z = sqrt(max(1 - s, 0)) .* (C0 + C2.*s + C4.*s.^2); % full thickness

    V = trapz(rho, 2*pi*rho.*Z);

    z = Z / 2; % one face (graph)
    dz = gradient(z, rho);
    A = 2 * trapz(rho, 2*pi*rho.*sqrt(1 + dz.^2)); % both faces
end
