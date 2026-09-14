function [X, Y, Z] = genMeshRBC(r, t_min, t_max, d, N)
%GENMESHRBC Surface mesh for an Evans–Fung biconcave RBC.
%
%   [X,Y,Z] = GENMESHRBC(r,t_min,t_max,d,N) returns a structured mesh of
%   size (2N) x (N+1). Inputs r,t_min,t_max,d are in metres; outputs X,Y,Z
%   are in micrometres.
%
%   The surface is the Evans–Fung graph
%       z = ± (1/2) Z(rho),
%       Z(rho) = sqrt(1-(rho/r)^2)*(C0+C2*(rho/r)^2+C4*(rho/r)^4),
%   with coefficients from RBC_EVANS_COEFFS so that Z(0)=t_min, Z(d)=t_max
%   and Z'(d)=0.
%
%   See also rbc_evans_coeffs, rbc_volume_area, RBC.

    if nargin < 5 || isempty(N)
        N = 64;
    end

    Theta = linspace(0, pi, N+1);
    Phi   = linspace(0, 2*pi, 2*N);
    [Theta, Phi] = meshgrid(Theta, Phi);

    % Work in micrometres (matches historical OTGO RBC mesh units)
    r = r * 1e6;
    t_min = t_min * 1e6;
    t_max = t_max * 1e6;
    d = d * 1e6;

    [C0, C2, C4] = rbc_evans_coeffs(r, t_min, t_max, d);

    rho = r .* sin(Theta);
    s = (rho ./ r).^2;
    Zfull = sqrt(max(1 - s, 0)) .* (C0 + C2.*s + C4.*s.^2);

    X = rho .* cos(Phi);
    Y = rho .* sin(Phi);
    Z = 0.5 .* Zfull .* sign(cos(Theta));
    % Poles: sign(cos(Theta))=0 at Theta=pi/2 belt is equator (Z=0);
    % at Theta=0/pi, sign is ±1.
    Z(cos(Theta) == 0) = 0;
end
