function [C0, C2, C4] = rbc_evans_coeffs(r, t_min, t_max, d)
%RBC_EVANS_COEFFS Evans–Fung coefficients for a biconcave RBC.
%
%   [C0,C2,C4] = RBC_EVANS_COEFFS(r,t_min,t_max,d) returns coefficients of
%   the full-thickness profile
%
%       Z(rho) = sqrt(1-(rho/r)^2) * (C0 + C2*(rho/r)^2 + C4*(rho/r)^4)
%
%   with C0 = t_min, Z(d) = t_max and dZ/drho|(rho=d) = 0.
%   All lengths must share the same unit.
%
%   See also genMeshRBC, rbc_volume_area, rbc_evans_halfthickness.

    C0 = t_min;
    s = (d./r).^2;
    s = s(:).';
    C0v = C0(:).';
    tmax = t_max(:).';

    C2 = zeros(size(C0v));
    C4 = zeros(size(C0v));
    for k = 1:numel(C0v)
        sk = s(k);
        fun = @(x) localResidual(x, C0v(k), sk, tmax(k));
        x0 = [7.83, -4.39];
        x = fsolve(fun, x0, optimset('Display', 'off'));
        C2(k) = x(1);
        C4(k) = x(2);
    end

    C2 = reshape(C2, size(C0));
    C4 = reshape(C4, size(C0));
end

function res = localResidual(x, C0, s, t_max)
    C2 = x(1);
    C4 = x(2);
    p = C0 + C2*s + C4*s^2;
    Z = sqrt(1 - s) * p;
    % Z'(rho)=0 iff -p/2 + (1-s) p' = 0, p' = C2 + 2 C4 s
    res = [Z - t_max; -0.5*p + (1 - s)*(C2 + 2*C4*s)];
end
