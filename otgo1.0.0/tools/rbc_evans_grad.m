function [zx, zy] = rbc_evans_grad(x, y, r, C0, C2, C4)
%RBC_EVANS_GRAD Partial derivatives of Evans–Fung half-thickness z(x,y).
%
%   z = (1/2) * sqrt(1-s) * P(s),  s=(x^2+y^2)/r^2,  P=C0+C2 s+C4 s^2.

    s = (x.^2 + y.^2) ./ r.^2;
    if s >= 1 - 1e-12
        zx = 0;
        zy = 0;
        return;
    end
    P = C0 + C2*s + C4*s^2;
    Pp = C2 + 2*C4*s;
    sqrt1s = sqrt(1 - s);
    % dz/ds = 0.5 * [ -P/(2 sqrt(1-s)) + sqrt(1-s) Pp ]
    dzds = 0.5 * (-P/(2*sqrt1s) + sqrt1s*Pp);
    dsdx = 2*x / r^2;
    dsdy = 2*y / r^2;
    zx = dzds * dsdx;
    zy = dzds * dsdy;
end
