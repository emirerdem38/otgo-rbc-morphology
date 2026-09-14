function z = rbc_evans_halfthickness(x, y, r, C0, C2, C4)
%RBC_EVANS_HALFTHICKNESS Half-thickness |z| of the Evans–Fung RBC surface.
%
%   z = RBC_EVANS_HALFTHICKNESS(x,y,r,C0,C2,C4) returns the absolute
%   surface height at (x,y) for the graph
%
%       z = ± (1/2) * sqrt(1-s) * (C0 + C2 s + C4 s^2),  s = (x^2+y^2)/r^2.
%
%   Outside the disk x^2+y^2 > r^2 the result is 0.

    s = (x.^2 + y.^2) ./ r.^2;
    inside = s <= 1;
    z = zeros(size(s), 'like', s);
    if any(inside(:))
        si = s(inside);
        p = C0 + C2.*si + C4.*si.^2;
        z(inside) = 0.5 .* sqrt(1 - si) .* p;
    end
end
