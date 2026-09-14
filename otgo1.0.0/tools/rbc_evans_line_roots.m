function rts = rbc_evans_line_roots(a1, a2, a3, b1, b2, b3, r, C0, C2, C4)
%RBC_EVANS_LINE_ROOTS Real roots t for ray ∩ Evans–Fung RBC.
%
%   Line: (x,y,z) = (a1,a2,a3) + t (b1,b2,b3).
%   Surface half-thickness z_s = (1/2) sqrt(1-s) (C0 + C2 s + C4 s^2),
%   s = (x^2+y^2)/r^2, so 4 z^2 = (1-s) P(s)^2 with P = C0+C2 s+C4 s^2.

    D1 = (b1^2 + b2^2) / r^2;
    D2 = (2*(a1*b1 + a2*b2)) / r^2;
    D3 = (a1^2 + a2^2) / r^2;
    s_poly = [D1, D2, D3];                         % s(t)

    s2 = conv(s_poly, s_poly);                      % s^2
    P = poly_add(C4*s2, poly_add(C2*[0, s_poly], C0*[0, 0, 1]));
    P2 = conv(P, P);                               % P(s(t))^2
    one_minus_s = poly_add(-s_poly, [0, 0, 1]);    % 1 - s
    RHS = conv(one_minus_s, P2);                   % (1-s) P^2
    LHS = 4 * conv([b3, a3], [b3, a3]);            % 4 (a3+b3 t)^2
    eq = poly_add(RHS, -LHS);

    % Drop leading zeros for roots
    nz = find(abs(eq) > 1e-14*max(1, max(abs(eq))), 1, 'first');
    if isempty(nz)
        rts = [];
        return;
    end
    eq = eq(nz:end);

    try
        rts = roots(eq);
        rts = sort(real(rts(abs(imag(rts)) < 1e-8 & real(rts) > -1e-3)));
    catch
        rts = [];
    end
end

function c = poly_add(a, b)
    a = a(:).';
    b = b(:).';
    n = max(numel(a), numel(b));
    a = [zeros(1, n - numel(a)), a];
    b = [zeros(1, n - numel(b)), b];
    c = a + b;
end
