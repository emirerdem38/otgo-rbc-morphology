function eqSym =  LegendreSym(n, m)

% EQSYM - Returns the symbolic expressions for legendre functions of order n
% and degree m.
%

syms g u

switch n
    case 0
        if m == 0
            g = 1;
        end
    
    case 1
        if m == 0
            eq = g;
        elseif m == 1
            eq = -sqrt(1-g^2);
        end

    case 2
        if m == 0
            eq = 1.5 * g^2 - 0.5;
        elseif m == 1
            eq = -3 * g * sqrt(1 - u^2);
        elseif m == 2
            eq = 3 * (1 - g^2);
        end

    case 3
        if m == 0
            eq = 0.5 * (5*g^3 - 3*g);
        elseif m == 1
            eq = 1.5 * (1 - 5*g^2) * sqrt(1 - g^2);
        elseif m == 2
            eq = 15 * g * (1 - g^2);
        elseif m == 3
            eq = -15 * (1 - g^2)^(3/2);
        end

end

eqSym = subs(eq,g,cos(u));

end