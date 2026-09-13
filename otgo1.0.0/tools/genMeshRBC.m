function [X, Y, Z] = genMeshRBC(r, t_min, t_max, d, N)
    % Generates a surface mesh and stores it in obj.cart
    % N controls resolution (similar to 'range' in plot)

    Theta = linspace(0, pi, N+1);
    Phi   = linspace(0, 2*pi, 2*N);
    [Theta, Phi] = meshgrid(Theta, Phi);

    r = r * 1e+6;
    t_min = t_min * 1e+6;
    t_max = t_max * 1e+6;
    d = d * 1e+6;

    % Shape coefficients
    C0 = t_min / (2*r);
    C1 = ((4*r^2)/(2*d^2)) * (-4*C0 + ((t_max*abs(5*d^2-16*r^2))/((4*r^2-d^2)^(3/2))));
    C2 = ((16*r^2)/(2*d^4)) * (2*C0 + ((t_max*(16*r^2-3*d^2))/(sign(5*d^2-16*r^2)*(16*r^2-d^2)^(3/2))));

    % Coordinates
    X = r * sin(Theta) .* cos(Phi);
    Y = r * sin(Theta) .* sin(Phi);
    Z = sqrt((cos(Theta).^2) .* ...
        (C0 + C1*(1-cos(Theta).^2) + C2*(1-cos(Theta).^2).^2)) ...
        .* sign(cos(Theta));
end