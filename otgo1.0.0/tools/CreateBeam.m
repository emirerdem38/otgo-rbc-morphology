function r = CreateBeam(nm, f, NA, Ex0, Ey0, w0, Nphi, Nr, power)

L = f*NA/nm; % Iris aperture [m]

bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

r = Ray.beam2focused(bg,f);

for i = 1:size(r,1)
    for j = 1:size(r,2)

        ordX = floor(log10(abs(r.v.X(i,j)))); % Order of Magnitude
        ordY = floor(log10(abs(r.v.Y(i,j))));
        ordZ = floor(log10(abs(r.v.Z(i,j))));
        if ordX < 0 && ordY < 0 && ordX < 0
            A = 10^(max([ordX,ordY,ordZ]));
        else
            A = 1;
        end

        % lncX = r.v.Vx(i,j)*A;
        % lncY = r.v.Vy(i,j)*A;
        % lncZ = r.v.Vz(i,j)*A;

        r.v.Vx(i,j) = r.v.Vx(i,j)*A;
        r.v.Vy(i,j) = r.v.Vy(i,j)*A;
        r.v.Vz(i,j) = r.v.Vz(i,j)*A;

    end
end

end

