function P = paper_optical_params()
%PAPER_OPTICAL_PARAMS Optical / material parameters used in the paper.
%
% Beam: 5 mW Gaussian, NA = 1.3, ray set 20 x 40, up to 10 scattering events.
% Medium / cell refractive indices: 1.33 / 1.38.

P = struct();
P.nm = 1.33;
P.np = 1.38;
P.f = 5e-6;
P.NA = 1.3;
P.Ex0 = 1e4;
P.Ey0 = 1i * 1e4;
P.w0 = 5e-6;
P.Nphi = 20;   % azimuthal divisions
P.Nr = 40;     % radial divisions
P.power = 5e-3;
P.scat = 10;
P.mag = -6;
P.vol = 'vol';  % use natural enclosed volume of the parametric surface

% Morphology -> recommended surface tessellation density p (paper production values)
P.p_of = containers.Map( ...
    {'dumbbell','fourBump','eightBump','neck','oblate75'}, ...
    {64, 64, 80, 80, 80});
end
