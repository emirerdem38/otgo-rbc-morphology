function particle = make_particle(shp, p, c, or1, or2, or3, P)
%MAKE_PARTICLE Construct a ParticleCELL with paper defaults.
%
%   particle = make_particle(shp, p, c, or1, or2, or3, P)
%   If or* are omitted, body axes coincide with the laboratory axes.

if nargin < 4 || isempty(or1)
    or1 = Vector(0, 0, 0, 1, 0, 0);
    or2 = Vector(0, 0, 0, 0, 1, 0);
    or3 = Vector(0, 0, 0, 0, 0, 1);
end
if nargin < 3 || isempty(c)
    c = Point(0, 0, 0);
end

particle = ParticleCELL(p, shp, c, or1, or2, or3, P.vol, P.mag, P.nm, P.np);
end
