% Geometrical Optics
% Version 1.0.0
%
% Object to define a light ray
%   Ray                 - Ray : Set of rays in 3D
%
% Objects to define optically trapped particles
%   Particle            - (Abstract) < Shape : Optically trappable particle
%   ParticleSpherical   - < <a href="matlab:help Particle">Particle</a> : Spherical optically trappable particle
%   ParticleEllipsoidal - < <a href="matlab:help Particle">Particle</a> : Ellipsoidal optically trappable particle
%   ParticleCylindrical - < <a href="matlab:help Particle">Particle</a> : Cylindrical optically trappable particle
%
% Examples
%   example_ray                 - Example to demonstrate the use of <a href="matlab:help Ray">Ray</a>
%   example_particlespherical   - Example to demonstrate the use of <a href="matlab:help ParticleSpherical">ParticleSpherical</a>
%   example_particleellipsoidal - Example to demonstrate the use of <a href="matlab:help ParticleEllipsoidal">ParticleEllipsoidal</a>
%   example_particlecylindrical - Example to demonstrate the use of <a href="matlab:help ParticleCylindrical">ParticleCylindrical</a>
%
% See also OTS, Transform, shapes, beams.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01


clc
help go