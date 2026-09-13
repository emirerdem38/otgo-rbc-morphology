classdef Particle < Shape
    % Particle (Abstract) < Shape : Optically trappable particle
    %   A particle is an object that can be optically trapped, i.e., when
    %   an electromagnetic field impinges on it, scattering occurs and
    %   forces and torques are produced.
    % 
    % Particle abstract methods:
    %   plot            -   (Abstract) plots shape set in 3D < Shape
    %   disp            -   (Abstract) prints shape set < Shape
    %   translate       -   (Abstract) translation < Shape
    %   xrotation       -   (Abstract) rotation around x-axis < Shape
    %   yrotation       -   (Abstract) rotation around y-axis < Shape
    %   zrotation       -   (Abstract) rotation around z-axis < Shape
    %   numel           -   (Abstract) number of shapes < Shape
    %   size            -   (Abstract) size of shape set < Shape
    %   barycenter      -   (Abstract) particle center of mass
    %   scattering      -   (Abstract) scattered rays
    %   force           -   (Abstract) force due to a set of rays
    %   torque          -   (Abstract) torque due to a set of rays
    %
    % See also Shape, ParticleSpherical, ParticleEllipsoidal, ParticleCylindrical, Ray.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Giovanni Volpe
    %   Version: 1.0.0
    %   Date: 2014/01/01
    
    
    methods (Abstract)
        barycenter(par)  % particle center of mass
        scattering(par,r)  % scattered rays
        force(par,r)  % force due to a set of rays
        torque(par,r)  % torque due to a set of rays
    end
end