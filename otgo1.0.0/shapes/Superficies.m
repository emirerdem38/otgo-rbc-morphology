classdef Superficies < Shape
    % Superficies (Abstract) < Shape : 2D surface embedded in 3D space
    %   A surface is a shape for which is possible to define a normal line
    %   and a tangential plane at each point.
    % 
    % Superficies abstract methods:
    %   intersectionpoint   -   (Abstract) intersection point set with line/vector set
    %   perpline            -   (Abstract) perpendicular line at point
    %   tangentplane        -   (Abstract) tangent plane set passing by point set
    %
    % See also Shape, Point, Vector, SLine, Plane, Spherical, Ellipsoidal, Cylindrical, RBC.
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
        intersectionpoint(sup,d,n)  % intersection point set with line/vector set
        perpline(sup,p)  % perpendicular line at point
        tangentplane(sup,p)  % tangent plane set passing by point set
    end
end