classdef BeamGauss < Beam
    % BeamGauss < Beam : Paraxial Gaussian beam
    %   A Gaussian beam is defined by its waist w0 and its electric field
    %   amplitudes Ex0 and Ey0 along the x and y axis at its center.
    %
    % BeamGauss properties:
    %   r       -   (read only) radial coordinate [m] < Beam
    %   phi     -   (read only) azimuthal coordinate [m] < Beam
    %   Er      -   radial electric field [V/m] < Beam
    %   Ephi    -   azimuathal electric field [V/m] < Beam
    %   er      -   relative dielectric permettivity < Beam
    %   mr      -   relative magnetic permeability < Beam
    %   lambda0 -   vacuum wavelength [m] < Beam
    %   w0      -   beam waist [m]
    %
    % BeamGauss methods:
    %   BeamGauss   -   constructor
    %   plot        -   plots a beam < Beam
    %   uplus       -   +b (= b) < Beam
    %   uminus      -   -b (inverts fields) < Beam
    %   plus        -   b1+b2 (sums fields) < Beam
    %   minus       -   b1-b2 (subtracts fields) < Beam
    %   mtimes      -   a*b (fields by a scalar) < Beam
    %   intensity   -   intensity [W/m^2] < Beam
    %   power       -   power [W] < Beam
    %   attenuate   -   attenuator < Beam
    %   normalize   -   normalization in power < Beam
    %   expand      -   expander < Beam
    %   iris        -   iris / aperture stop < Beam
    %
    % See also Beam, BeamHG, BeamLG.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Giovanni Volpe
    %   Version: 1.0.0
    %   Date: 2014/01/01

    
    properties
        w0 % beam waist [m]
        Ex0 % x electric field amplitude [V/m]
        Ey0 % y electric field amplitude [V/m]
    end
    methods
        function obj = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,varargin)
            % BEAMGAUSS(Ex0,Ey0,w0,R,Nphi,Nr) constructs a Gaussian beam
            %   with electric field components Ex0 and Ey0 [V/m] and waist
            %   w0 [m]. The radial coordinate goes up to R with Nr
            %   divisions and the azimuthal coordinate has Nphi divisions.
            %
            % BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also BeamGauss, Beam.
            
            Check.isnumeric('Ex0 must be a number',Ex0)
            Check.isnumeric('Ey0 must be a number',Ey0)
            Check.isreal('w0 must be a positive real number',w0,'>',0)

            obj = obj@Beam(R,Nphi,Nr,varargin);
            
            obj.w0 = w0;
            
            Ex = Ex0*exp(-obj.r.^2/w0^2);
            Ey = Ey0*exp(-obj.r.^2/w0^2);
            [obj.Ephi,obj.Er] = Transform.Car2PolVector(obj.phi,Ex,Ey);
        end
    end
end