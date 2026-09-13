classdef ParticleCELL < Particle
    % ParticleCELL < Particle : Cell optically trappable particle
    %   This object can model a vesicles and cells with different shapes.
    %
    % ParticleCELL properties:
    %   cell  - particle (single CELL)
    %   nm  - medium refractive index
    %   np  - particle refractive index
    %
    % ParticleCELL methods:
    %   ParticleCELL        -   constructor
    %   plot                -   plots particle in 3D
    %   disp                -   prints particle
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of particle (=1)
    %   size                -   size of particle set (=[1 1])
    %   barycenter          -   particle center of mass
    %   scattering          -   scattered rays
    %   force               -   force due to a set of rays
    %   torque              -   torque due to a set of rays
    %
    % See also Particle, Cell, Ray.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Emir Erdem
    %   Version: 1.0.0
    %   Date: 2023/08/23


    properties
        cell    % particle (single CELL)
        nm      % medium refractive index
        np      % particle refractive index
    end
    methods
        function obj = ParticleCELL(p,shape,c,or1,or2,or3,vol,mag,nm,np)
            % PARTICLECELL(p,shape,c,or,mag,vol,nm,np) constructs a CELL particle
            %   with number of discretization points p, shape description shape, 
            %   center c, orientation or, order of magnitude mag,
            %   medium refractive index nm and particle refractive index np.
            %   Note that p must be a positive integer, shape must be a valid shape type, 
            %   c must be a single Point, or must be a single vector, mag must be integer,
            %   nm and np must be real numbers.
            %
            % See also Cell, Point.
            
            Check.isreal('p must be a positive real number',p,'>',0)
            Check.isa('c must be a single Point',c,'Point')
            Check.isa('or1 must be a single Vector',or1,'Vector')
            Check.isa('or2 must be a single Vector',or2,'Vector')
            Check.isa('or3 must be a single Vector',or3,'Vector')
            Check.isnumeric('mag must be a number',mag)
            Check.isnumeric('nm must be a number',nm)
            Check.isnumeric('np must be a number',np)
            
            obj.cell = Cell(p,shape,c,or1,or2,or3,vol,mag);
            obj.nm = nm;
            obj.np = np;
        end
        function h = plot(particle,varargin)
            % PLOT Plots CELL particle in 3D
            %
            % H = PLOT(PARTICLE) plots the CELL particle in 3D. It
            %   returns a graphic handler to the plotted particle.
            %
            % H = PLOT(PARTICLE,'Range',N) sets the numberof discretization points
            %   to be plotted to N. N = cell.p (default) corresponds to the
            %   default number of discretization points when the cell is
            %   constructed.
            %
            % H = PLOT(PARTICLE,'Scale',S) rescales the CELL by S 
            %   before plotting it. S=1 by default. 
            %
            % H = PLOT(PARTICLE,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(PARTICLE,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Cell, surf.

            h = particle.cell.plot(varargin{:});
        end
        function disp(particle)
            % DISP Prints CELL particle
            %
            % DISP(PARTICLE) prints CELL particle PARTICLE.
            %
            % See also Cell.

            disp('<a href="matlab:help ParticleCELL">ParticleCELL</a>');
            disp(particle.cell)
        end
        function particle_t = translate(particle,dp)
            % TRANSLATE 3D translation of CELL particle
            %
            % PARTICLEt = TRANSLATE(PARTICLE,dP) translates CELL particle by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Cell, Vector, Point.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            particle_t = particle;
            particle_t.cell = particle_t.cell.translate(dp);
        end
        function particle_r = xrotation(particle,phi)
            % XROTATION Rotation around x-axis of CELL particle
            %
            % PARTICLEr = XROTATION(PARTICLE,phi) rotates CELL particle 
            %   around x-axis by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.cell = particle_r.cell.xrotation(phi);
        end
        function particle_r = yrotation(particle,phi)
            % XROTATION Rotation around y-axis of CELL particle
            %
            % PARTICLEr = YROTATION(PARTICLE,phi) rotates CELL particle 
            %   around y-axis by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.cell = particle_r.cell.yrotation(phi);
        end
        function particle_r = zrotation(particle,phi)
            % XROTATION Rotation around z-axis of CELL particle
            %
            % PARTICLEr = ZROTATION(PARTICLE,phi) rotates CELL particle 
            %   around z-axis by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.cell = particle_r.cell.zrotation(phi);
        end
        function n = numel(particle)
            % NUMEL Number of particle (=1)
            %
            % N = NUMEL(PARTICLE) number of particles in BEAD (=1).
            %
            % See also Cell.

            n = particle.cell.numel();
        end
        function s = size(particle,varargin)
            % SIZE Size of the particle set (=[1 1])
            % 
            % S = SIZE(PARTICLE) returns a two-element row vector with the number 
            %   of rows and columns in the particle set PARTICLE (=[1 1]).
            %
            % S = SIZE(PARTICLE,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the particle set PARTICLE (=1).
            %
            % See also Cell.

            if ~isempty(varargin)
                s = particle.cell.size(varargin{1});
            else
                s = particle.cell.size();
            end
        end        
        function p = barycenter(particle)
            % BARYCENTER CELL particle center of mass
            %
            % P = BARYCENTER(PARTICLE) returns the point P representing the
            %   center of mass of the CELL particle PARTICLE.
            %
            % See also CELL, Point.
            
            p = particle.cell.c;
        end
        function r_vec = scattering(particle,r,err,N)
            % SCATTERING Scattered rays
            %
            % S = SCATTERING(PARTICLE,R) calculates the set of scattered rays S
            %   due to the scattering of the set of rays R on the CELL
            %   particle PARTICLE.
            %   S is a structure indexed on the scattering events. S(n).r is
            %   the n-th reflected set of rays and S(n).t is the n-th
            %   transmitted set of rays.
            %
            % S = SCATTERING(PARTICLE,R,ERR) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12].
            %
            % S = SCATTERING(PARTICLE,R,ERR,N) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12] or
            %   the number of iterations is N [default N=10].
            %
            % See also Cell, Ray.
            
            if nargin<4
                N = 10;
            end
            
            if nargin<3
                err = 1e-12;
            end
            
            Check.isa('R must be a Ray',r,'Ray')
            Check.isreal('The relative error ERR must be a non-negative real number',err,'>=',0)
            Check.isinteger('The maximum number of itrations N must be a positive integer',N,'>',0)
            
            [r_vec(1).r,r_vec(1).t] = r.snellslaw(particle.cell,particle.nm,particle.np,1);
            
            [r_vec(2).r,r_vec(2).t] = r_vec(1).t.snellslaw(particle.cell,particle.np,particle.nm,2);

            for n = 2:1:N

                [r_vec(n+1).r,r_vec(n+1).t] = r_vec(n).r.snellslaw(particle.cell,particle.np,particle.nm,2);
               
                if r_vec(n+1).r.P < r.P*err | isnan(r_vec(n+1).r.P)
                    break;
                end
            end
            disp('Scattering event completed.')
            
        end
        function f = force(particle,r_vec,r)
            % FORCE Force due to rays
            %
            % F = FORCE(PARTICLE,R_VEC,R) calculates the force due to the scattering 
            %   of the  set of rays R on the CELL particle PARTICLE.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the CELL particle.
            %
            % See also CELL, Ray, Vector.

            cm = PhysConst.c0/particle.nm; % speed of light in medium [m/s]

            fi = (r.P/cm).*r.versor(); % Incoming momentum

            r_r1 = r_vec(1).r;
            fe = (r_r1.P/cm).*r_r1.versor(); % first reflection
            
            for n = 2:1:length(r_vec) % transmissions
                r_t = r_vec(n).t;
                df = (r_t.P/cm).*r_t.versor();
                df.Vx(isnan(r_t.P)) = 0;
                df.Vy(isnan(r_t.P)) = 0;
                df.Vz(isnan(r_t.P)) = 0;
                fe = fe + df;
            end
                       
            f = fi-fe;
            f.X = particle.cell.c.X.*ones(size(f));
            f.Y = particle.cell.c.Y.*ones(size(f));
            f.Z = particle.cell.c.Z.*ones(size(f));
        end
        function T = torque(particle,r_vec,r)
            % TORQUE Torque due to rays
            %
            % T = TORQUE(PARTICLE,R_VEC,R) calculates the torque due to the scattering 
            %   of the  set of rays R_VEC on the CELL particle PARTICLE.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the CELL particle.
            %
            % See also CELL, Ray, Vector.

            cm = PhysConst.c0/particle.nm; % speed of light in medium [m/s]

            C = particle.barycenter(); % Barycenter
            C = Point(C.X*ones(size(r)),C.Y*ones(size(r)),C.Z*ones(size(r)));

            P0 = Point(r_vec(1).r.v.X,r_vec(1).r.v.Y,r_vec(1).r.v.Z); % Application point incoming beam
            CP0 = SLine(C,P0).tovector();
            
            mi = (r.P/cm).*r.versor();
            Ti = CP0*mi; % Incoming angular momentum
            
            r_r1 = r_vec(1).r;
            me = (r_r1.P/cm).*r_r1.versor();
            Te = CP0*me; % first reflection
            
            for n = 2:1:length(r_vec) % transmissions
                r_t = r_vec(n).t;
                Pn = Point(r_vec(n).t.v.X,r_vec(n).t.v.Y,r_vec(n).t.v.Z); % Application point of the n-th transmitted beam
                CPn = SLine(C,Pn).tovector();
                me = (r_t.P/cm).*r_t.versor();
                dT = CPn*me;
                dT.Vx(isnan(r_t.P)) = 0;
                dT.Vy(isnan(r_t.P)) = 0;
                dT.Vz(isnan(r_t.P)) = 0;
                Te = Te + dT;
            end
            
            T = Ti - Te;
            T.X = particle.cell.c.X.*ones(size(T));
            T.Y = particle.cell.c.Y.*ones(size(T));
            T.Z = particle.cell.c.Z.*ones(size(T));
        end
    end
end