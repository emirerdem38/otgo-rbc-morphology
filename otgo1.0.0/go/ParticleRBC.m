classdef ParticleRBC < Particle
    % ParticleRBC < Particle : RBC optically trappable particle
    %   This object can model the bi-concave red blood cells with different 
    %   geometrical properties expressed by shape coefficients that can be optically trapped.
    %
    % ParticleRBC properties:
    %   rbc  - particle (single RBC)
    %   nm   - medium refractive index
    %   np   - particle refractive index
    %
    % ParticleRBC methods:
    %   ParticleRBC         -   constructor
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
    % See also Particle, RBC, Ray.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Emir Erdem
    %   Version: 1.0.0
    %   Date: 2023/20/01

    properties
        rbc     % particle (single RBC)
        nm      % medium refractive index
        np      % particle refractive index
    end
    methods
        function obj = ParticleRBC(c,or1,or2,or3,r,t_min,t_max,d,nm,np,vol)
            % PARTICLERBC(c,or1,or2,or3,r,t_min,t_max,d,nm,np,vol) construct a RBC particle
            %   with center at point c, orientation or, radius r, minimal thickness t_min,
            %   maximal thickness t_max, distance of maximal thickness d,
            %   medium refractive index nm and particle refractive index np.
            %   vol is optional: 'vol' (default) keeps the natural volume, or pass a
            %   positive target volume [µm^3] to scale the geometry uniformly.
            %   Note that c must be a single Point, or must be a single vector,
            %   r, t_min, t_max, d must be positive numbers.
            %
            % See also ParticleRBC, Point, RBC.

            if nargin < 11
                vol = 'vol';
            end

            Check.isa('c must be a single Point',c,'Point')
            Check.isa('or1 must be a single Vector',or1,'Vector')
            Check.isa('or2 must be a single Vector',or2,'Vector')
            Check.isa('or3 must be a single Vector',or3,'Vector')
            Check.isreal('r must be a positive real number',r,'>',0)
            Check.isreal('t_min must be a positive real number',t_min,'>',0)
            Check.isreal('t_max must be a positive real number',t_max,'>',0)
            Check.isreal('d must be a positive real number',d,'>',0)
            Check.isnumeric('nm must be a number',nm)
            Check.isnumeric('np must be a number',np)
            Check.samesize('c,or1,or2,or3,r,t_min,t_max and d must be of size 1',c,or1,or2,or3,r,t_min,t_max,d,1)
            if ischar(vol) || isstring(vol)
                if ~strcmp(vol, 'vol')
                    error('vol must be either ''vol'' or a positive target volume [µm^3].')
                end
            else
                Check.isreal('vol must be real number greater than 0',vol,'>',0)
            end
            
            obj.rbc = RBC(c,or1,or2,or3,r,t_min,t_max,d,vol);
            obj.nm = nm;
            obj.np = np;
           end
        function h = plot(particle,varargin)
            % PLOT Plots RBC particle in 3D
            %
            % H = PLOT(PARTICLE) plots the RBC particle RBC in 3D. It
            %   returns a graphic handler to the plotted particle.
            %
            % H = PLOT(PARTICLE,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(PARTICLE,'Scale',S) rescales the RBC by S 
            %   before plotting it. S=1 by default. 
            %
            % H = PLOT(PARTICLE,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(PARTICLE,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also ParticleRBC, RBC, surf.

            h = particle.rbc.plot(varargin{:});
        end
        function disp(particle)
            % DISP Prints RBC particle
            %
            % DISP(PARTICLE) prints RBC particle PARTICLE.
            %
            % See also ParticleRBC.

            disp(['<a href="matlab:help ParticleRBC">ParticleRBC</a> (r=' num2str(particle.rbc.r) ', vol=' num2str(particle.rbc.vol) ' µm^3, nm=' num2str(particle.nm) ', np=' num2str(particle.np) ') : x=' num2str(particle.rbc.c.X) ' y=' num2str(particle.rbc.c.Y) ' z=' num2str(particle.rbc.c.Z) ' or1.Vx=' num2str(particle.rbc.or1.Vx) ' or1.Vy=' num2str(particle.rbc.or1.Vy) ' or1.Vz=' num2str(particle.rbc.or1.Vz) ' or2.Vx=' num2str(particle.rbc.or2.Vx) ' or2.Vy=' num2str(particle.rbc.or2.Vy) ' or2.Vz=' num2str(particle.rbc.or2.Vz) ' or3.Vx=' num2str(particle.rbc.or3.Vx) ' or3.Vy=' num2str(particle.rbc.or3.Vy) ' or3.Vz=' num2str(particle.rbc.or3.Vz) ' t_min=' num2str(particle.rbc.t_min) ' t_max=' num2str(particle.rbc.t_max) ' d=' num2str(particle.rbc.d)]);
        end
        function particle_t = translate(particle,dp)
            % TRANSLATE 3D translation of RBC particle
            %
            % PARTICLEt = TRANSLATE(PARTICLE,dP) translates RBC particle PARTICLE by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also ParticleRBC, Vector, Point, RBC.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            particle_t = particle;
            particle_t.rbc = particle_t.rbc.translate(dp);
        end
        function particle_r = xrotation(particle,phi)
            % XROTATION Rotation around x-axis of RBC particle
            %
            % PARTICLEr = XROTATION(PARTICLE,phi) rotates RBC particle PARTICLE 
            %   around x-axis by an angle phi [rad].
            %
            % See also ParticleRBC, RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.rbc = particle_r.rbc.xrotation(phi);
        end
        function particle_r = yrotation(particle,phi)
            % XROTATION Rotation around y-axis of RBC particle
            %
            % PARTICLEr = YROTATION(PARTICLE,phi) rotates RBC particle PARTICLE 
            %   around y-axis by an angle phi [rad].
            %
            % See also ParticleRBC, RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.rbc = particle_r.rbc.yrotation(phi);
        end
        function particle_r = zrotation(particle,phi)
            % XROTATION Rotation around z-axis of RBC particle
            %
            % PARTICLEr = ZROTATION(PARTICLE,phi) rotates RBC particle PARTICLE 
            %   around z-axis by an angle phi [rad].
            %
            % See also ParticleRBC, RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            particle_r = particle;
            particle_r.rbc = particle_r.rbc.zrotation(phi);
        end
        function n = numel(particle)
            % NUMEL Number of particle (=1)
            %
            % N = NUMEL(PARTICLE) number of particles in PARTICLE (=1).
            %
            % See also ParticleRBC.

            n = particle.rbc.numel();
        end
        function s = size(particle,varargin)
            % SIZE Size of the particle set (=[1 1])
            % 
            % S = SIZE(PARTICLE) returns a two-element row vector with the number 
            %   of rows and columns in the particle PARTICLE (=[1 1]).
            %
            % S = SIZE(PARTICLE,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the particle PARTICLE (=1).
            %
            % See also ParticleRBC.

            if ~isempty(varargin)
                s = particle.rbc.size(varargin{1});
            else
                s = particle.rbc.size();
            end
        end        
        function p = barycenter(particle)
            % BARYCENTER RBC particle center of mass
            %
            % P = BARYCENTER(PARTICLE) returns the point P representing the
            %   center of mass of the RBC particle PARTICLE.
            %
            % See also ParticleRBC, Point, RBC.
            
            p = particle.rbc.c;
        end
        function r_vec = scattering(particle,r,err,N)
            % SCATTERING Scattered rays
            %
            % S = SCATTERING(PARTICLE,R) calculates the set of scattered rays S
            %   due to the scattering of the set of rays R on the RBC
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
            % See also ParticleRBC, Ray.
            
            if nargin<4
                N = 10;
            end
            
            if nargin<3
                err = 1e-12;
            end
            
            Check.isa('R must be a Ray',r,'Ray')
            Check.isreal('The relative error ERR must be a non-negative real number',err,'>=',0)
            Check.isinteger('The maximum number of itrations N must be a positive integer',N,'>',0)
            
            [r_vec(1).r,r_vec(1).t] = r.snellslaw(particle.rbc,particle.nm,particle.np,1);
           
            [r_vec(2).r,r_vec(2).t] = r_vec(1).t.snellslaw(particle.rbc,particle.np,particle.nm,2);
          
            for n = 2:1:N

                [r_vec(n+1).r,r_vec(n+1).t] = r_vec(n).r.snellslaw(particle.rbc,particle.np,particle.nm,2);
               
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
            %   R_VEC of the set of rays R on the RBC particle PARTICLE.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the RBC particle.
            %
            % See also ParticleRBC, Ray, Vector.

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
            f.X = particle.rbc.c.X.*ones(size(f));
            f.Y = particle.rbc.c.Y.*ones(size(f));
            f.Z = particle.rbc.c.Z.*ones(size(f));
        end
        function T = torque(particle,r_vec,r)
            % TORQUE Torque due to rays
            %
            % T = TORQUE(PARTICLE,R_VEC,R) calculates the torque due to the scattering 
            %   R_VEC of the set of rays R on the RBC particle PARTICLE.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the RBC particle.
            %
            % See also ParticleRBC, Ray, Vector.

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
            T.X = particle.rbc.c.X.*ones(size(T));
            T.Y = particle.rbc.c.Y.*ones(size(T));
            T.Z = particle.rbc.c.Z.*ones(size(T));
        end
    end
end