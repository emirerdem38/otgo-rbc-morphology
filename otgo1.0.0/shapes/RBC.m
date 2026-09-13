classdef RBC < Superficies
    % RBC < Superficies : Set of Red Blood Cells (RBC) in 3D
    %   A Red Blood Cell (RBC) is defined by its center c, orientation or, radius r,
    %   minimal thickness t_min, maximal thickness t_max and distance of maximal thickness d.
    %   c must be a Point; or must be a vector; r, t_min, t_max and d must be real 
    %   scalar matrices with the same size.
    %
    % RBC properties:
    %   c     - centers (Point)
    %   or    - orientations (Vector)
    %   r     - radii (matrix)
    %   t_min - minimal thicknesses (matrix)
    %   t_max - maximal thicknesses (matrix)
    %   d     - distances of maximal thicknesses (matrix)
    %   cart    - cartesian coordinates of the RBC
    %   vol     - volume of the RBC
    %   surArea - surface area of the RBC
    %
    % RBC methods:
    %   RBC                 -   constructor
    %   plot                -   plots RBC set in 3D
    %   disp                -   prints RBC set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of RBCs
    %   size                -   size of RBC set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line set at point set
    %   tangentplane        -   tangent plane set passing by point set
    %   generateMesh        -   generates a surface mesh for the RBC
    %
    % See also Shape, Superficies, Point, Vector, SLine, Plane.
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
        c       % centers (Point)
        or1     % reference axis of x direction with respect to the particle reference frame (Vector)
        or2     % reference axis of y direction with respect to the particle reference frame (Vector)
        or3     % reference axis of z direction with respect to the particle reference frame (Vector)
        r       % radii (matrix)
        t_min   % minimal thicknesses (matrix)
        t_max   % maximal thicknesses (matrix)
        d       % distances of maximal thicknesses (matrix)
        cart    % cartesian coordinates of the RBC
        vol     % volume of the RBC
        surArea % surface area of the RBC
    end
    methods
        function obj = RBC(c,or1,or2,or3,r,t_min,t_max,d,vol)
            % RBC(c,or1,or2,or3,r,t_min,t_max,d,vol) constructs a set of RBCs 
            %   with centers c, orientations or, radii r, minimal thicknesses t_min, maximal
            %   thicknesses t_max and distance of maximal thicknesses d.
            %   c must be a Point; or must be a vector; r, t_min, t_max and d must be real 
            %   scalar matrices with the same size.
            %   vol must be either the string 'vol' (default geometry volume) or a
            %   positive real target volume [µm^3]; in the latter case all linear
            %   dimensions are scaled uniformly to match the requested volume.
            %
            % See also RBC, Point, Vector.

            if nargin < 9
                vol = 'vol';
            end

            Check.isa('c must be a Point',c,'Point')
            Check.isa('or1 must be a Vector',or1,'Vector')
            Check.isa('or2 must be a Vector',or2,'Vector')
            Check.isa('or3 must be a Vector',or3,'Vector')
            Check.isreal('r must be real matrix greater than 0',r,'>',0)
            Check.isreal('t_min must be real matrix greater than 0',t_min,'>',0)
            Check.isreal('t_max must be real matrix greater than 0',t_max,'>',0)
            Check.isreal('d must be real matrix greater than 0',d,'>',0)
            Check.samesize('c,r,t_min,t_max and d must have the same size',c,r,t_min,t_max,d)

            [X, Y, Z] = genMeshRBC(r, t_min, t_max, d, 64);
            V = volume_from_mesh(X(:), Y(:), Z(:));

            if ischar(vol) || isstring(vol)
                if ~strcmp(vol, 'vol')
                    error('vol must be either ''vol'' or a positive real number.')
                end
                k = 1;
                obj.vol = V;
            else
                Check.isreal('vol must be real number greater than 0',vol,'>',0)
                k = (vol/V)^(1/3);
                obj.vol = vol;
                r = r * k;
                t_min = t_min * k;
                t_max = t_max * k;
                d = d * k;
                [X, Y, Z] = genMeshRBC(r, t_min, t_max, d, 64);
            end

            obj.c = c;
            obj.or1 = Vector(0,0,0,or1.Vx,or1.Vy,or1.Vz).normalize;
            obj.or2 = Vector(0,0,0,or2.Vx,or2.Vy,or2.Vz).normalize;
            obj.or3 = Vector(0,0,0,or3.Vx,or3.Vy,or3.Vz).normalize;
            obj.r = r;
            obj.t_min = t_min;
            obj.t_max = t_max;
            obj.d = d;
            obj.cart = d3Vec([X(:), Y(:), Z(:)]); 
            obj.surArea = surface_area_from_mesh(X(:), Y(:), Z(:));
        end
        function h = plot(rbc,varargin)
            % PLOT Plots RBC set in 3D
            %
            % H = PLOT(RBC) plots the set of RBCs RBC in 3D. It returns a
            %   graphic handler to the plotted set of RBCs.
            %
            % H = PLOT(RBC,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(RBC,'Scale',S) rescales the coordinates of
            %   by S before plotting the set of RBCs. S=1 by default. 
            %
            % H = PLOT(RBC,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(RBC,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also RBC, surf.

            % Range to be plotted
            N = 24;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    N = varargin{n+1};
                end
            end
            Theta = 0:pi/N:pi;
            Phi = -2*pi:pi/N:2*pi;
            
            % RBC model does not work properly in micro scale, so it
            % should first be scaled up to unit size, then brought back to
            % its original sizes. This is same for intersectionpoint and
            % perpline methods.
            ord = floor(log10(abs(rbc.r))); % Order of Magnitude of RBC size

            % Scale RBC set to unit sizes
            rbc.r = rbc.r .* 10.^(ord*(-1));
            rbc.t_min = rbc.t_min .* 10.^(ord*(-1));
            rbc.t_max = rbc.t_max .* 10.^(ord*(-1));
            rbc.d = rbc.d .* 10.^(ord*(-1));
            rbc.c = rbc.c .* 10.^(ord*(-1));
        
            % Scale
            S = 10.^ord;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = S*varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                end
            end
            
            % Color level
            C = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'colorlevel')
                    C = varargin{n+1};
                    Check.isreal('The scaling factor must be a real number',C)
                end
            end

            % Plots

            % Shape coefficients from 'Theoretical investigation of erythrocytes optical trapping in ray optics approximation'
            C0 = rbc.t_min ./ (2*rbc.r);
            C1 = ((4*rbc.r.^2)./(2*rbc.d.^2)) .* (-4*C0 + ((rbc.t_max.*abs(5*rbc.d.^2-16.*rbc.r.^2))./((4*rbc.r.^2-rbc.d.^2).^(3/2))));
            C2 = ((16*rbc.r.^2)./(2*rbc.d.^4)) .* (2*C0 + ((rbc.t_max.*(16*rbc.r.^2-3*rbc.d.^2))./(sign(5*rbc.d.^2-16*rbc.r.^2).*(16*rbc.r.^2-rbc.d.^2).^(3/2))));
            
            % Different shape coefficients can also be used, below shape coefficients 
            % corresponds to experimental findings of Evan Evans and Yuan-Cheng Fung on
            % 'Improved Measurements of the Erythrocyte Geometry'. To change the shape 
            % coefficients that OTGO operates, it should be changed on all relevant methods: 
            % plot, intersectionpoint, perpline.
            
            % C0 = 0.81*ones(size(rbc));
            % C1 = 7.83*ones(size(rbc));
            % C2 = -4.39*ones(size(rbc));

            ht = zeros(rbc.size());
            for m = 1:1:rbc.size(1)
                for n = 1:1:rbc.size(2)

                    % Mathematical expression of RBC describes a RBC centered at the origin 
                    % with a default orientation; Hence, after the nodal points are calculated, they
                    % must be translated and rotated with regard to the configuration of the RBC.

                    % Points to be plotted
                    X = rbc.r(m,n)*cos(Theta')*cos(Phi);
                    Y = rbc.r(m,n)*sin(Theta')*cos(Phi);
                    Z = ones(size(Theta')) * (((1-cos(Phi).^2) .* (C0(m,n) + C1(m,n)*cos(Phi).^2 + C2(m,n)*cos(Phi).^4)).^(1/2).* sign(Phi));
                    
                    R = [rbc.or1.Vx rbc.or2.Vx rbc.or3.Vx;
                         rbc.or1.Vy rbc.or2.Vy rbc.or3.Vy;
                         rbc.or1.Vz rbc.or2.Vz rbc.or3.Vz];
            
                    % Calculate R using the following way if the cell will have a
                    % single reference vector (along +z) denoting its orientation.

                    % angles = orientation(RBC(Point(rbc.c.X(m,n),rbc.c.Y(m,n),rbc.c.Z(m,n)),Vector(0,0,0,rbc.or.Vx(m,n),rbc.or.Vy(m,n),rbc.or.Vz(m,n)),rbc.r(m,n),rbc.t_min(m,n),rbc.t_max(m,n),rbc.d(m,n)),false);
                    % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
                    % 
                    % %Set of rotation matrices
                    % Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
                    % Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
                    % Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];
                    % 
                    % R = Rz * Ry * Rx;
                    
                    X_rot = zeros(size(X,1),size(X,2));
                    Y_rot = zeros(size(Y,1),size(Y,2));
                    Z_rot = zeros(size(Z,1),size(Z,2));
                    
                    % Rotation
                    for i = 1:size(X,1)
                        for j = 1:size(X,2)
                            vertex = [X(i,j); Y(i,j); Z(i,j)];
                            rot_vertex = R * vertex;
                            X_rot(i,j) = rot_vertex(1);
                            Y_rot(i,j) = rot_vertex(2);
                            Z_rot(i,j) = rot_vertex(3);
                        end
                    end

                    % Translation
                    X = rbc.c.X(m,n) + X_rot;
                    Y = rbc.c.Y(m,n) + Y_rot;
                    Z = rbc.c.Z(m,n) + Z_rot;
                    
                    % FaceAlpha slightly below Cell: thin biconcave shell shows two
                    % semi-transparent faces, which otherwise reads more saturated.
                    
                    % ht(m,n) = mesh(S(m,n)*X,S(m,n)*Y,S(m,n)*Z,C*ones(size(X)), ...
                    %     'EdgeColor','k', 'EdgeAlpha', 0.25, 'FaceColor', 'r', 'FaceAlpha', 0.25, ...
                    %     'FaceLighting', 'none', 'EdgeLighting', 'none');

                    % red - classic
                    ht(m,n) = mesh(S(m,n)*X,S(m,n)*Y,S(m,n)*Z,C*ones(size(X)), ...
                        'EdgeColor','k', 'EdgeAlpha', 1, 'FaceColor', 'r', 'FaceAlpha', 1, ...
                        'FaceLighting', 'none', 'EdgeLighting', 'none');
                end
                
                % draw the center (optional; suppress with 'ShowCenter', false)
                showCenter = false;
                for narg = 1:2:length(varargin)
                    if strcmpi(varargin{narg}, 'showcenter')
                        showCenter = logical(varargin{narg+1});
                    end
                end
                if showCenter
                    set(gca,'SortMethod','childorder')
                    plot3(S.*rbc.c.X, S.*rbc.c.Y, S.*rbc.c.Z, 'o', ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceColor', [1 1 0], ...
                        'MarkerSize', 10, ...
                        'LineWidth', 1.2);
                end

                % --- draw z-axis (beam direction) ---
                % L = max([range(X(:)*S), range(Y(:)*S), range(Z(:)*S)]);
                % 
                % quiver3(0, 0, 0, 0, 0, 2*L, ...
                %     'k', ...
                %     'LineWidth', 3, ...
                %     'MaxHeadSize', 0.2);
                % 
                % % Put label almost ON the arrow tip
                % text(0, 0, 2.3*L, '+z', ...
                %     'FontSize', 14, ...
                %     'FontWeight', 'bold', ...
                %     'HorizontalAlignment', 'center', ...
                %     'VerticalAlignment', 'top');

                % xlabel('x')
                % ylabel('y')
                % zlabel('z')
                % xlim([-4e-6 6e-6])
                % ylim([-5e-6 5e-6])
                % zlim([-2e-5 2e-5])
                % --- Axes style (Fig. 1 / publication) ---
                grid on
                axis equal
                axis on
                xlabel('x')
                ylabel('y')
                zlabel('z')
                %axis off
                view(3)
                % axis on   % enable for interactive exploration
                
                % Hide ticks and labels but keep grid
                % ax = gca;
                % ax.XColor = 'none';
                % ax.YColor = 'none';
                % ax.ZColor = 'none';
                % ax.GridColor = [0 0 0];   % or whatever color you want
                % ax.GridAlpha = 0.3;
                % grid(ax, 'on')
            end
            
            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'range') && ~strcmpi(varargin{n},'scale') && ~strcmpi(varargin{n},'colorlevel') && ~strcmpi(varargin{n},'showcenter')
                    set(ht,varargin{n},varargin{n+1});
                end
            end
            
            % Output if needed
            if nargout>0
                h = ht;
            end
        end
        function disp(rbc)
            % DISP Prints RBC set
            %
            % DISP(RBC) prints set of RBCs RBC.
            %
            % See also RBC.
            
            disp(['<a href="matlab:help RBC">RBC</a> [' int2str(rbc.size) '] : X Y Z OR1.VX OR1.VY OR1.VZ OR2.VX OR2.VY OR2.VZ OR3.VX OR3.VY OR3.VZ R T_MIN T_MAX D']);
            disp([reshape(rbc.c.X,1,rbc.numel());reshape(rbc.c.Y,1,rbc.numel());reshape(rbc.c.Z,1,rbc.numel());...
                reshape(rbc.or1.Vx,1,rbc.numel());reshape(rbc.or1.Vy,1,rbc.numel());reshape(rbc.or1.Vz,1,rbc.numel());...
                reshape(rbc.or2.Vx,1,rbc.numel());reshape(rbc.or2.Vy,1,rbc.numel());reshape(rbc.or2.Vz,1,rbc.numel());...
                reshape(rbc.or3.Vx,1,rbc.numel());reshape(rbc.or3.Vy,1,rbc.numel());reshape(rbc.or3.Vz,1,rbc.numel());...
                reshape(rbc.r,1,rbc.numel());reshape(rbc.t_min,1,rbc.numel());reshape(rbc.t_max,1,rbc.numel());reshape(rbc.d,1,rbc.numel())]);
        end
        function rbc_t = translate(rbc,dp)
            % TRANSLATE 3D translation of RBC set
            %
            % RBCt = TRANSLATE(RBC,dP) translates set of RBCs RBC by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also RBC, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            rbc_t = rbc;
            rbc_t.c = rbc.c.translate(dp);
        end
        function rbc_r = xrotation(rbc,phi)
            % XROTATION Rotation around x-axis of RBC set
            %
            % RBCr = XROTATION(RBC,phi) rotates set of RBCs RBC around x-axis 
            %   by an angle phi [rad].
            %
            % See also RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            rbc_r = rbc;
            rbc_r.c = rbc.c.xrotation(phi);
            rbc_r.or1 = rbc.or1.xrotation(phi);
            rbc_r.or2 = rbc.or2.xrotation(phi);
            rbc_r.or3 = rbc.or3.xrotation(phi);
        end
        function rbc_r = yrotation(rbc,phi)
            % YROTATION Rotation around y-axis of RBC set
            %
            % RBCr = YROTATION(RBC,phi) rotates set of RBCs RBC around y-axis 
            %   by an angle phi [rad].
            %
            % See also RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            rbc_r = rbc;
            rbc_r.c = rbc.c.yrotation(phi);
            rbc_r.or1 = rbc.or1.yrotation(phi);
            rbc_r.or2 = rbc.or2.yrotation(phi);
            rbc_r.or3 = rbc.or3.yrotation(phi);
        end
        function rbc_r = zrotation(rbc,phi)
            % ZROTATION Rotation around z-axis of RBC set
            %
            % RBCr = ZROTATION(RBC,phi) rotates set of RBCs RBC around z-axis 
            %   by an angle phi [rad].
            %
            % See also RBC.

            Check.isreal('The rotation angle phi must be a real number',phi)

            rbc_r = rbc;
            rbc_r.c = rbc.c.zrotation(phi);
            rbc_r.or1 = rbc.or1.zrotation(phi);
            rbc_r.or2 = rbc.or2.zrotation(phi);
            rbc_r.or3 = rbc.or3.zrotation(phi);
        end
        function n = numel(rbc)
            % NUMEL Number of RBCs
            %
            % N = NUMEL(RBC) number of RBCs in set RBC.
            %
            % See also RBC.

            n = numel(rbc.c);
        end
        function s = size(rbc,varargin)
            % SIZE Size of the RBC set
            % 
            % S = SIZE(RBC) returns a two-element row vector with the number 
            %   of rows and columns in the RBC set RBC.
            %
            % S = SIZE(RBC,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the RBC set RBC.
            %
            % See also RBC.

            if ~isempty(varargin)
                s = rbc.c.size(varargin{1});
            else
                s = rbc.c.size();
            end
        end     
        function p = intersectionpoint(rbc,d,n)
            % INTERSECTIONPOINT Intersection point between RBC and line/vector/ray set
            %
            % P = INTERSECTIONPOINT(RBC,D,N) calculates intersection points 
            %   between a set of lines (or vectors/rays) D and a RBC. There can
            %   be up to 6 intersection points depending on the shape coefficients, 
            %   varargin can take values {1,2,3,4,5,6}
            %   
            %   If D does not intersect RBC, the coordinates of P are NaN.
            % 
            % See also RBC, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            Check.isinteger('A must be either 1,2,3,4,5 or 6',n,'>=',1,'<=',6)

            if isa(d,'SLine')
                ln_set = d;
            else
                ln_set = d.toline();                
            end
            
            % Order of Magnitude of RBC size
            ord = floor(log10(abs(rbc.r))); 
            
            % Scale RBC to unit sizes
            rbc.r = rbc.r * 10^(ord*(-1));
            rbc.t_min = rbc.t_min * 10^(ord*(-1));
            rbc.t_max = rbc.t_max * 10^(ord*(-1));
            rbc.d = rbc.d * 10^(ord*(-1));
            rbc.c = rbc.c * 10^(ord*(-1));
            
            S = 10^ord; %Scaling factor
            
            % % Calculate the corresponding orientation and translation of RBC from the origin. 
            % % This will then be applied to the line to map it to a RBC
            % % which is centered at the origin.
            % angles = orientation(rbc,true);
            % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
            
            % Preset matrices for calculated intersection points
            pX = NaN(size(ln_set,1),size(ln_set,2));
            pY = NaN(size(ln_set,1),size(ln_set,2));
            pZ = NaN(size(ln_set,1),size(ln_set,2));

            % Matrices for line set
            p1X = NaN(size(ln_set,1),size(ln_set,2));
            p1Y = NaN(size(ln_set,1),size(ln_set,2));
            p1Z = NaN(size(ln_set,1),size(ln_set,2));
            p2X = NaN(size(ln_set,1),size(ln_set,2));
            p2Y = NaN(size(ln_set,1),size(ln_set,2));
            p2Z = NaN(size(ln_set,1),size(ln_set,2));
            
            % Shape coefficients
            C0 = rbc.t_min / (2*rbc.r);
            C1 = ((4*rbc.r^2)/(2*rbc.d^2)) * (-4*C0 + ((rbc.t_max*abs(5*rbc.d^2-16*rbc.r^2))/((4*rbc.r^2-rbc.d^2)^(3/2))));
            C2 = ((16*rbc.r^2)/(2*rbc.d^4)) * (2*C0 + ((rbc.t_max*(16*rbc.r^2-3*rbc.d^2))/(sign(5*rbc.d^2-16*rbc.r^2)*(16*rbc.r^2-rbc.d^2)^(3/2))));
            
            %C0 = 0.81;
            %C1 = 7.83;
            %C2 = -4.39;
            
            % RBC model
            %rbc_eq = @(x, y) abs((1 - (x^2 + y^2)/rbc.r^2) * (C0 + C1*((x^2 + y^2)/rbc.r^2) + C2*((x^2 + y^2)^2/rbc.r^4)))^(1/2);

            rbc_eq = @(x, y) abs((1 - (x.^2 + y.^2)/rbc.r^2) .* (C0 + C1*((x.^2 + y.^2)/rbc.r^2) + C2*((x.^2 + y.^2).^2/rbc.r^4))).^(1/2);
            
            % Scale up the lines to unit size
            ordx = floor(log10(abs(ln_set.p1.X)));
            ordy = floor(log10(abs(ln_set.p1.Y)));
            ordz = floor(log10(abs(ln_set.p1.Z)));
            if max([ordx,ordy,ordz]) <= -6
                ln_set.p1 = ln_set.p1 * 10^6;
                ln_set.p2 = ln_set.p2 * 10^6;
            end
            
            R = [rbc.or1.Vx rbc.or2.Vx rbc.or3.Vx;
                 rbc.or1.Vy rbc.or2.Vy rbc.or3.Vy;
                 rbc.or1.Vz rbc.or2.Vz rbc.or3.Vz];

            % Adjust the line so that it intersects with the cell
            % that is centered at the origin and with a default
            % orientation (Because this is what the mathematical expressions denote)
            ln_set = ln_set.translate(-1*rbc.c);

            for i = 1:numel(ln_set)
                p1 = R' * [ln_set.p1.X(i); ln_set.p1.Y(i); ln_set.p1.Z(i)];
                p2 = R' * [ln_set.p2.X(i); ln_set.p2.Y(i); ln_set.p2.Z(i)];
                p1X(i) = p1(1); p1Y(i) = p1(2); p1Z(i) = p1(3);
                p2X(i) = p2(1); p2Y(i) = p2(2); p2Z(i) = p2(3);
            end

            ln_set = SLine(Point(p1X,p1Y,p1Z),Point(p2X,p2Y,p2Z));
            lnc = normalize(ln_set.p2 - ln_set.p1);

            % Use the below code if RBC orientation is expressed by a
            % single vector (typically along +z)

            % % Adjust the line set so that it intersects with the RBC
            % % that is centered at the origin and with a default
            % % orientation (because this is what the mathematical expression denotes).
            % 
            % ln_set = ln_set.translate(-1*rbc.c);
            % ln_set = ln_set.xrotation(x_angle);
            % ln_set = ln_set.yrotation(y_angle);
            % ln_set = ln_set.zrotation(z_angle);
            % 
            % lnc = normalize(ln_set.p2 - ln_set.p1);
            
            % Coefficients to simplify the analytical expression for the intersection points
            a1 = ln_set.p1.X; a2 = ln_set.p1.Y; a3 = ln_set.p1.Z;
            b1 = lnc.X; b2 = lnc.Y; b3 = lnc.Z;
            D1 = (b1.^2+b2.^2) / (rbc.r^2);
            D2 = (2*(a1.*b1+a2.*b2)) / (rbc.r^2);
            D3 = ((a1.^2)+(a2.^2)) / (rbc.r^2);

            E1 = C2*D1.^2;
            E2 = 2*C2*D1.*D2;
            E3 = 2*C2*D1.*D3 + C2*D2.^2 + C1*D1;
            E4 = 2*C2*D2.*D3 + C1*D2;
            E5 = C0 + C1*D3 + C2*D3.^2;

            F1 = D1.*E1;
            F2 = D1.*E2 + D2.*E1;
            F3 = D1.*E3 + D2.*E2 + D3.*E1 - E1;
            F4 = D1.*E4 + D2.*E3 + D3.*E2 - E2;
            F5 = D1.*E5 + D2.*E4 + D3.*E3 - E3 + b3.^2;
            F6 = D2.*E5 + D3.*E4 - E4 + 2*a3.*b3;
            F7 = D3.*E5 - E5 + a3.^2;

            % Iterate for all lines
            parfor i = 1:numel(ln_set)
            
                ln = SLine(Point(ln_set.p1.X(i),ln_set.p1.Y(i),ln_set.p1.Z(i)),Point(ln_set.p2.X(i),ln_set.p2.Y(i),ln_set.p2.Z(i)));
                lnc_new = Point(lnc.X(i),lnc.Y(i),lnc.Z(i));

                % Default polynomial solver of MATLAB
                eq = [F1(i) F2(i) F3(i) F4(i) F5(i) F6(i) F7(i)];
                try
                    rts = roots(eq);
                    rts = sort(rts(imag(rts) == 0 & real(rts) > -1e-3));
                catch
                    rts = [NaN; NaN; NaN; NaN; NaN; NaN];
                end

                % Use the below solver of MATLAB Symbollic toolbox for
                % more accuracy (computational time increases drastically). 
                
                % syms t
                % eqn_intersection = F1(i,j)*t^6 + F2(i,j)*t^5 + F3(i,j)*t^4 + F4(i,j)*t^3 + F5(i,j)*t^2 + F6(i,j)*t + F7(i,j) == 0;
                % t_sol = double(solve(eqn_intersection,t));
                % t_real_sol = double(t_sol(imag(t_sol) == 0));
                % rts = double(t_real_sol(real(t_real_sol) > -1e-3));

                pts = [];

                for k = 1:size(rts,1)
                    point_x = ln.p1.X+rts(k)*lnc_new.X; 
                    point_y = ln.p1.Y+rts(k)*lnc_new.Y;
                    point_z = ln.p1.Z+rts(k)*lnc_new.Z;

                    if (point_x^2+point_y^2) < rbc.r^2 + 1e-6
                        if point_z >= 0
                            point_z = rbc_eq(point_x,point_y);
                        else
                            point_z = -rbc_eq(point_x,point_y);          
                        end
                        pts = [pts,[point_x; point_y; point_z]];
                    end
                end

                if size(pts,2) >= n
                    pX(i) = pts(1,n);
                    pY(i) = pts(2,n);
                    pZ(i) = pts(3,n);
                end
            end

            int_point = Point(pX,pY,pZ);

            % Reverse-rotate and Reverse-translate the Intersection Point
            for i = 1:numel(int_point)
                vertex = [int_point.X(i); int_point.Y(i); int_point.Z(i)];
                rot_vertex = R * vertex;
                pX(i) = rot_vertex(1);
                pY(i) = rot_vertex(2);
                pZ(i) = rot_vertex(3);
            end

            int_point = Point(pX,pY,pZ).translate(rbc.c);
            
            % Scale back to original size
            int_point.X = int_point.X * S;
            int_point.Y = int_point.Y * S;
            int_point.Z = int_point.Z * S;
            
            p = int_point; % desired intersection poit
        end
        function ln = perpline(rbc,p_set)
            % PERPLINE Line perpendicular to RBC passing by point set
            % 
            % LN = PERPLINE(RBC,P-SET) calculates the line set LN perpendicular 
            %   to the RBC set RBC and passing by the point set P.
            % 
            % See also RBC, Point, SLine.

            Check.isa('P must be a Point',p_set,'Point')
        
            % Order of Magnitude of RBC size
            ord = floor(log10(abs(rbc.r))); 

            % Scale RBC to unit sizes
            rbc.r = rbc.r * 10^(ord*(-1));
            rbc.t_min = rbc.t_min * 10^(ord*(-1));
            rbc.t_max = rbc.t_max * 10^(ord*(-1));
            rbc.d = rbc.d * 10^(ord*(-1));
            rbc.c = rbc.c * 10^(ord*(-1));
            
            % Scale point set to unit sizes
            p_set = p_set * (10^-ord);

            R = [rbc.or1.Vx rbc.or2.Vx rbc.or3.Vx;
                 rbc.or1.Vy rbc.or2.Vy rbc.or3.Vy;
                 rbc.or1.Vz rbc.or2.Vz rbc.or3.Vz];

            % Use the below code if the rbc orientation is expres using a
            % single axis (typically along +z)

            % angles = orientation(rbc,true);
            % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
            % 
            % % Adjust the point set so that it intersects with the RBC
            % % that is centered at the origin and with a default
            % % orientation (because this is what the mathematical expression denotes).
            % p_set = p_set.translate(-1*rbc.c);
            % p_set = p_set.xrotation(x_angle);
            % p_set = p_set.yrotation(y_angle);
            % p_set = p_set.zrotation(z_angle);

            p_set = p_set.translate(-1*rbc.c);

            % Matrices for calculated point set
            pX = NaN(size(p_set,1),size(p_set,2));
            pY = NaN(size(p_set,1),size(p_set,2));
            pZ = NaN(size(p_set,1),size(p_set,2));

            for i = 1:numel(p_set)
                vertex = [p_set.X(i); p_set.Y(i); p_set.Z(i)];
                rot_vertex = R' * vertex;
                pX(i) = rot_vertex(1);
                pY(i) = rot_vertex(2);
                pZ(i) = rot_vertex(3);
            end

            p_set = Point(pX,pY,pZ);
            
            % Shape Coefficients
            C0 = rbc.t_min / (2*rbc.r);
            C1 = ((4*rbc.r^2)/(2*rbc.d^2)) * (-4*C0 + ((rbc.t_max*abs(5*rbc.d^2-16*rbc.r^2))/((4*rbc.r^2-rbc.d^2)^(3/2))));
            C2 = ((16*rbc.r^2)/(2*rbc.d^4)) * (2*C0 + ((rbc.t_max*(16*rbc.r^2-3*rbc.d^2))/(sign(5*rbc.d^2-16*rbc.r^2)*(16*rbc.r^2-rbc.d^2)^(3/2))));

            % C0 = 0.81;
            % C1 = 7.83;
            % C2 = -4.39;
            
            % Partial derivatives of RBC model (necessary for gradient calculations)
            rbc_x = @(x, y) (1/2)*(abs((1-((x^2+y^2)/rbc.r^2))*(C0+C1*((x^2+y^2)/rbc.r^2)+C2*((x^2+y^2)^2/rbc.r^4)))^(-1/2)) * ( ((-2*x)/rbc.r^2)*(C0+C1*((x^2+y^2)/rbc.r^2)+C2*((x^2+y^2)^2/rbc.r^4)) + (((2*C1*x)/rbc.r^2) + ((4*C2*x^3)/rbc.r^4) + ((4*C2*y^2*x)/rbc.r^4) )*(1-((x^2+y^2)/(rbc.r^2))));
            rbc_y = @(x, y) (1/2)*(abs((1-((x^2+y^2)/rbc.r^2))*(C0+C1*((x^2+y^2)/rbc.r^2)+C2*((x^2+y^2)^2/rbc.r^4)))^(-1/2)) * ( ((-2*y)/rbc.r^2)*(C0+C1*((x^2+y^2)/rbc.r^2)+C2*((x^2+y^2)^2/rbc.r^4)) + (((2*C1*y)/rbc.r^2) + ((4*C2*y^3)/rbc.r^4) + ((4*C2*x^2*y)/rbc.r^4) )*(1-((x^2+y^2)/(rbc.r^2))));
            
            % Preset matrices for calculated perpendicular vectors
            ln_setX = NaN(size(p_set,1),size(p_set,2));
            ln_setY = NaN(size(p_set,1),size(p_set,2));
            ln_setZ = NaN(size(p_set,1),size(p_set,2));
            ln_setVX = NaN(size(p_set,1),size(p_set,2));
            ln_setVY = NaN(size(p_set,1),size(p_set,2));
            ln_setVZ = NaN(size(p_set,1),size(p_set,2));

            for i = 1:size(p_set,1)
                for j = 1:size(p_set,2)       
                    p = Point(p_set.X(i,j),p_set.Y(i,j),p_set.Z(i,j));
                    
                    % Calculate the normal vector by evaulating the
                    % gradient at the point of interest.
                    normal_vector = Vector(p.X,p.Y,p.Z,-rbc_x(p.X,p.Y),-rbc_y(p.X,p.Y),sign(p.Z));

                    ln_setX(i,j) = normal_vector.X; ln_setY(i,j) = normal_vector.Y; ln_setZ(i,j) = normal_vector.Z;
                    ln_setVX(i,j) = normal_vector.Vx; ln_setVY(i,j) = normal_vector.Vy; ln_setVZ(i,j) = normal_vector.Vz;
                end
            end
            
            ln_set = Vector(ln_setX,ln_setY,ln_setZ,ln_setVX,ln_setVY,ln_setVZ);
            
            ln_set = SLine(Point(ln_set.X,ln_set.Y,ln_set.Z),Point(ln_set.X+ln_set.Vx,ln_set.Y+ln_set.Vy,ln_set.Z+ln_set.Vz));

            % Matrices for line set
            p1X = NaN(size(ln_set,1),size(ln_set,2));
            p1Y = NaN(size(ln_set,1),size(ln_set,2));
            p1Z = NaN(size(ln_set,1),size(ln_set,2));
            p2X = NaN(size(ln_set,1),size(ln_set,2));
            p2Y = NaN(size(ln_set,1),size(ln_set,2));
            p2Z = NaN(size(ln_set,1),size(ln_set,2));
            
            % Reverse-rotate and Reverse-translate the perpendicular line set
            for i = 1:numel(ln_set)
                p1 = R * [ln_set.p1.X(i); ln_set.p1.Y(i); ln_set.p1.Z(i)];
                p2 = R * [ln_set.p2.X(i); ln_set.p2.Y(i); ln_set.p2.Z(i)];
                p1X(i) = p1(1); p1Y(i) = p1(2); p1Z(i) = p1(3);
                p2X(i) = p2(1); p2Y(i) = p2(2); p2Z(i) = p2(3);
            end

            ln_set = SLine(Point(p1X,p1Y,p1Z),Point(p2X,p2Y,p2Z));
            ln_set = ln_set.translate(rbc.c);

            % Scale back to original size
            ln_set.p1.X = ln_set.p1.X * (10^ord);
            ln_set.p1.Y = ln_set.p1.Y * (10^ord);
            ln_set.p1.Z = ln_set.p1.Z * (10^ord);
            ln_set.p2.X = ln_set.p2.X * (10^ord);
            ln_set.p2.Y = ln_set.p2.Y * (10^ord);
            ln_set.p2.Z = ln_set.p2.Z * (10^ord);
    
            % Perpendicular line
            ln = ln_set;
        end
        function pl = tangentplane(rbc,p_set)
            % TANGENTPLANE Plane tangent to RBC passing by point set
            % 
            % PL = TANGENTPLANE(RBC,P_SET) calculates plane set PL tangent to 
            %    RBC and passing by point set P_SET.
            % 
            % See also RBC, Point, Plane.          

            Check.isa('P must be a Point',p_set,'Point')

            pl = Plane.perpto(rbc.perpline(p_set),p_set);
        end
        % function angles = orientation(rbc,reverse)
        %     % ORIENTATION set of angles that denote the orientation of a RBC
        %     % 
        %     % ANGLES = ORIENTATION(RBC,REVERSE) calculates three angles that RBC 
        %     %  makes with each axis. If reverse is true, it will return the
        %     %  angles that a default RBC makes with each reference axis of
        %     %  the given RBC.
        %     % 
        %     % See also RBC.
        % 
        %     if reverse
        %         or_init = [rbc.or.Vx rbc.or.Vy rbc.or.Vz]; % RBC orientation
        %         or_final = [0 0 1]; % Defult orientation 
        %     else
        %         or_init = [0 0 1]; 
        %         or_final = [rbc.or.Vx rbc.or.Vy rbc.or.Vz];
        %     end
        % 
        %     rot_angle = acos(dot(or_init, or_final)); % Angle of rotation
        % 
        %     k = cross(or_init, or_final); % Rotation axis
        % 
        %     q = [cos(rot_angle/2), k * sin(rot_angle/2)]; % Rotation quaternion
        % 
        %     x_angle = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)*q(2) + q(3)*q(3)));
        %     y_angle = asin(2*(q(1)*q(3) - q(4)*q(2)));
        %     z_angle = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(3)*q(3) + q(4)*q(4)));
        % 
        %     % Rotation angles
        %     angles = [x_angle y_angle z_angle];
        % end
    end 
end