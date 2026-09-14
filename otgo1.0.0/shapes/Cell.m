classdef Cell < Superficies
    % Cell < Superficies : Arbitrarily generated Cell or Vesicle in 3D
    %   A Cell is defined by its number of discretization points p, shape shp, 
    %   center c, orientations or1, or2, or3 and order of magnitude mag.
    %   p must be a real number, shp must be a valid shape type, c must be a point,
    %   or1, or2, or3 must be vectors, vol must be a real number and mag must be a real number.
    %
    % Cell properties:
    %   p       - number of discretization points
    %   shp     - description of the shape
    %   c       - center (Point)
    %   or1     - reference axis of x direction with respect to the particle reference frame (Vector)
    %   or2     - reference axis of y direction with respect to the particle reference frame (Vector)
    %   or3     - reference axis of z direction with respect to the particle reference frame (Vector)
    %   cart    - cartesian coordinates of the cell 
    %   shc     - spherical harmonic coefficients
    %   geoProp - geometrical properties of the cell such as divergence,
    %             gradient, normal vector, linear curvature, laplacian,
    %             mean curvature, bending and tension
    %   vol     - volume of the cell
    %   surArea - surface area of the cell
    %   mag     - order of magnitude
    %
    % Cell methods:
    %   Cell                -   constructor
    %   plot                -   plots cell in 3D
    %   disp                -   prints cella
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of cells (=1)
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   tangent plane set passing by point set
    %   orientation         -   angle of orientation of the cell
    %   t_limits            -   limits for the intersection point on the line
    %   Volume              -   volume of the cell
    %   intpointaccuracy    -   how accurate the intersection point computation is
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
    %   Date: 2023/08/23

    properties
        p        % number of discretization points
        shp      % description of the shape
        c        % center (Point)
        or1      % reference axis of x direction with respect to the particle reference frame (Vector)
        or2      % reference axis of y direction with respect to the particle reference frame (Vector)
        or3      % reference axis of z direction with respect to the particle reference frame (Vector)
        cart     % cartesian coordinates of the cell 
        shc      % spherical harmonic coefficients
        geoProp  % geometrical properties of the cell
        vol      % volume of the cell
        surArea  % surface area of the cell
        mag      % order of magnitude of cell dimensions
    end
    methods
        function obj = Cell(p,shp,c,or1,or2,or3,vol,mag)
            % CELL(p,shp,c,or1,or2,or3,vol,mag) constructs a Cell
            %   with number of discretization points p, cell type shp,
            %   center c, orientation or and a order of magnitude mag.
            %   p must be a real positive number, shp must be a valid shape
            %   type, c must be a point, or1, or2, or3 must be a vectors, 
            %   vol must be a positive real number and mag must be an integer.
            %
            % See also Cell, Point, Vector.

            Check.isreal('p must be real number greater than 0',p,'>',0)
            Check.isa('c must be a Point',c,'Point')
            Check.isa('or1 must be a Vector',or1,'Vector')
            Check.isa('or2 must be a Vector',or2,'Vector')
            Check.isa('or3 must be a Vector',or3,'Vector')
            Check.isreal('The order of magnitude must be a real number',mag)      
                
            try
                shape = boundary(p,shp);

                V = volume_from_mesh(shape.cart.x, shape.cart.y, shape.cart.z);

                if strcmp(vol, 'vol')
                    k = 1;
                    obj.vol = V;
                else
                    Check.isreal('vol must be real number greater than 0',vol,'>',0)
                    k = (vol/V)^(1/3);
                    obj.vol = vol;
                end
                   
                obj.p = p;
                obj.shp = shp;
                obj.c = c;
                obj.or1 = Vector(0,0,0,or1.Vx,or1.Vy,or1.Vz).normalize;
                obj.or2 = Vector(0,0,0,or2.Vx,or2.Vy,or2.Vz).normalize;
                obj.or3 = Vector(0,0,0,or3.Vx,or3.Vy,or3.Vz).normalize;
                obj.cart = d3Vec(shape.cart.x*k,shape.cart.y*k,shape.cart.z*k);
                obj.shc = shape.shc;
                obj.geoProp = shape.geoProp;
                obj.surArea = surface_area_from_mesh(shape.cart.x, shape.cart.y, shape.cart.z) * (k^2);
                obj.mag = mag;
            catch
                error('Error providing valid cell properties!')
            end
        end
        function h = plot(cell,varargin)
            % PLOT Plots Cell in 3D
            %
            % H = PLOT(CELL) plots the cell in 3D. It returns a
            %   graphic handler to the plotted cell.
            %
            % H = PLOT(CELL,'Range',N) sets the number of discretization points 
            %   of the plotting. N = cell.p (default) corresponds to default 
            %   number of discretization points of the shape.
            %
            % H = PLOT(CELL,'Scale',S) rescales the coordinates of
            %   by S before plotting the cell. S=1 by default. 
            %
            % H = PLOT(CELL,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(CELL,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Cell, surf.
            
            Xin = vecForm(cell.cart); % Cartesian coordinates of the shape

            N = (round(sqrt(2*size(Xin,1)/3+1))-1)/2; % Number of discretization points
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    N = varargin{n+1};
                end
            end
            
            S = 10^cell.mag; % Scale
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

            Xin = reshape(Xin,[], 3 * size(Xin,2));  

            X = Xin(:,1:3);
    
            % Putting the caps
            [u, v] = parDomain(N);

            u = reshape(u, N+1,2*N);
            v = reshape(v, N+1,2*N);
            u = [pi*ones(1,2*N); u; 0*ones(1,2*N)];
            v = [v; v(1:2,:)];
        
            for i = 1:3
              XCap(:,i) = sumBasis( shAna(X(:,i)), 'Ynm', true, u(:), v(:));
            end

            x = reshape(XCap(:,1), N+3, 2*N);
            y = reshape(XCap(:,2), N+3, 2*N);
            z = reshape(XCap(:,3), N+3, 2*N);

            % Connect two ends of the grid for a continious surface
            x = [x x(:,1)];
            y = [y y(:,1)];
            z = [z z(:,1)];
            
            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];
                    
            X_rot = zeros(size(x,1),size(x,2));
            Y_rot = zeros(size(y,1),size(y,2));
            Z_rot = zeros(size(z,1),size(z,2));
            
            % Rotation
            for i = 1:size(x,1)*size(x,2)
                vertex = [x(i); y(i); z(i)];
                rot_vertex = R * vertex;
                X_rot(i) = rot_vertex(1);
                Y_rot(i) = rot_vertex(2);
                Z_rot(i) = rot_vertex(3);
            end

            % Translation
            X = cell.c.X + X_rot*S;
            Y = cell.c.Y + Y_rot*S;
            Z = cell.c.Z + Z_rot*S;

            % ht = mesh(X, Y, Z, C*ones(size(X)), ...
            %     'EdgeColor', 'k', 'EdgeAlpha', 0.4, 'FaceColor', 'r', 'FaceAlpha', 0.45, ...
            %     'FaceLighting', 'none', 'EdgeLighting', 'none');

            ht = mesh(X, Y, Z, C*ones(size(X)), ...
                'EdgeColor', 'k', 'EdgeAlpha', 1, 'FaceColor', 'r', 'FaceAlpha', 1);
                
            hold on;

            % draw the center
            % set(gca, 'SortMethod', 'childorder')
            % plot3(cell.c.X, cell.c.Y, cell.c.Z, 'o', ...
            %     'MarkerEdgeColor', 'k', ...
            %     'MarkerFaceColor', [1 1 0], ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 1.2);

            % ── Axis tripod ──────────────────────────────────────────────
            % L = 0.4 * max([range(X(:)), range(Y(:)), range(Z(:))]);
            % 
            % cx = cell.c.X;  cy = cell.c.Y;  cz = cell.c.Z;
            % 
            % % X arrow (blue)
            % xOffset = 0.7 * max(X(:)) - cx;
            % quiver3(cx + xOffset, cy, cz, L, 0, 0, 0, 'Color', [0.00 0.30 0.80], ...
            %     'LineWidth', 2.5, 'MaxHeadSize', 0.6);
            % text(cx + xOffset + 1.2*L, cy, cz, 'x', 'FontSize', 26, 'FontWeight', 'bold', ...
            %     'Color', [0.00 0.30 0.80], 'HorizontalAlignment', 'center');
            % 
            % % Y arrow (green)
            % yOffset = 1 * max(Y(:)) - cy;
            % quiver3(cx, cy + yOffset, cz, 0, L, 0, 0, 'Color', [0.10 0.60 0.10], ...
            %     'LineWidth', 2.5, 'MaxHeadSize', 0.6);
            % text(cx, cy + yOffset + 1.4*L, cz, 'y', 'FontSize', 26, 'FontWeight', 'bold', ...
            %     'Color', [0.10 0.60 0.10], 'HorizontalAlignment', 'center');
            % 
            % % Z arrow (yellow)
            % zOffset = 0.6 * max(Z(:)) - cz;
            % quiver3(cx, cy, cz + zOffset, 0, 0, L, 0, 'Color', [0.85 0.75 0.00], ...
            %     'LineWidth', 2.5, 'MaxHeadSize', 0.6);
            % text(cx, cy, cz + zOffset + 1.4*L, 'z', 'FontSize', 26, 'FontWeight', 'bold', ...
            %     'Color', [0.85 0.75 0.00], 'HorizontalAlignment', 'center');
            
            % ── Omega arc arrow (rotation about x-axis) ──────────────────
            % arcR   = 1.2 * L;          % increased from 0.85
            % arcPhi = linspace(pi/4, 5*pi/4, 60);
            % arc_y  = cy + arcR * cos(arcPhi);
            % arc_z  = cz + arcR * sin(arcPhi);
            % arc_x  = (cx + 1.6*L) * ones(size(arcPhi));  % increased from 1.2
            % 
            % plot3(arc_x, arc_y, arc_z, 'k-', 'LineWidth', 2.0);
            % 
            % % Arrowhead at the end of the arc
            % dphi  = arcPhi(end) - arcPhi(end-1);
            % dy    = -arcR * sin(arcPhi(end)) * dphi;
            % dz    =  arcR * cos(arcPhi(end)) * dphi;
            % quiver3(arc_x(end), arc_y(end), arc_z(end), 0, dy, dz, 3, ...
            %     'k', 'LineWidth', 2.0, 'MaxHeadSize', 1.5);
            % 
            % % Omega label
            % text(cx + 1.6*L, cy + arcR*cos(arcPhi(30)), cz + arcR*sin(arcPhi(30)) + 0.6*L, ...
            %     '\Omega', 'FontSize', 26, 'FontWeight', 'bold', ...
            %     'Interpreter', 'tex', 'HorizontalAlignment', 'center');
                        
            % ── Clean up axes (match RBC / Fig. 1) ───────────────────────
            grid off
            axis equal
            axis off
            view(3)
            % Interactive exploration (commented for clean paper insets):
            % xlabel('x')
            % ylabel('y')
            % zlabel('z')
            % axis on

            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'precision') && ~strcmpi(varargin{n},'scale') && ~strcmpi(varargin{n},'colorlevel') && ~strcmpi(varargin{n},'range')
                    set(ht,varargin{n},varargin{n+1});
                end
            end
            
            % Output if needed
            if nargout>0
                h = ht;
            end
        end
        function disp(cell)
            % DISP Prints the cell
            %
            % See also Cell.

            fprintf('Cell properties:\n');
            fprintf('Precision: %g\n', cell.p);
            fprintf('Center:\n');
            disp(cell.c);
            fprintf('Reference axis (X axis of particle space):\n');
            disp(cell.or1);
            fprintf('Reference axis (Y axis of particle space):\n');
            disp(cell.or2);
            fprintf('Reference axis (Z axis of particle space):\n');
            disp(cell.or3);
            fprintf('Cartesian coordinates of points:\n');
            disp(cell.cart);
            fprintf('Spherical harmonics coefficients:\n');
            disp(cell.shc);
            fprintf('Geometrical properties:\n');
            disp(cell.geoProp);
            fprintf('Order of Magnitude: %g\n', cell.mag);
            fprintf('Cell Volume: %g\n', cell.vol);
        end
        function cell_t = translate(cell,dp)
            % TRANSLATE 3D translation of the cell
            %
            % CELL_T = TRANSLATE(CELL,dP) translates the cell by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Cell, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            cell_t = cell;
            cell_t.c = cell_t.c.translate(dp);
        end
        function cell_r = xrotation(cell,phi)
            % XROTATION Rotation around x-axis of the cell
            %
            % CELL_R = XROTATION(CELL,phi) rotates the cell around x-axis 
            %   by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            cell_r = cell;
            cell_r.c = cell_r.c.xrotation(phi);
            cell_r.or1 = cell_r.or1.xrotation(phi);
            cell_r.or2 = cell_r.or2.xrotation(phi);
            cell_r.or3 = cell_r.or3.xrotation(phi);
        end
        function cell_r = yrotation(cell,phi)
            % XROTATION Rotation around y-axis of the cell
            %
            % CELL_R = XROTATION(CELL,phi) rotates the cell around y-axis 
            %   by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            cell_r = cell;
            cell_r.c = cell_r.c.yrotation(phi);
            cell_r.or1 = cell_r.or1.yrotation(phi);
            cell_r.or2 = cell_r.or2.yrotation(phi);
            cell_r.or3 = cell_r.or3.yrotation(phi);
        end
        function cell_r = zrotation(cell,phi)
            % XROTATION Rotation around z-axis of the cell
            %
            % CELL_R = XROTATION(CELL,phi) rotates the cell around x-axis 
            %   by an angle phi [rad].
            %
            % See also Cell.

            Check.isreal('The rotation angle phi must be a real number',phi)

            cell_r = cell;
            cell_r.c = cell_r.c.zrotation(phi);
            cell_r.or1 = cell_r.or1.zrotation(phi);
            cell_r.or2 = cell_r.or2.zrotation(phi);
            cell_r.or3 = cell_r.or3.zrotation(phi);
        end
        function n = numel(cell)
            % NUMEL Number of cell set
            %
            % N = NUMEL(CELL) number of cell(s).
            %
            % See also Cell.

            n = numel(cell.c);
        end
        function s = size(cell,varargin)
            % SIZE Size of the cell set
            % 
            % S = SIZE(CELL) returns a two-element row vector with the number 
            %   of rows and columns in the cell set. Default value is (1,1)
            %   if not provided as a set.
            %
            % S = SIZE(CELL,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the cell set.
            %
            % See also Cell, RBC.
         
            if ~isempty(varargin)
                s = cell.c.size(varargin{1});
            else
                s = cell.c.size();
            end
        end

        % Use the below algoithm if a completely numerical soultion is
        % desired. Typically very fast and within some reasonable accuracy
        % level. Suitable with any cell that can be expressed with grid
        % poitns. (Adjusted to only work with rays)

        function p = intersectionpoint(cell,d,n)
            % INTERSECTIONPOINT Intersection point between the cell and line/vector/ray
            %
            % P = INTERSECTIONPOINT(CELL,D,N,NPHI) calculates intersection points 
            %   between a set of lines (or vectors) D and a cell. The
            %   numeric algorithm allows to calculate the interactions
            %   where the line enters the cell for the first time and the second one where the ray
            %   leaves the cell. So n can either be 1 or 2. Nphi is the
            %   number of radial divisions in ray discretization, this
            %   value is used to gemoterically merge the two rays that has
            %   same direction but different polarization.
            % 
            %   Number of grid points on the ray in each of which, the distance with
            %   the cell is evaluated can be changed. However it is not
            %   recommended to be lower than 20. It's also recorded that after
            %   some level, increasing this value further do not increase
            %   accuracy. Default value is 25 if not changed.
            %   
            %   If D does not intersect with RBC, the coordinates of P are NaN.
            % 
            % See also Cell, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            Check.isinteger('A must be either 1 or 2',n,'>=',1,'<=',2)

            if isa(d,'SLine')
                ln_set = d;
            else
                ln_set = d.toline();                
            end

            size_r = size(ln_set);

            min_max = [min(cell.cart.x) min(cell.cart.y) min(cell.cart.z);
                       max(cell.cart.x) max(cell.cart.y) max(cell.cart.z)];

            cart_vec = [cell.cart.x'; cell.cart.y'; cell.cart.z'];
            cart_size = size(cart_vec,2);
            mg = cell.mag;
            center = [cell.c.X; cell.c.Y; cell.c.Z];

            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];

            % Adjust the line so that it intersects with the cell
            % that is centered at the origin and with a default
            % orientation (Because this is what the mathematical expressions denote)
            offset = center*(10^-mg);
            
            Nphi = 20;
            n_div = 20;
            tol = 0.5; % The limit of the unit distance between the cell and the int point 

            % Matrices for calculated intersection points
            pX = NaN(size_r); pY = NaN(size_r); pZ = NaN(size_r);

            % Set the line set to unit scale
            maxord = max([floor(log10(abs(ln_set.p1.X(:)))); 
              floor(log10(abs(ln_set.p1.Y(:)))); 
              floor(log10(abs(ln_set.p1.Z(:))))]);
            if max(maxord) <= -6
                ln_set.p1 = ln_set.p1 * 10^6;
                ln_set.p2 = ln_set.p2 * 10^6;
            end
            
            % Stack the coordinates of p1 and p2 into 3xN matrices
            P1 = [ln_set.p1.X(:)' - offset(1); 
                  ln_set.p1.Y(:)' - offset(2); 
                  ln_set.p1.Z(:)' - offset(3)];

            P2 = [ln_set.p2.X(:)' - offset(1); 
                  ln_set.p2.Y(:)' - offset(2); 
                  ln_set.p2.Z(:)' - offset(3)];

            % Rotate
            P1_rot = R' * P1;
            P2_rot = R' * P2;

            % Extract components as row vectors
            p1X = P1_rot(1, :); 
            p1Y = P1_rot(2, :); 
            p1Z = P1_rot(3, :);

            p2X = P2_rot(1, :); 
            p2Y = P2_rot(2, :); 
            p2Z = P2_rot(3, :);

            idx_list = [];
            for start = 1:(2*Nphi):numel(ln_set)
                idx_list = [idx_list, start:start+Nphi-1];
            end

            for i = idx_list

                ln_p1 = [p1X(i); p1Y(i); p1Z(i)];
                ln_p2 = [p2X(i); p2Y(i); p2Z(i)];
                lnc = ln_p2 - ln_p1;

                [t_lower, t_upper] = t_limits(ln_p1,lnc,min_max);

                if n == 1
                    t_upper = (t_upper + t_lower) / 2; 
                elseif n == 2
                    % t_lower = (t_upper + t_lower) / 2; 
                    t_lower = 0.25;
                end

                t = linspace(t_lower,t_upper,n_div);
                d = ray_cell_distance(t, ln_p1, lnc, cart_size, cart_vec, n_div);
                [min_dists, ~] = min(d, [], 2);

                try
                    idx = find(islocalmin(min_dists),1);
                    %[~, k] = min(min_dists(idx));
                    %idx = idx(k);

                    if idx >= 2; t_lower = t(idx-1); end
                    if idx <= n_div-1; t_upper = t(idx+1); end

                    t = linspace(t_lower,t_upper,n_div);
                    d = ray_cell_distance(t, ln_p1, lnc, cart_size, cart_vec, n_div);
                    [min_dists, ~] = min(d, [], 2);
                    second = true;
                catch
                end
                
                if true
                    [~, closest_point_index] = min(min_dists);
    
                    if closest_point_index == 1; closest_point_index = 2; end
                    if closest_point_index == n_div; closest_point_index = n_div-1; end
    
                    t = linspace(t(closest_point_index-1),t(closest_point_index+1),n_div);
                    d = ray_cell_distance(t, ln_p1, lnc, cart_size, cart_vec, n_div);
    
                    [min_dists, ~] = min(d, [], 2);
                end

                [val, closest_point_index] = min(min_dists);

                try
                    if val < tol
                        t = t(closest_point_index);
                        int_p = [ln_p1(1) + t * lnc(1); 
                                 ln_p1(2) + t * lnc(2); 
                                 ln_p1(3) + t * lnc(3)];
                    else
                        int_p = [NaN; NaN; NaN];
                    end
                catch
                    int_p = [NaN; NaN; NaN];
                end

                int_p = (R * int_p) + offset;
                scaled = int_p * 10^mg;
                pX(i) = scaled(1);
                pX(i+Nphi) = scaled(1);
                pY(i) = scaled(2);
                pY(i+Nphi) = scaled(2);
                pZ(i) = scaled(3);
                pZ(i+Nphi) = scaled(3);
            end

            p = Point(pX, pY, pZ); 
        end
        function ln = perpline(cell,p_set)
            % PERPLINE Line perpendicular to cell passing by point set
            % 
            % LN = PERPLINE(CELL,P-SET) calculates the line set LN perpendicular 
            %   to the cell CELL and passing by the point set P.
            % 
            % See also RBC, Point, SLine.

            Check.isa('p_set must be a Point',p_set,'Point')

            % Create X_mesh and nor_mesh matrices
            cell_mesh = [cell.cart.x'; cell.cart.y'; cell.cart.z'];
            nor_mesh = [cell.geoProp.nor.x'; cell.geoProp.nor.y'; cell.geoProp.nor.z'];

            Nphi = 20;
            mg = cell.mag;
            center = [cell.c.X; cell.c.Y; cell.c.Z];
            size_p = size(p_set);

            % Preset matrices for calculated perpendicular vectors
            p1X = NaN(size_p); p1Y = NaN(size_p); p1Z = NaN(size_p);
            p2X = NaN(size_p); p2Y = NaN(size_p); p2Z = NaN(size_p);
                
            % Scale point set to unit sizes

            % Adjust the point set so that they are on the cell
            % that is centered at the origin and with a default
            % orientation (because this is what the mathematical expression denotes).
            P = [(p_set.X(:)'-center(1))*(10^-mg);
                 (p_set.Y(:)'-center(2))*(10^-mg); 
                 (p_set.Z(:)'-center(3))*(10^-mg)];

            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];

            P_rot = R' * P;

            idx_list = [];
            for start = 1:(2*Nphi):numel(p_set)
                idx_list = [idx_list, start:start+Nphi-1];
            end

            for i = idx_list

                query_point = [P_rot(1,i); P_rot(2,i); P_rot(3,i)];

                % Calculate the distances from the query_point to all mesh points.
                distances = sqrt(sum((cell_mesh - query_point).^2, 1));

                % Find the indices of the nearest mesh points.
                [~, nearest_indices] = sort(distances);

                num_nearest = 10;
                nearest_indices = nearest_indices(1:num_nearest);

                % Use the nearest mesh points for interpolation (e.g., using weighted average).
                weights = 1 ./ distances(nearest_indices);
                normalized_weights = weights / sum(weights);

                % Interpolate the normal vectors using the weights.
                interpolated_normal = sum(nor_mesh(:, nearest_indices) .* normalized_weights, 2);

                % Normalize the interpolated normal vector.
                normal = interpolated_normal / norm(interpolated_normal);

                % ln_setX(i) = query_point(1); ln_setY(i) = query_point(2); ln_setZ(i) = query_point(3);
                % ln_setVX(i) = normal(1); ln_setVY(i) = normal(2); ln_setVZ(i) = normal(3);

                p1 = (R * [query_point(1); query_point(2); query_point(3)]) * 10^mg + center;
                p2 = (R * [query_point(1)+normal(1); query_point(2)+normal(2); query_point(3)+normal(3)]) * 10^mg + center;

                p1X(i) = p1(1); p1Y(i) = p1(2); p1Z(i) = p1(3);
                p1X(i+Nphi) = p1(1); p1Y(i+Nphi) = p1(2); p1Z(i+Nphi) = p1(3);
                p2X(i) = p2(1); p2Y(i) = p2(2); p2Z(i) = p2(3);
                p2X(i+Nphi) = p2(1); p2Y(i+Nphi) = p2(2); p2Z(i+Nphi) = p2(3);
            end
    
            % Perpendicular line
            ln = SLine(Point(p1X,p1Y,p1Z),Point(p2X,p2Y,p2Z));
        end
        function pl = tangentplane(cell,p_set)
            % TANGENTPLANE Plane tangent to cell passing by point set
            % 
            % PL = TANGENTPLANE(CELL,P_SET) calculates plane set PL tangent to 
            %    cell and passing by point set P_SET.
            % 
            % See also Cell, Point, Plane.          

            Check.isa('P must be a Point',p_set,'Point')

            pl = Plane.perpto(cell.perpline(p_set),p_set);
        end

        function angles = orientation(cell,reverse)
            % ORIENTATION set of angles that denote the orientation of a cell
            % 
            % ANGLES = ORIENTATION(CELL,REVERSE) calculates three angles that cell 
            % makes with each axis. If reverse is true, it will return the
            % angles that a default cell makes with each reference axis of
            % the cell CELL.
            % 
            % See also Cell.

            if reverse
                or_init = [cell.or3.Vx cell.or3.Vy cell.or3.Vz]; % Cell orientation
                or_final = [0 0 1]; % Default orientation
            else
                or_init = [0 0 1];
                or_final = [cell.or3.Vx cell.or3.Vy cell.or3.Vz];
            end
            
            rot_angle = acos(dot(or_init, or_final)); % Angle of rotation
           
            k = cross(or_init, or_final); % Rotation axis
            
            q = [cos(rot_angle/2), k * sin(rot_angle/2)]; % Rotation quaternion
            
            x_angle = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)*q(2) + q(3)*q(3)));
            y_angle = asin(2*(q(1)*q(3) - q(4)*q(2)));
            z_angle = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(3)*q(3) + q(4)*q(4)));
            
            % Rotation angles
            angles = [x_angle y_angle z_angle];
        end
        function vol = Volume(cell)
            [~, vol] = convhull(cell.cart.x,cell.cart.y,cell.cart.z);
        end
        function t_sol = newton_solve(eq, t_initial, tol, max_iter, lim_lower, lim_upper)
            t_sol = t_initial;
            
            for i = 1:max_iter
                f_val = eq(t_sol);  % Evaluate the function at the current guess
                f_prime = (eq(t_sol + tol) - f_val) / tol;  % Numerical derivative
                
                % Update t_sol using the Newton-Raphson formula
                t_sol_new = t_sol - f_val / f_prime;
                
                % Check if the new solution is within the specified range
                if t_sol_new < lim_lower
                    t_sol_new = lim_lower;  % Clamp to lower bound
                elseif t_sol_new > lim_upper
                    t_sol_new = lim_upper;  % Clamp to upper bound
                end
                
                % Check if the solution has converged
                if abs(t_sol_new - t_sol) < tol
                    break;
                end
                
                % Update the current guess
                t_sol = t_sol_new;
            end
            
            % If max iterations are reached without convergence, display a warning
            if i == max_iter
                warning('Maximum number of iterations reached without convergence.');
            end
        end
    end
end