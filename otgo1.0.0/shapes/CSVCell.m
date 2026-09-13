classdef CSVCell < Superficies
    % CSVCell < Superficies : Numeric cell loaded from a csv file.
    %   A CSVCell is defined by its center c, orientations or1, or2, or3, volume 
    %   vol and order of magnitude mag.
    %   c must be a point, or1, or2, or3 must be vectors, vol must be a real number 
    %   and mag must be a real number.
    %
    % CSVCell properties:
    %   c       - center (Point)
    %   or1     - reference axis of x direction with respect to the particle reference frame (Vector)
    %   or2     - reference axis of y direction with respect to the particle reference frame (Vector)
    %   or3     - reference axis of z direction with respect to the particle reference frame (Vector)
    %   cart    - cartesian coordinates of the cell
    %   alpha   - alpha radius
    %   vol     - volume of the cell
    %   mag     - order of magnitude
    %
    % CSVCell methods:
    %   CSVCell             -   constructor
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
    %   Volume              -   volume of the cell
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
    %   Date: 2023/12/13

    properties
        c        % center (Point)
        or1      % reference axis of x direction with respect to the particle reference frame (Vector)
        or2      % reference axis of y direction with respect to the particle reference frame (Vector)
        or3      % reference axis of z direction with respect to the particle reference frame (Vector)
        cart     % cartesian coordinates of the cell 
        ashp     % alpha shape
        vol      % volume of the cell
        mag      % order of magnitude of cell dimensions
    end
    methods
        function obj = CSVCell(c,or1,or2,or3,vol,mag,data)
            % CSVCELL(c,or1,or2,or3,vol,mag) constructs a Cell
            %   with number center c, orientations or1, or2, or3, order of 
            %   magnitude mag and volume vol. c must be a point, or1, or2, 
            %   or3 must be a vectors, mag must be an integer and vol must be 
            %   a positive real number (if not "vol" which creates a cell without 
            %   changing the original volume) 
            %
            % See also CSVCell, Point, Vector.

            Check.isa('c must be a Point',c,'Point')
            Check.isa('or1 must be a Vector',or1,'Vector')
            Check.isa('or2 must be a Vector',or2,'Vector')
            Check.isa('or3 must be a Vector',or3,'Vector')
            Check.isreal('The order of magnitude must be a real number',mag)      
                
            pos = table2array(readtable(data));
            pos = pos(:,24:26);
            
            shp = alphaShape(pos,0.3);
            %alpha = criticalAlpha(shp,'one-region') + 0.1;
            %shp = alphaShape(pos,alpha);

            try
                k = (vol/shp.volume)^(1/3);
            catch
                k = 1;
            end

            try
                obj.c = c;
                obj.or1 = Vector(0,0,0,or1.Vx,or1.Vy,or1.Vz).normalize;
                obj.or2 = Vector(0,0,0,or2.Vx,or2.Vy,or2.Vz).normalize;
                obj.or3 = Vector(0,0,0,or3.Vx,or3.Vy,or3.Vz).normalize;

                x_mean = mean(d3Vec(pos(:,1)*k,pos(:,2)*k,pos(:,3)*k).x);

                obj.cart = d3Vec(pos(:,1)*k-x_mean,pos(:,2)*k,pos(:,3)*k);
                obj.ashp = shp;
                obj.vol = vol;
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
            % See also CCSVell, surf.
            
            S = 10^cell.mag; % Scale
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = S*varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                end
            end

            % Color level
            % C = 0;
            % for n = 1:2:length(varargin)
            %     if strcmpi(varargin{n},'colorlevel')
            %         C = varargin{n+1};
            %         Check.isreal('The scaling factor must be a real number',C)
            %     end
            % end

            x = cell.cart.x;
            y = cell.cart.y;
            z = cell.cart.z;

            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];

            % Calculate R using the following way if the cell will have a
            % single reference vector (along +z) denoting its orientation.

            % angles = orientation(cell,false);
            % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 

            % Set of rotation matrices
            % Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
            % Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
            % Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];

            %R = Rz * Ry * Rx;
                 
            X_rot = zeros(size(x,1),size(x,2));
            Y_rot = zeros(size(y,1),size(y,2));
            Z_rot = zeros(size(z,1),size(z,2));

            % Rotation
            for i = 1:size(x,1)*size(x,2)
                vertex = [x(i,:); y(i,:); z(i,:)];
                rot_vertex = R * vertex;
                X_rot(i) = rot_vertex(1);
                Y_rot(i) = rot_vertex(2);
                Z_rot(i) = rot_vertex(3);
            end

            % Translation
            X = cell.c.X + X_rot*S;
            Y = cell.c.Y + Y_rot*S;
            Z = cell.c.Z + Z_rot*S;
            
            shp = alphaShape(X,Y,Z,0.45*S);
            %shp = alphaShape([X Y Z],criticalAlpha(shp,'one-region') + 0.1);
        
            %ht = plot(shp, 'FaceAlpha', 0, 'EdgeColor', 'red');
            %ht = plot(shp, 'FaceAlpha', 1, 'EdgeColor', 'red','FaceColor','white');
            ht = plot(shp, 'FaceAlpha', 1, 'EdgeColor', 'black', 'FaceColor', 'red');

            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis('equal')
            grid on
            view(3)

            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'scale')
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
            fprintf('Center:\n');
            disp(cell.c);
            fprintf('Reference axis (X axis of particle space):\n');
            disp(cell.or1);
            fprintf('Reference axis (Y axis of particle space):\n');
            disp(cell.or2);
            fprintf('Reference axis (Z axis of particle space):\n');
            disp(cell.or3);
            fprintf('Order of Magnitude: %g\n', cell.mag);
            fprintf('Volume: %g\n', cell.vol);
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
            % See also CSVCell.
         
            if ~isempty(varargin)
                s = cell.c.size(varargin{1});
            else
                s = cell.c.size();
            end
        end

        % Use the below solver if a completely numeric solution is desired
        % function p = intersectionpoint(cell,d,n,varargin)
        %     % INTERSECTIONPOINT Intersection point between the cell and line/vector/ray
        %     %
        %     % P = INTERSECTIONPOINT(CELL,D,N,FIRST) calculates intersection points 
        %     %   between a set of lines (or vectors) D and a cell. The
        %     %   numeric algorithm allows to calculate the interactions
        %     %   where the line enters the cell for the first time and the second one where the ray
        %     %   leaves the cell for the first time. So n can either be 1 or 2. 
        %     % 
        %     %   Number of grid points on the ray where the distance with
        %     %   the cell is evaluated can be changed. However it is not
        %     %   recommended to be lower than 20. Its also recorded after
        %     %   some level, increasing this value further do not increase
        %     %   accuracy. Default value is 25 if not explicitly provided.
        %     %   
        %     %   If D does not intersect CELL, the coordinates of P are NaN.
        %     % 
        %     % See also CSVCell, Point, Vector, SLine, Ray.
        % 
        %     Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
        %     Check.isinteger('A must be either 1 or 2',n,'>=',1,'<=',2)
        % 
        %     if isa(d,'SLine')
        %         ln_set = d;
        %     else
        %         ln_set = d.toline();                
        %     end
        % 
        %     % Matrices for calculated intersection points
        %     pX = NaN(size(ln_set,1),size(ln_set,2));
        %     pY = NaN(size(ln_set,1),size(ln_set,2));
        %     pZ = NaN(size(ln_set,1),size(ln_set,2));
        % 
        %     % Matrices for line set
        %     p1X = NaN(size(ln_set,1),size(ln_set,2));
        %     p1Y = NaN(size(ln_set,1),size(ln_set,2));
        %     p1Z = NaN(size(ln_set,1),size(ln_set,2));
        %     p2X = NaN(size(ln_set,1),size(ln_set,2));
        %     p2Y = NaN(size(ln_set,1),size(ln_set,2));
        %     p2Z = NaN(size(ln_set,1),size(ln_set,2));
        % 
        %     % Set the line set to unit scale
        %     ordx = floor(log10(abs(ln_set.p1.X)));
        %     ordy = floor(log10(abs(ln_set.p1.Y)));
        %     ordz = floor(log10(abs(ln_set.p1.Z)));
        %     if max([ordx,ordy,ordz]) <= -6
        %         ln_set.p1 = ln_set.p1 * 10^6;
        %         ln_set.p2 = ln_set.p2 * 10^6;
        %     end
        % 
        %     R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
        %          cell.or1.Vy cell.or2.Vy cell.or3.Vy;
        %          cell.or1.Vz cell.or2.Vz cell.or3.Vz];
        % 
        %     % Calculate R using the following way if the cell will have a
        %     % single reference vector (along +z) denoting its orientation.
        % 
        %     % angles = orientation(cell,false);
        %     % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
        %     % 
        %     % % Set of rotation matrices
        %     % Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
        %     % Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
        %     % Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];
        %     % 
        %     % R = Rz * Ry * Rx;
        % 
        %     % Adjust the line so that it intersects with the cell
        %     % that is centered at the origin and with a default
        %     % orientation (Because this is what the mathematical expressions denote)
        %     ln_set = ln_set.translate(-1*cell.c*(10^-cell.mag));
        % 
        %     for i = 1:numel(ln_set)
        %         p1 = R' * [ln_set.p1.X(i); ln_set.p1.Y(i); ln_set.p1.Z(i)];
        %         p2 = R' * [ln_set.p2.X(i); ln_set.p2.Y(i); ln_set.p2.Z(i)];
        %         p1X(i) = p1(1); p1Y(i) = p1(2); p1Z(i) = p1(3);
        %         p2X(i) = p2(1); p2Y(i) = p2(2); p2Z(i) = p2(3);
        %     end
        % 
        %     ln_set = SLine(Point(p1X,p1Y,p1Z),Point(p2X,p2Y,p2Z));
        %     lnc_set = ln_set.p2 - ln_set.p1;
        % 
        %     n_div = 25; % default
        %     tol = 0.1;
        % 
        %     for i = 1:length(varargin)
        %         n_div = varargin{i};
        %     end
        % 
        %     cart_size = size(cell.cart.x,1);
        % 
        %     for i = 1:numel(ln_set)
        % 
        %         ln = SLine(Point(ln_set.p1.X(i),ln_set.p1.Y(i),ln_set.p1.Z(i)),Point(ln_set.p2.X(i),ln_set.p2.Y(i),ln_set.p2.Z(i)));
        %         lnc = Point(lnc_set.X(i),lnc_set.Y(i),lnc_set.Z(i));
        % 
        %         [t_lower, t_upper] = t_limits(cell,ln);
        % 
        %         if n == 1; t_upper = (t_upper + t_lower) / 2; end
        % 
        %         t = linspace(t_lower,t_upper,n_div);
        % 
        %         % Coordinates corresponding to different t values
        %         line_x = ln.p1.X + t' .* lnc.X;
        %         line_y = ln.p1.Y + t' .* lnc.Y;
        %         line_z = ln.p1.Z + t' .* lnc.Z;
        % 
        %         % Duplicate vectors to the appropriate size
        %         line_x = repmat(line_x, 1, cart_size);
        %         line_y = repmat(line_y, 1, cart_size);
        %         line_z = repmat(line_z, 1, cart_size);
        % 
        %         % Bring cartesian coordiantes of the cell to the
        %         % appropriate size
        %         cell_x = repmat(cell.cart.x, n_div, 1);
        %         cell_y = repmat(cell.cart.y, n_div, 1);
        %         cell_z = repmat(cell.cart.z, n_div, 1);
        % 
        %         dx = line_x - cell_x;
        %         dy = line_y - cell_y;
        %         dz = line_z - cell_z;
        % 
        %         d = (dx.^2 + dy.^2 + dz.^2).^(1/2);
        %         [min_dists, columnIndices] = min(d, [], 2);
        % 
        %         if n == 2
        %             localmax = find(islocalmax(min_dists));
        %             if isempty(localmax); localmax = 1; end
        %             rowIndex = min_dists == min(min_dists(localmax:end));
        %             if min(min_dists(localmax:end)) > tol
        %                 rowIndex = -1;
        %             end
        %         else
        %             rowIndex = min_dists == min(min_dists);
        %             if min(min_dists) > tol
        %                 rowIndex = -1;
        %             end
        %         end
        % 
        %         try
        %             minIndex = columnIndices(rowIndex);
        %             p_closest = Point(cell.cart.x(minIndex),cell.cart.y(minIndex),cell.cart.z(minIndex));
        % 
        %             % Tangent plane at the closest point
        %             tan_plane = Plane.perpto(perpline(cell,p_closest),p_closest);
        % 
        %             int_point = tan_plane.intersectionpoint(ln);
        %         catch
        %             int_point = NaN;
        %         end
        % 
        %         try                   
        %             pX(i) = int_point.X;
        %             pY(i) = int_point.Y;
        %             pZ(i) = int_point.Z;
        %         catch
        %             pX(i) = NaN;
        %             pY(i) = NaN;
        %             pZ(i) = NaN;
        %         end
        %     end
        %         int_point = Point(pX,pY,pZ);
        % 
        %         % Reverse-rotate and Reverse-translate the Intersection Point
        %         for i = 1:numel(int_point)
        %             vertex = [int_point.X(i); int_point.Y(i); int_point.Z(i)];
        %             rot_vertex = R * vertex;
        %             pX(i) = rot_vertex(1);
        %             pY(i) = rot_vertex(2);
        %             pZ(i) = rot_vertex(3);
        %         end
        % 
        %         int_point = Point(pX,pY,pZ);
        % 
        %         % Use if one reference axis is used to denote the shape
        %         % int_point = int_point.zrotation(-z_angle);
        %         % int_point = int_point.yrotation(-y_angle);
        %         % int_point = int_point.xrotation(-x_angle);
        % 
        %         int_point = int_point.translate(cell.c*(10^-cell.mag));
        % 
        %         % Scale back to original size
        %         int_point.X = int_point.X * 10^cell.mag;
        %         int_point.Y = int_point.Y * 10^cell.mag;
        %         int_point.Z = int_point.Z * 10^cell.mag;
        % 
        %         p = int_point; 
        % end

        % Ray triangle intersection with Möller–Trumbore algorithm
        function p = intersectionpoint(cell,d,n)
            % INTERSECTIONPOINT Intersection point between the cell and line/vector/ray
            %
            % P = INTERSECTIONPOINT(CELL,D,N,FIRST) calculates intersection points 
            %   between a set of lines (or vectors) D and a cell. The
            %   Möller–Trumbore algorithm used to find multiple
            %   intersection points.
            %   
            %   If D does not intersect CELL, the coordinates of P are NaN.
            % 
            % See also CSVCell, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            Check.isinteger('N must be either 1 or 2',n,'>=',1,'<=',2)

            if isa(d,'SLine')
                ln_set = d;
            else
                ln_set = d.toline();                
            end

            % Matrices for calculated intersection points
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

            % Set the line set to unit scale
            ordx = floor(log10(abs(ln_set.p1.X)));
            ordy = floor(log10(abs(ln_set.p1.Y)));
            ordz = floor(log10(abs(ln_set.p1.Z)));
            if max([ordx,ordy,ordz]) <= -6
                ln_set.p1 = ln_set.p1 * 10^6;
                ln_set.p2 = ln_set.p2 * 10^6;
            end
            
            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];

            % Calculate R using the following way if the cell will have a
            % single reference vector (along +z) denoting its orientation.

            % angles = orientation(cell,false);
            % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
            % 
            % % Set of rotation matrices
            % Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
            % Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
            % Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];
            % 
            % R = Rz * Ry * Rx;

            % Adjust the line so that it intersects with the cell
            % that is centered at the origin and with a default
            % orientation (Because this is what the mathematical expressions denote)
            ln_set = ln_set.translate(-1*cell.c*(10^-cell.mag));

            for i = 1:numel(ln_set)
                p1 = R' * [ln_set.p1.X(i); ln_set.p1.Y(i); ln_set.p1.Z(i)];
                p2 = R' * [ln_set.p2.X(i); ln_set.p2.Y(i); ln_set.p2.Z(i)];
                p1X(i) = p1(1); p1Y(i) = p1(2); p1Z(i) = p1(3);
                p2X(i) = p2(1); p2Y(i) = p2(2); p2Z(i) = p2(3);
            end

            ln_set = SLine(Point(p1X,p1Y,p1Z),Point(p2X,p2Y,p2Z));
            lnc_set = ln_set.p2 - ln_set.p1;
           
            vertices = [cell.cart.x cell.cart.y cell.cart.z];
            faces = cell.ashp.boundaryFacets;
            vert1 = vertices(faces(:,1),:);
            vert2 = vertices(faces(:,2),:);
            vert3 = vertices(faces(:,3),:);

            for i = 1:numel(ln_set)

                ln = SLine(Point(ln_set.p1.X(i),ln_set.p1.Y(i),ln_set.p1.Z(i)),Point(ln_set.p2.X(i),ln_set.p2.Y(i),ln_set.p2.Z(i)));
                lnc = Point(lnc_set.X(i),lnc_set.Y(i),lnc_set.Z(i)).normalize;
                
                orig = [ln.p1.X ln.p1.Y ln.p1.Z];
                dir = [lnc.X lnc.Y lnc.Z];

                [intersect,~,~,~,xcoor] = TriangleRayIntersection(orig, dir, ...
                    vert1, vert2, vert3, 'lineType' , 'ray', 'border', 'inclusive');
                
                int_x = xcoor(intersect,1);
                int_y = xcoor(intersect,2);
                int_z = xcoor(intersect,3);
                
                t = zeros(size(int_x,1),1);

                for j = 1:size(int_x,1)
                    if round(lnc.X,10) ~= 0
                        tx = (int_x(j)-ln.p1.X)/lnc.X;
                    else
                        tx = NaN;
                    end
                    if round(lnc.Y,10) ~= 0
                        ty = (int_y(j)-ln.p1.Y)/lnc.Y;
                    else
                        ty = NaN;
                    end
                    if round(lnc.Z,10) ~= 0
                        tz = (int_z(j)-ln.p1.Z)/lnc.Z;
                    else
                        tz = NaN;
                    end
                    ts = [tx ty tz];
                    t(j) = mean(ts(~isnan(ts)));
                end

                try 
                    if n == 1
                        [~, minIndex] = min(t);
                        int_point = Point(int_x(minIndex),int_y(minIndex),int_z(minIndex));
                    else
                        [~, maxIndex] = max(t);
                        int_point = Point(int_x(maxIndex),int_y(maxIndex),int_z(maxIndex));
                    end
                    pX(i) = int_point.X;
                    pY(i) = int_point.Y;
                    pZ(i) = int_point.Z;

                catch
                    pX(i) = NaN;
                    pY(i) = NaN;
                    pZ(i) = NaN;
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
       
                int_point = Point(pX,pY,pZ);

                % Use if one reference axis is used to denote the shape
                % int_point = int_point.zrotation(-z_angle);
                % int_point = int_point.yrotation(-y_angle);
                % int_point = int_point.xrotation(-x_angle);

                int_point = int_point.translate(cell.c*(10^-cell.mag));

                % Scale back to original size
                int_point.X = int_point.X * 10^cell.mag;
                int_point.Y = int_point.Y * 10^cell.mag;
                int_point.Z = int_point.Z * 10^cell.mag;

                p = int_point; 
        end

        function ln = perpline(cell,p_set)
            % PERPLINE Line perpendicular to cell passing by point set
            % 
            % LN = PERPLINE(CELL,P-SET) calculates the line set LN perpendicular 
            %   to the cell CELL and passing by the point set P.
            % 
            % See also RBC, Point, SLine.

            Check.isa('P must be a Point',p_set,'Point')
                
            % Matrices for calculated point set
            pX = NaN(size(p_set,1),size(p_set,2));
            pY = NaN(size(p_set,1),size(p_set,2));
            pZ = NaN(size(p_set,1),size(p_set,2));

            % Calculate R using the following way if the cell will have a
            % single reference vector (along +z) denoting its orientation.

            % angles = orientation(cell,false);
            % x_angle = angles(1); y_angle = angles(2); z_angle = angles(3); 
            % 
            % % Set of rotation matrices
            % Rx = [[1, 0, 0]; [0, cos(x_angle), -sin(x_angle)]; [0, sin(x_angle), cos(x_angle)]];
            % Ry = [[cos(y_angle), 0, sin(y_angle)]; [0, 1, 0]; [-sin(y_angle), 0, cos(y_angle)]];
            % Rz = [[cos(z_angle), -sin(z_angle), 0]; [sin(z_angle), cos(z_angle), 0]; [0, 0, 1]];
            % 
            % R = Rz * Ry * Rx;

            R = [cell.or1.Vx cell.or2.Vx cell.or3.Vx;
                 cell.or1.Vy cell.or2.Vy cell.or3.Vy;
                 cell.or1.Vz cell.or2.Vz cell.or3.Vz];

            % Scale point set to unit sizes
            p_set = p_set * (10^-cell.mag);

            % Adjust the point set so that it intersects with the cell
            % that is centered at the origin and with a default
            % orientation (because this is what the mathematical expression denotes).
            p_set = p_set.translate(-1*cell.c);

            for i = 1:numel(p_set)
                vertex = [p_set.X(i); p_set.Y(i); p_set.Z(i)];
                rot_vertex = R' * vertex;
                pX(i) = rot_vertex(1);
                pY(i) = rot_vertex(2);
                pZ(i) = rot_vertex(3);
            end

            p_set = Point(pX,pY,pZ);
            
            % Preset matrices for calculated perpendicular vectors
            ln_setX = NaN(size(p_set,1),size(p_set,2));
            ln_setY = NaN(size(p_set,1),size(p_set,2));
            ln_setZ = NaN(size(p_set,1),size(p_set,2));
            ln_setVX = NaN(size(p_set,1),size(p_set,2));
            ln_setVY = NaN(size(p_set,1),size(p_set,2));
            ln_setVZ = NaN(size(p_set,1),size(p_set,2));

            for i = 1:numel(p_set)
                
                vertices = [cell.cart.x cell.cart.y cell.cart.z];
                faces = cell.ashp.boundaryFacets;
                all = vertices(faces,:);
                
                verticesx = vertices(:,1);
                verticesy = vertices(:,2);
                verticesz = vertices(:,3);

                allx = verticesx(faces,:);
                ally = verticesy(faces,:);
                allz = verticesz(faces,:);

                Xs = reshape(allx, size(faces,1), size(faces,2));
                Ys = reshape(ally, size(faces,1), size(faces,2));
                Zs = reshape(allz, size(faces,1), size(faces,2));

                % Xs = all(1:size(faces,1),:);
                % Ys = all(size(faces,1)+1:2*size(faces,1),:);
                % Zs = all(2*size(faces,1)+1:end,:);

                PXs = repmat([p_set.X(i) p_set.X(i) p_set.X(i)],size(faces,1),1);
                PYs = repmat([p_set.Y(i) p_set.Y(i) p_set.Y(i)],size(faces,1),1);
                PZs = repmat([p_set.Z(i) p_set.Z(i) p_set.Z(i)],size(faces,1),1);

                [~, minIndex] = min(sum(((PXs-Xs).^2+(PYs-Ys).^2+(PZs-Zs).^2).^(1/2),2));
                
                % vert1 = vertices(faces(minIndex,:),1)';
                % vert2 = vertices(faces(minIndex,:),2)';
                % vert3 = vertices(faces(minIndex,:),3)';
                vert1 = vertices(faces(minIndex,1),:);
                vert2 = vertices(faces(minIndex,2),:);
                vert3 = vertices(faces(minIndex,3),:);

                % pl = Plane(Point(vertices(faces(minIndex,1),1),vertices(faces(minIndex,1),2),vertices(faces(minIndex,1),3)), ...
                %     Point(vertices(faces(minIndex,2),1),vertices(faces(minIndex,2),2),vertices(faces(minIndex,2),3)), ...
                %     Point(vertices(faces(minIndex,3),1),vertices(faces(minIndex,3),2),vertices(faces(minIndex,3),3)));
                % 
                % plot(pl.p1)
                % plot(pl.p2)
                % plot(pl.p3)

                faceNormal = cross(vert2-vert1, vert3-vert1,2);

                %ln = perpline(pl,Point(p_set.X(i),p_set.Y(i),p_set.Z(i)));
                
                ln_setX(i) = p_set.X(i); ln_setY(i) = p_set.Y(i); ln_setZ(i) = p_set.Z(i);
                %ln_setVX(i) = ln.p2.X-ln.p1.X; ln_setVY(i) = ln.p2.Y-ln.p1.Y; ln_setVZ(i) = ln.p2.Z-ln.p1.Z;
                ln_setVX(i) = faceNormal(1); ln_setVY(i) = faceNormal(2); ln_setVZ(i) = faceNormal(3);
            end

            ln_set = Vector(ln_setX,ln_setY,ln_setZ,ln_setVX,ln_setVY,ln_setVZ).normalize;
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
            ln_set = ln_set.translate(cell.c);

            % Scale back to original size
            ln_set.p1.X = ln_set.p1.X * (10^cell.mag);
            ln_set.p1.Y = ln_set.p1.Y * (10^cell.mag);
            ln_set.p1.Z = ln_set.p1.Z * (10^cell.mag);
            ln_set.p2.X = ln_set.p2.X * (10^cell.mag);
            ln_set.p2.Y = ln_set.p2.Y * (10^cell.mag);
            ln_set.p2.Z = ln_set.p2.Z * (10^cell.mag);
    
            % Perpendicular line
            ln = ln_set;
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

        function [t_lower, t_upper] = t_limits(cell,ln)
            % TLIMITS lower and upper limits denote the limits of the variable of line
            %  equtions t that intersects with the cell.
            % 
            % [TLOWER, TUPPER] = TLIMITS(CELL,LN) calculates lower and
            %  upper limit of t values of the line LN. This interval of t
            %  values denotes the portion of line LN that possibly intersets
            %  with the cell CELL.
            %
            % See also Cell.

            lnc = ln.p2 - ln.p1;
            if lnc.X > 0
                tx_lower = (min(cell.cart.x) - ln.p1.X) / lnc.X;
                if tx_lower < 0 || round(lnc.X,6) == 0
                    tx_lower = 0;
                end
                tx_upper = (max(cell.cart.x) - ln.p1.X) / lnc.X;
                if tx_upper < 0
                    tx_upper = Inf;
                end
            else
                tx_lower = (max(cell.cart.x) - ln.p1.X) / lnc.X;
                if tx_lower < 0 || round(lnc.X,6) == 0
                    tx_lower = 0;
                end
                tx_upper = (min(cell.cart.x) - ln.p1.X) / lnc.X;
                if tx_upper < 0
                    tx_upper = Inf;
                end
            end
            
            if lnc.Y > 0
                ty_lower = (min(cell.cart.y) - ln.p1.Y) / lnc.Y;
                if ty_lower < 0 || round(lnc.Y,6) == 0
                    ty_lower = 0;
                end
                ty_upper = (max(cell.cart.y) - ln.p1.Y) / lnc.Y;
                if ty_upper < 0
                    ty_upper = Inf;
                end
            else
                ty_lower = (max(cell.cart.y) - ln.p1.Y) / lnc.Y;
                if ty_lower < 0 || round(lnc.Y,6) == 0
                    ty_lower = 0;
                end
                ty_upper = (min(cell.cart.y) - ln.p1.Y) / lnc.Y;
                if ty_upper < 0
                    ty_upper = Inf;
                end
            end
            
            if lnc.Z > 0
                tz_lower = (min(cell.cart.z) - ln.p1.Z) / lnc.Z;
                if tz_lower < 0 || round(lnc.Z,6) == 0
                    tz_lower = 0;
                end
                tz_upper = (max(cell.cart.z) - ln.p1.Z) / lnc.Z;
                if tz_upper < 0
                    tz_upper = Inf;
                end
            else
                tz_lower = (max(cell.cart.z) - ln.p1.Z) / lnc.Z;
                if tz_lower < 0 || round(lnc.Z,6) == 0
                    tz_lower = 0;
                end
                tz_upper = (min(cell.cart.z) - ln.p1.Z) / lnc.Z;
                if tz_upper < 0
                    tz_upper = Inf;
                end
            end
            t_lower = max([tx_lower,ty_lower,tz_lower]);
            t_upper = min([tx_upper,ty_upper,tz_upper]);
        end
        function vol = Volume(cell)
            [~, vol] = convhull(cell.cart.x,cell.cart.y,cell.cart.z);
        end
    end
end
