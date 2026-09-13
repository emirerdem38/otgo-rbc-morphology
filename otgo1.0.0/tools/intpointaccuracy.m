function acc = intpointaccuracy(d, justshp, c, varargin)
            
    % ACC = INTPOINTACCURACY(D,JUSTSHP,VARARGIN) returns how accurate the 
    % intersection points calcualtions are by computing the intersection 
    % points of the given ray set D and on the shape JUSTSHP, using both
    % the numerical algorithm (whose accuracy will be assesed) and a
    % theoretical approach. Then two results are compared to return
    % the accuracy. Currently these shapes are sphere and healthy-steady 
    % RBC. 
    %
    % Using varagin, a radius for the sphere can be specified,
    % meanwhile the RBC will always have default dimensions.
    %
    % Given ray set should be focused at the origin for the best
    % results.

    Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')

    if isa(d,'SLine')
        ln_set = d;
    else
        ln_set = d.toline();                
    end   

    if ~isempty(varargin)
        radius = varargin{1};
    else
        radius = 4*1e-6;
    end

    if strcmp(justshp, 'sphere')
        vol = (4/3)*pi*radius^3; % 10^-18 m^3
        shape_reference = Spherical(c,radius);
        shape_actual = Cell(64, 'sphere', c, Vector(0,0,0,1,0,0), ...
            Vector(0,0,0,0,1,0), Vector(0,0,0,0,0,1), vol*1e+18, -6);
        
    elseif strcmp(justshp, 'RBC')
        shape_reference = RBC(c,Vector(0,0,0,1,0,0),3.91*1e-6,0.81*1e-6,2.52*1e-6,2.76*1e-6);
        shape_actual = Cell(64, 'RBC', c, Vector(0,0,0,1,0,0), ...
            Vector(0,0,0,0,1,0), Vector(0,0,0,0,0,1), 'vol', -6);
    end

    pset_reference_I1 = intersectionpoint(shape_reference,ln_set,1);
    X_reference_I1 = pset_reference_I1.X;
    Y_reference_I1 = pset_reference_I1.Y;
    Z_reference_I1 = pset_reference_I1.Z;

    pset_reference_I2 = intersectionpoint(shape_reference,ln_set,2);
    X_reference_I2 = pset_reference_I2.X;
    Y_reference_I2 = pset_reference_I2.Y;
    Z_reference_I2 = pset_reference_I2.Z;

    pset_actual_I1 = intersectionpoint(shape_actual,ln_set,1);
    X_actual_I1 = pset_actual_I1.X;
    Y_actual_I1 = pset_actual_I1.Y;
    Z_actual_I1 = pset_actual_I1.Z;

    pset_actual_I2 = intersectionpoint(shape_actual,ln_set,2);
    X_actual_I2 = pset_actual_I2.X;
    Y_actual_I2 = pset_actual_I2.Y;
    Z_actual_I2 = pset_actual_I2.Z;

    figure;
    title('Reference Cell')
    hold on;
    plot(shape_reference)
    %plot(ln_set)
    plot(pset_reference_I1)
    plot(pset_reference_I2)

    figure;
    title('Actual Cell')
    hold on;
    plot(shape_actual)
    %plot(ln_set)
    plot(pset_actual_I1)
    plot(pset_actual_I2)

    not_calculated1 = sum(sum(isnan(X_actual_I1)));
    not_calculated2 = sum(sum(isnan(X_actual_I2)));
    fprintf('%d first intersection points were not calculated!\n', not_calculated1);
    fprintf('%d second intersection points were not calculated!\n', not_calculated2);

    nanmask1 = ~isnan(X_actual_I1);
    nanmask2 = ~isnan(X_actual_I2);

    X_remaning_reference1 = X_reference_I1(nanmask1);
    Y_remaning_reference1 = Y_reference_I1(nanmask1);
    Z_remaning_reference1 = Z_reference_I1(nanmask1);
    X_remaning_actual1 = X_actual_I1(nanmask1);
    Y_remaning_actual1 = Y_actual_I1(nanmask1);
    Z_remaning_actual1 = Z_actual_I1(nanmask1);

    dists1 = ((X_remaning_reference1-X_remaning_actual1).^2+ ...
       (Y_remaning_reference1-Y_remaning_actual1).^2+(Z_remaning_reference1-Z_remaning_actual1).^2).^(1/2);

    X_remaning_reference2 = X_reference_I2(nanmask2);
    Y_remaning_reference2 = Y_reference_I2(nanmask2);
    Z_remaning_reference2 = Z_reference_I2(nanmask2);
    X_remaning_actual2 = X_actual_I2(nanmask2);
    Y_remaning_actual2 = Y_actual_I2(nanmask2);
    Z_remaning_actual2 = Z_actual_I2(nanmask2);

    dists2 = ((X_remaning_reference2-X_remaning_actual2).^2+ ...
       (Y_remaning_reference2-Y_remaning_actual2).^2+(Z_remaning_reference2-Z_remaning_actual2).^2).^(1/2);

    disp('*------------------------------------------------------------------------*')
    fprintf('1st intersection analysis was done with %d points\n', numel(ln_set)-not_calculated1)
    fprintf('Mean Distance: %g\n', mean(dists1));
    fprintf('Total Distance: %g\n', sum(dists1));
    fprintf('Max Distance: %g\n', max(dists1));
    fprintf('Root Mean Square Error: %g\n', sqrt(sum(dists1.^2)/numel(dists1)));

    disp('*------------------------------------------------------------------------*')
    fprintf('2nd intersection analysis was done with %d points\n', numel(ln_set)-not_calculated2)
    fprintf('Mean Distance: %g\n', mean(dists2));
    fprintf('Total Distance: %g\n', sum(dists2));
    fprintf('Max Distance: %g\n', max(dists2));
    fprintf('Root Mean Square Error: %g\n', sqrt(sum(dists2.^2)/numel(dists2)));

    acc = [mean(dists1) sum(dists1) max(dists1) sqrt(sum(dists1.^2)/numel(dists1));
           mean(dists2) sum(dists2) max(dists2) sqrt(sum(dists2.^2)/numel(dists2))];

end