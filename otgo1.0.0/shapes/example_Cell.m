% Series of examples to demonstrate the use of Cell.
%
% See also Shape, Cell.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Emir Erdem
%   Version: 1.0.0
%   Date: 2023/08/29

example('Use of Cell')

%% DEFINITION OF CELL
exampletitle('DEFINITION OF CELL')

examplecode('p = 32;')
examplecode('shp = ''dumbbell'';')
examplecode('c = Point(randn*1e-6,randn*1e-6,randn*1e-6);')
examplecode('or = Vector(c.X,c.Y,c.Z,rand,rand,rand);')
examplecode('mag = -6;')

examplecode('cell = Cell(p,shp,c,or,mag);')
examplewait()

%% PLOTTING OF CELL
exampletitle('PLOTTING OF CELL')

figure
title('Cell')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cell.plot(''range'',32);')
examplecode('cell.c.plot(''marker'',''.'',''color'',''k'');')
examplecode('Vector(cell.c.X,cell.c.Y,cell.c.Z,cell.or.Vx*1e-6,cell.or.Vy*1e-6,cell.or.Vz*1e-6).plot(''color'',''r'');')
examplewait()

%% VECTORS INTERSECTING CELL
exampletitle('VECTORS INTERSECTING CELL')

figure
title('Vectors intersecting Cell')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cell = Cell(32,''dumbbell'',Point(0,0,0),Vector(0,0,0,0,0,1),-6);')
examplecode('cell.plot();')

examplecode('v = Vector(randn(3,2)*1e-6,randn(3,2)*1e-6,randn(3,2)*1e-6,1*ones(3,2)*1e-6,0*ones(3,2)*1e-6,0*ones(3,2)*1e-6);')
examplecode('v.plot();')
examplecode('v.toline().plot(''range'',[-5*1e-6 5*1e-6]);')

examplecode('p1 = cell.intersectionpoint(v,1,true);')
examplecode('p1.plot(''marker'',''.'',''markersize'',20,''markerfacecolor'',''k'');')

examplecode('p2 = cell.intersectionpoint(v,2,false);')
examplecode('p2.plot(''marker'',''.'',''markersize'',20,''markerfacecolor'',''k'');')
examplewait()

%% SLINES PERPENDICULAR TO CELL & PLANES TANGENT TO CELL
exampletitle('SLINES PERPENDICULAR TO CELL & PLANES TANGENT TO CELL')

figure
title('SLines perpendicular to Cell & Planes tangent to Cell')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cell.plot();')

examplecode('p1.plot();')

examplecode('ln = cell.perpline(p1)')
examplecode('ln.plot();')

examplecode('pl = cell.tangentplane(p1)')
examplecode('pl.plot();')