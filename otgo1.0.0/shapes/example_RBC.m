% Series of examples to demonstrate the use of Red Blood Cell (RBC).
%
% See also Shape, RBC.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Emir Erdem
%   Version: 1.0.0
%   Date: 2023/08/23

example('Use of RBC')

%% DEFINITION OF RBCs
exampletitle('DEFINITION OF RBCs')

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(5*randn(m,n)*1e-6,5*randn(m,n)*1e-6,5*randn(m,n)*1e-6);')
examplecode('or = Vector(zeros(m,n),zeros(m,n),zeros(m,n),rand(m,n),rand(m,n),rand(m,n));')
examplecode('r = ones(m,n).*3.91*1e-6;')
examplecode('t_min = ones(m,n).*0.81*1e-6;')
examplecode('t_max = ones(m,n).*2.52*1e-6;')
examplecode('d = ones(m,n).*2.76*1e-6;')

examplecode('rbc = RBC(p,or,r,t_min,t_max,d);')
examplewait()

%% PLOTTING OF RBCs
exampletitle('PLOTTING OF RBCs')

figure
title('Red Blood Cells (RBC)')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('rbc.plot(''range'',32);')
examplecode('rbc.c.plot(''marker'',''.'',''color'',''k'');')
examplecode('Vector(rbc.c.X,rbc.c.Y,rbc.c.Z,rbc.or.Vx*1e-6,rbc.or.Vy*1e-6,rbc.or.Vz*1e-6).plot(''color'',''r'');')
examplewait()

%% VECTORS INTERSECTING RBC
exampletitle('VECTORS INTERSECTING RBC')

figure
title('Vectors intersecting RBC')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('rbc = RBC(Point(0,0,0),Vector(0,0,0,0,0,1),3.91*(1e-6),0.81*(1e-6),2.52*(1e-6),(2.76)*(1e-6));')
examplecode('rbc.plot();')

examplecode('v = Vector(randn(m,n)*1e-6,randn(m,n)*1e-6,randn(m,n)*1e-6,1*ones(m,n)*1e-6,0*ones(m,n)*1e-6,0*ones(m,n)*1e-6);')
examplecode('v.plot();')
examplecode('v.toline().plot(''range'',[-5*1e-6 5*1e-6]);')

examplecode('p1 = rbc.intersectionpoint(v,1);')
examplecode('p1.plot(''marker'',''.'',''markersize'',20,''markerfacecolor'',''k'');')

examplecode('p2 = rbc.intersectionpoint(v,2);')
examplecode('p2.plot(''marker'',''.'',''markersize'',20,''markerfacecolor'',''k'');')
examplewait()

%% SLINES PERPENDICULAR TO RBC & PLANES TANGENT TO RBC
exampletitle('SLINES PERPENDICULAR TO RBC & PLANES TANGENT TO RBC')

figure
title('SLines perpendicular to RBC & Planes tangent to RBC')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('rbc.plot();')

examplecode('p1.plot();')

examplecode('ln = rbc.perpline(p1)')
examplecode('ln.plot();')

examplecode('pl = rbc.tangentplane(p1)')
examplecode('pl.plot();')