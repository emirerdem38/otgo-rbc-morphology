% Series of examples to demonstrate the use of Vector.
%
% See also Shape, Vector.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01

example('Use of Vector')

%% DEFINITION OF POINTS
exampletitle('DEFINITION OF VECTORS')

examplecode('vx = Vector(0,0,0,1,0,0)')
examplecode('vy = Vector(0,0,0,0,1,0)')
examplecode('vz = Vector(0,0,0,0,0,1)')
examplewait()

%% OPERATIONS ON VECTORS
exampletitle('OPERATIONS ON VECTORS')

examplecode('v1 = 2*vx+3*vy')
examplecode('v2 = vy*vx')
examplewait()

%% PLOTTING OF VECTORS
exampletitle('PLOTTING OF VECTORS')


figure
title('VECTORS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('vx.plot(''color'',''k'');')
examplecode('vy.plot(''color'',''k'');')
examplecode('vz.plot(''color'',''k'');')
examplecode('v1.plot(''color'',''r'');')
examplecode('v2.plot(''color'',''g'');')