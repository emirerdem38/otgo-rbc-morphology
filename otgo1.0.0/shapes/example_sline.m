% Series of examples to demonstrate the use of SLine.
%
% See also Shape, SLine.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01


example('Use of SLine')

%% DEFINITION OF SLINES
exampletitle('DEFINITION OF SLINES')

examplecode('p0 = Point(0,0,0);')
examplecode('px = Point(1,0,0);')
examplecode('ln1 = SLine(p0,px)')
examplewait()

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('ln2 = SLine(p,p+px)')
examplewait()

%% PLOTTING OF SLINES
exampletitle('PLOTTING OF SLINES')

figure
title('SLINES')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('ln1.plot(''color'',''k'');')
examplecode('ln2.plot(''color'',''r'');')
examplecode('ln2.p1.plot(''marker'',''x'',''color'',''b'');')
examplecode('ln2.p2.plot(''marker'',''.'',''color'',''b'');')