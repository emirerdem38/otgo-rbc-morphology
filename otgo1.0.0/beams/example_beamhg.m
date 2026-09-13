% Series of examples to demonstrate the use of BeamHG.
%
% See also Beam, BeamHG.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01


example('Use of BeamHG')

%% DEFINITION OF BEAMHG
exampletitle('DEFINITION OF BEAMHG')

examplecode('m = 2;')
examplecode('n = 1;')
examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 1i;')
examplecode('R = 10e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bhg = BeamHG(m,n,Ex0,Ey0,w0,R,Nphi,Nr)')
examplewait()

%% PLOTTING OF BEAMHG
exampletitle('PLOTTING OF BEAMHG')

figure
title('BEAMHG')

examplecode('bhg.plot();')