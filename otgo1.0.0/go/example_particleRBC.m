% Series of examples to demonstrate the use of ParticleRBC.
%
% See also Particle, ParticleRBC.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Emir Erdem
%   Version: 1.0.0
%   Date: 2023/08/23


example('Use of ParticleRBC')

%% Definition of ParticleRBC
exampletitle('Definition of PARTICLERBC')

examplecode('c = Point(1*1e-6,1*1e-6,1*1e-6);')
examplecode('or = Vector(0,0,0,0,0,1);')
examplecode('r = 3.91*1e-6;')
examplecode('t_min = 0.81*1e-6;')
examplecode('t_max = 2.52*1e-6;')
examplecode('d = 2.76*1e-6;')
examplecode('nm = 1.33;')
examplecode('np = 1.5;')
examplecode('particle = ParticleRBC(c,or,r,t_min,t_max,d,nm,np)')
examplewait()

%% Plotting of ParticleRBC
exampletitle('Plotting of PARTICLERBC')

figure
title('Particle RBC')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('particle.plot();')
examplewait()

%% Scattering
exampletitle('Scattering')

examplecode('mr = 3;')
examplecode('nr = 2;')
examplecode('v = Vector(zeros(mr,nr)*1e-6,zeros(mr,nr)*1e-6,zeros(mr,nr)*1e-6,rand(mr,nr)*1e-6,rand(mr,nr)*1e-6,rand(mr,nr)*1e-6);')
examplecode('P = ones(mr,nr);')
examplecode('pol = Vector(zeros(mr,nr)*1e-6,zeros(mr,nr)*1e-6,zeros(mr,nr)*1e-6,ones(mr,nr)*1e-6,ones(mr,nr)*1e-6,ones(mr,nr)*1e-6); pol = v*pol;')
examplecode('r = Ray(v,P,pol);')
examplewait()

examplecode('r.plot(''color'',''k'');')
examplewait()

examplecode('r_vec = particle.scattering(r)')
examplewait()

examplecode('rr = r_vec(1).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(1).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(2).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(2).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(3).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(3).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(4).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(4).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(5).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(5).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

%% FORCE
exampletitle('FORCE')

examplecode('F = particle.force(r_vec,r) % fN')
examplewait()

%% TORQUE
exampletitle('TORQUE')

examplecode('T = particle.torque(r_vec,r) % fN*nm')