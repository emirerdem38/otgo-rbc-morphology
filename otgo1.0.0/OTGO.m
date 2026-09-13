% Optical Tweezers in the Geometrical Optics Approximation
% Version 1.0.0
%
% This script loads all software packages.
%
% Optical Tweezers Software packages
%   <a href="matlab:help utility">utility</a>   - (folder) : general utility functions (to be loaded always)
%   <a href="matlab:help tools">tools</a>       - (folder) : common tools
%   <a href="matlab:help shapes">shapes</a>     - (folder) : 3D geometrical shapes
%   <a href="matlab:help beams">beams</a>       - (folder) : optical beams
%   <a href="matlab:help go">go</a>             - (folder) : Geometrical Optics
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01

clc

format long

addpath(cd)

fprintf('\n')
fprintf('Optical Tweezers in Geometrical Optics\n')
fprintf('version 1.0.0\n')
fprintf('\n')

fprintf('loading <a href="matlab:help utility">utility</a> - general utility functions ...\n')
addpath([cd filesep 'utility'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help tools">tools</a> - common tools ...\n')
addpath([cd filesep 'tools'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help shapes">shapes</a> - 3D geometrical shapes ...\n')
addpath([cd filesep 'shapes'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help beams">beams</a> - optical beams ...\n')
addpath([cd filesep 'beams'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help go">go</a> - geometrical optics ...\n')
addpath([cd filesep 'go'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help numeric_shapes">numeric_shapes</a> - Numeric Cells ...\n')
addpath([cd filesep 'numeric_shapes'])
fprintf('loaded\n')
fprintf('\n')