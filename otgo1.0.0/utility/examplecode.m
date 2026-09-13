function examplecode(codeline,pausetime)
% EXAMPLECODE   Auxiliary function for examples - Execute code
%
% EXAMPLECODE(CODELINE,PAUSETIME) executes and displays CODELINE 
%   and then waits for PAUSETIME seconds.
%
% EXAMPLECODE(CODELINE) executes and displays CODELINE 
%   and then waits for 1 second.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01


if nargin<2
    pausetime = 1;
end

fprintf(['>> ' codeline '\n'])
evalin('caller',codeline)
pause(pausetime)