function setup_paths()
%SETUP_PATHS Add OTGO 1.0.0 packages to the MATLAB path (relative to this folder).
%
% Run once at the start of each demo:
%   setup_paths;

here = fileparts(mfilename('fullpath'));
otgoRoot = fullfile(here, '..', 'otgo1.0.0');
if ~exist(fullfile(otgoRoot, 'OTGO.m'), 'file')
    error('Could not find otgo1.0.0/OTGO.m next to paper_examples/.');
end

addpath(otgoRoot);
addpath(fullfile(otgoRoot, 'utility'));
addpath(fullfile(otgoRoot, 'tools'));
addpath(fullfile(otgoRoot, 'shapes'));
addpath(fullfile(otgoRoot, 'beams'));
addpath(fullfile(otgoRoot, 'go'));
addpath(fullfile(otgoRoot, 'numeric_shapes'));
end
