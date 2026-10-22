function report = verify()
%IPM.VERIFY Check the standalone wall-integral gauge and short restart path.
%   Run from this release's root with IPM.VERIFY(). The checks use a reduced
%   fixed mesh and five short steps; no production run is started.

root = fileparts(fileparts(mfilename('fullpath')));
savedPath = path;
pathCleanup = onCleanup(@()path(savedPath)); %#ok<NASGU>
addpath(fullfile(root,'tests'),'-begin');
assert(startsWith(which('ipm.solve'),[root,filesep]), ...
    'Use the standalone release whose +ipm package is on the MATLAB path.');
assert(startsWith(which('check_wall_integral'),[root,filesep]), ...
    'The standalone verification file is missing or shadowed.');
report = check_wall_integral();
end
