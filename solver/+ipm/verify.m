function report = verify(suite,stabilityMode,originalRoot)
%IPM.VERIFY Run a verification suite without permanently changing the path.
%   IPM.VERIFY() runs the maintained baseline checks. Suite choices are
%   'baseline', 'fourth', 'sixth', 'quadrant', 'all', and 'equivalence'. The optional
%   stability mode is 'quick' (default) or 'heavy' for fourth-order checks.
%   Equivalence alone needs the original solver; pass its directory as the
%   third argument or retain the sibling sixth_order_integration directory.

if nargin < 1 || isempty(suite)
    suite = 'baseline';
end
if nargin < 2 || isempty(stabilityMode)
    stabilityMode = 'quick';
end
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 3 || isempty(originalRoot)
    originalRoot = fullfile(fileparts(root),'sixth_order_integration');
end
suite = validatestring(lower(char(string(suite))), ...
    {'baseline','fourth','sixth','quadrant','all','equivalence'},mfilename,'suite',1);
stabilityMode = validatestring(lower(char(string(stabilityMode))), ...
    {'quick','heavy'},mfilename,'stabilityMode',2);

savedPath = path;
pathCleanup = onCleanup(@()path(savedPath));
addpath(fullfile(root,'tests'));
started = tic;
layoutReport = ipmtests.layout();
switch suite
    case 'baseline'
        report = ipmtests.baseline.suite();
    case 'fourth'
        report = ipmtests.fourth.suite(stabilityMode);
    case 'sixth'
        report = ipmtests.sixth.suite();
    case 'quadrant'
        report = ipmtests.quadrant();
    case 'all'
        report = struct();
        report.baseline = ipmtests.baseline.suite();
        report.fourth = ipmtests.fourth.suite(stabilityMode);
        report.sixth = ipmtests.sixth.suite();
        report.quadrant = ipmtests.quadrant();
    case 'equivalence'
        report = ipmtests.equivalence(originalRoot);
end
report.layout = layoutReport;
report.passed = true;
report.elapsedSeconds = toc(started);
end
