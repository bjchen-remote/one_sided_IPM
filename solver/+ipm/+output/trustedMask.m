function [trusted,checks] = trustedMask(source,config)
%IPM.OUTPUT.TRUSTEDMASK Apply the shared numerical-quality acceptance contract.
%   MASK = IPM.OUTPUT.TRUSTEDMASK(RESULT) evaluates every recorded state in a
%   version-2 solver result. MASK = IPM.OUTPUT.TRUSTEDMASK(HISTORY,CONFIG) accepts
%   grouped version-2 history and configuration structures. CHECKS exposes
%   the individual logical conditions.

if nargin < 2
    result = ipm.output.validate(source);
    history = result.history;
    config = result.config;
else
    history = source;
end
if ~isstruct(history) || ~isscalar(history) || ...
        ~isstruct(config) || ~isscalar(config)
    error('ipm:TrustedMaskInput', ...
        'History and configuration must be scalar structures.');
end
require_fields(history,{'common','mesh'},'history');
require_fields(config,{'diagnostics','remesh'},'configuration');
values = history.common;
mesh = history.mesh;
diagnostics = config.diagnostics;
remesh = config.remesh;

requiredValues = {'wallPeakLocalMaxima','wallPeakTVRatio', ...
    'positiveWallNegativeRatio','physicalRangeViolation'};
require_fields(values,requiredValues,'history.common');
require_fields(mesh,{'maximumCellRatioX','maximumCellRatioY'}, ...
    'history.mesh');
require_fields(diagnostics,{'oscillationTVTolerance','rangeStopTolerance', ...
    'positiveWallNegativeTolerance'},'configuration.diagnostics');
require_fields(remesh,{'remeshMaximumCellRatio'},'configuration.remesh');

series = cell(1,numel(requiredValues)+2);
for index = 1:numel(requiredValues)
    series{index} = numeric_series(values.(requiredValues{index}), ...
        requiredValues{index});
end
series{end-1} = numeric_series(mesh.maximumCellRatioX, ...
    'maximumCellRatioX');
series{end} = numeric_series(mesh.maximumCellRatioY, ...
    'maximumCellRatioY');
count = numel(series{1});
if any(cellfun(@numel,series) ~= count)
    error('ipm:TrustedMaskLength', ...
        'All quality series must have the same record count.');
end
for index = 1:numel(requiredValues)
    values.(requiredValues{index}) = series{index};
end
mesh.maximumCellRatioX = series{end-1};
mesh.maximumCellRatioY = series{end};
if isfield(mesh,'safetyFactor')
    safetyFactor = numeric_series(mesh.safetyFactor,'safetyFactor');
    if numel(safetyFactor) ~= count
        error('ipm:TrustedMaskLength', ...
            'All quality series must have the same record count.');
    end
else
    % A physical run with no positive wall feature has no scale-resolution
    % gauge. Resolution is then not an acceptance condition.
    safetyFactor = zeros(count,1);
end

finiteConfig = {'oscillationTVTolerance','rangeStopTolerance'};
for index = 1:numel(finiteConfig)
    name = finiteConfig{index};
    validateattributes(diagnostics.(name),{'numeric'}, ...
        {'scalar','real','finite','nonnegative'},mfilename,name);
end
validateattributes(remesh.remeshMaximumCellRatio,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename, ...
    'remeshMaximumCellRatio');
validateattributes(diagnostics.positiveWallNegativeTolerance,{'numeric'}, ...
    {'scalar','real','nonnan','nonnegative'},mfilename, ...
    'positiveWallNegativeTolerance');

checks = struct();
checks.resolvedCfl = safetyFactor <= 1;
checks.densityRange = values.physicalRangeViolation < ...
    diagnostics.rangeStopTolerance;
checks.smoothWallPeak = values.wallPeakLocalMaxima <= 1 | ...
    values.wallPeakTVRatio <= diagnostics.oscillationTVTolerance;
checks.positiveWall = values.positiveWallNegativeRatio <= ...
    diagnostics.positiveWallNegativeTolerance;
cellRatioLimit = remesh.remeshMaximumCellRatio*(1+1e-10);
checks.meshX = mesh.maximumCellRatioX <= cellRatioLimit;
checks.meshY = mesh.maximumCellRatioY <= cellRatioLimit;

trusted = checks.resolvedCfl & checks.densityRange & ...
    checks.smoothWallPeak & checks.positiveWall & ...
    checks.meshX & checks.meshY;
end

function values = numeric_series(values,name)
if ~isnumeric(values) || ~isreal(values) || ~isvector(values) || ...
        isempty(values) || any(~isfinite(values))
    error('ipm:TrustedMaskSeries', ...
        'Quality field "%s" must be a nonempty finite numeric vector.',name);
end
values = values(:);
end

function require_fields(value,names,context)
for index = 1:numel(names)
    if ~isfield(value,names{index})
        error('ipm:TrustedMaskField', ...
            '%s is missing required field "%s".',context,names{index});
    end
end
end
