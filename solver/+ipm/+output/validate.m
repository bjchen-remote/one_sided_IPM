function result = validate(inputResult)
%IPM.OUTPUT.VALIDATE Load and validate one version-2 solver result.
%   RESULT = IPM.OUTPUT.VALIDATE(INPUT) accepts a scalar structure with
%   schemaVersion=2 or a MAT-file whose variable result has that contract.

raw = read_input(inputResult);
if ~isstruct(raw) || ~isscalar(raw)
    error('ipm:ResultScalar','A solver result must be a scalar structure.');
end

if ~isfield(raw,'schemaVersion') || ...
        ~isnumeric(raw.schemaVersion) || ~isreal(raw.schemaVersion) || ...
        ~isscalar(raw.schemaVersion) || ~isfinite(raw.schemaVersion) || ...
        raw.schemaVersion ~= 2
    error('ipm:ResultSchemaVersion', ...
        ['Only results with explicit schemaVersion=2 are supported; ' ...
        'legacy, unversioned, and other result schemas are rejected.']);
end
result = raw;
ipm.output.validateV2(result);
end

function raw = read_input(inputResult)
if ischar(inputResult) || (isstring(inputResult) && isscalar(inputResult))
    fileName = char(inputResult);
    if ~isfile(fileName)
        error('ipm:ResultFileMissing', ...
            'Result MAT-file does not exist: %s',fileName);
    end
    loaded = load(fileName);
    if ~isfield(loaded,'result')
        error('ipm:ResultFileContract', ...
            ['The MAT-file must contain a variable named result with ' ...
            'schemaVersion=2; legacy aliases and inferred variables are ' ...
            'not supported.']);
    end
    raw = loaded.result;
elseif isstruct(inputResult)
    raw = inputResult;
else
    error('ipm:ResultInput', ...
        'Pass a solver result structure or the path to a MAT-file.');
end
end
