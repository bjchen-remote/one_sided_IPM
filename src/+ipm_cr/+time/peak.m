function p=peak(f)
% Observe existing measurements only; missing diagnostics never reject a step.
p=struct('x',NaN,'cellIndex',NaN,'cellFraction',NaN, ...
    'tieNonunique',false,'curvature',NaN,'available',false);
if ~isstruct(f)||~isscalar(f)||~isfield(f,'measurement')|| ...
        ~isstruct(f.measurement)||~isfield(f.measurement,'peak')|| ...
        ~isstruct(f.measurement.peak)||~isscalar(f.measurement.peak)
    return;
end
source=f.measurement.peak;names={'x','cellIndex','cellFraction','tieNonunique','curvature'};
for k=1:numel(names)
    if isfield(source,names{k}),p.(names{k})=source.(names{k});end
end
p.available=true;
end
