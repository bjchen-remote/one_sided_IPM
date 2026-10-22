function policy = initialMeshObservationPolicy(input,choices)
%IPM.CONFIG.INITIALMESHOBSERVATIONPOLICY Fixed optional t0 analytic observer.
id='ipm:BadInitialMeshObservationPolicy';
assert(isstruct(input)&&isscalar(input),id,'Initial observation policy must be scalar.');
policy=struct('version',1,'enabled',true,'rule','primitive_k8_fixed_probe_v1', ...
    'innerHalfWidth',4,'innerSpacing',.025,'maximumTailRatio',1.05, ...
    'tailLogSlopeRampCells',20,'maximumObservationNodes',600000,'maximumFallbackAttempts',1);
names=fieldnames(input);assert(isempty(setdiff(names,fieldnames(policy))),id,'Unknown initial observation policy field.');
for k=1:numel(names)
    name=names{k};value=input.(name);expected=policy.(name);
    if strcmp(name,'enabled')
        assert((islogical(value)||isnumeric(value))&&isreal(value)&&isscalar(value)&& ...
            isfinite(value)&&any(value==[0,1]),id,'enabled must be logical or numeric zero/one.');
        policy.enabled=logical(value);
    elseif ischar(expected)
        if isstring(value)&&isscalar(value)&&~ismissing(value),value=char(value);end
        assert(ischar(value)&&isrow(value)&&strcmpi(value,expected),id,'Unknown fixed analytic observation rule.');
    else
        assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&& ...
            double(value)==expected,id,'Initial observation controls are fixed by their rule version.');
    end
end
if nargin<2 || ~policy.enabled,return;end
required={'autonomousMesh','initialCondition','degeneratePower','symmetryMode','transportAnchorX','xlim','ymax'};
assert(isstruct(choices)&&isscalar(choices)&&all(isfield(choices,required)),id,'Initial observer requires complete compatible choices.');
p=choices.autonomousMesh;
assert(p.enabled&&any(p.version == [2,3,4,5]),id, ...
    'The fixed initial observer requires enabled autonomous version two through five.');
assert(ischar(choices.initialCondition)&&strcmp(choices.initialCondition,'degenerate_primitive')&& ...
    choices.degeneratePower==8&&strcmp(choices.symmetryMode,'double_odd_omega')&& ...
    choices.transportAnchorX==1,id,'This rule supports only the original k8 primitive and anchor one.');
assert(choices.xlim(2)>4&&choices.ymax>4&&isequal(choices.xlim,[-choices.xlim(2),choices.xlim(2)]), ...
    id,'The fixed tail construction requires a symmetric box with H>4 and Ymax>4.');
end
