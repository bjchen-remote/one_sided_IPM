function c=config(overrides)
%CONFIG Standalone elliptic prototype; physical wall is y=0, psi=0.
c=struct('coreHalfWidth',1.2,'coreHeight',1.2,'coreDegree',[32,32], ...
    'rInner',.75,'rOuter',16,'polarDegree',16,'angularDegree',40, ...
    'logPanelWidth',log(2),'outerBoundary','decay', ...
    'linearTolerance',1e-10,'traceTolerance',1e-10, ...
    'overlapValueTolerance',1e-7,'overlapFluxTolerance',1e-6, ...
    'strongTolerance',1e-5,'requireAccuracy',true,'auditPoints',81);
if nargin>0
    assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:Config','Use scalar options.');
    names=fieldnames(overrides);
    for k=1:numel(names)
        assert(isfield(c,names{k}),'ipm_cr:Config','Unknown option %s.',names{k});
        c.(names{k})=overrides.(names{k});
    end
end
for name={'coreHalfWidth','coreHeight','rInner','rOuter','logPanelWidth', ...
        'linearTolerance','traceTolerance','overlapValueTolerance', ...
        'overlapFluxTolerance','strongTolerance'}
    validateattributes(c.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(c.coreDegree,{'numeric'},{'vector','numel',2,'integer','>=',6,'<=',192});
for name={'polarDegree','angularDegree','auditPoints'}
    validateattributes(c.(name{1}),{'numeric'},{'scalar','integer','>=',6,'<=',256});
end
assert(c.rInner<.95*min(c.coreHalfWidth,c.coreHeight), ...
    'ipm_cr:Overlap','A strict positive overlap is required.');
assert(c.rOuter>hypot(c.coreHalfWidth,c.coreHeight), ...
    'ipm_cr:OuterRadius','The annulus must contain all artificial core edges.');
c.outerBoundary=char(validatestring(c.outerBoundary,{'decay','dirichlet'}));
validateattributes(c.requireAccuracy,{'logical'},{'scalar'});
end
