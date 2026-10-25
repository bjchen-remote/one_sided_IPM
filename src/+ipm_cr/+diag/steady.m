function [tracker,record]=steady(s,ctx,policy,tracker)
%STEADY Complete native residuals and nonoverlapping candidate windows.
% Both active overlap copies are retained. No interpolation, shift, amplitude
% fit or exterior exclusion is used. A candidate is not an error certificate.
if nargin<3||isempty(policy),policy=struct();end
policy=resolve(policy);
assert(isstruct(s)&&isscalar(s)&&all(isfield(s,{'W','rhs','scale','step'}))&& ...
    isstruct(s.scale)&&isfield(s.scale,'canonicalTime'),'ipm_cr:SteadyState','Use a complete accepted state.');
if isa(ctx,'ipm_cr.context')
    assert(isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id)&&isfield(s,'ctxKey')&& ...
        isequaln(s.ctxKey,ctx.key),'ipm_cr:SteadyContext','State and prepared context must agree.');
    g=ctx.mesh;key=ctx.key;
else
    % Value meshes permit synthetic logic tests without constructing a PDE.
    assert(isstruct(ctx)&&isscalar(ctx)&&isfield(ctx,'mesh'), ...
        'ipm_cr:SteadyContext','Use a prepared context or an explicit value mesh.');
    g=ctx.mesh;key=[];
end
[masks,radius]=regions(g,policy);
if isempty(key)
    if isfield(g,'geometryKey'),key=g.geometryKey;
    else
        key=struct('coreX',g.core.X,'coreY',g.core.Y,'coreActive',g.core.active, ...
            'polarX',g.polar.X,'polarY',g.polar.Y,'polarActive',g.polar.active);
    end
end
assert(isnumeric(s.step)&&isreal(s.step)&&isscalar(s.step)&&isfinite(s.step)&& ...
    s.step>=0&&s.step==fix(s.step),'ipm_cr:SteadySequence','Accepted step must be a finite nonnegative integer.');
tau=s.scale.canonicalTime;
assert(isnumeric(tau)&&isreal(tau)&&isscalar(tau)&&isfinite(tau), ...
    'ipm_cr:SteadyFinite','Canonical time must be finite.');
data=snapshot(s,g);norms=struct();current=true;
for item={'near','body','global'}
    region=item{1};v=struct('radiusAtMost',radius.(region), ...
        'sampleCount',nnz(masks.(region).core)+nnz(masks.(region).polar), ...
        'duplicateActiveOverlapRetained',true);
    for name=fields()
        f=name{1};scale=f;if strcmp(f,'F'),scale='W';elseif strcmp(f,'Q'),scale='q';end
        v.(f)=metric(data.(f),{data.(scale)},g,masks.(region));
        v.(f).normalizationField=scale;
    end
    current=current&&v.F.mixedInfinity<=policy.residualTol&&v.Q.mixedInfinity<=policy.residualTol;
    norms.(region)=v;
end
record=struct('step',s.step,'tau',tau,'duplicate',false,'regions',norms, ...
    'rates',struct('names',{{'c_l','c_omega','c_r','delta','Uref'}},'values',data.rates), ...
    'currentResidualPassed',current,'window',empty_window(), ...
    'consecutive',0,'completedWindows',0,'candidate',false, ...
    'candidateIsAccuracyCertificate',false,'steadyStateCertified',false, ...
    'trajectoryAccuracyCertified',false,'gridConvergenceCertified',false, ...
    'scope','All native active nodes, including both overlap copies; no alignment or exterior mask.');
nowMax=residuals(norms);
if nargin<4||isempty(tracker)
    tracker=struct('kind','ipm_cr_native_steady_tracker_v1','policy',policy,'key',key, ...
        'anchor',data,'anchorStep',s.step,'anchorTau',tau,'windowResidualPassed',current, ...
        'residualMaxima',nowMax,'consecutive',0,'completedWindows',0,'observations',1, ...
        'last',data,'lastStep',s.step,'lastTau',tau,'lastRecord',record);
    return;
end
assert(isstruct(tracker)&&isscalar(tracker)&&all(isfield(tracker, ...
    {'kind','policy','key','anchor','anchorStep','anchorTau','windowResidualPassed', ...
    'residualMaxima','consecutive','completedWindows','observations','last','lastStep','lastTau','lastRecord'}))&& ...
    strcmp(tracker.kind,'ipm_cr_native_steady_tracker_v1'), ...
    'ipm_cr:SteadyTracker','Use an intact tracker returned by this function.');
assert(isequaln(tracker.policy,policy),'ipm_cr:SteadyPolicyChanged', ...
    'Reset the tracker explicitly before changing candidate policy.');
assert(isequaln(tracker.key,key),'ipm_cr:SteadySetupChanged', ...
    'Reset the tracker explicitly before changing the fixed mathematical setup.');
check_snapshot(tracker.anchor,g);check_snapshot(tracker.last,g);
values=[tracker.anchorStep,tracker.anchorTau,tracker.lastStep,tracker.lastTau, ...
    tracker.consecutive,tracker.completedWindows,tracker.observations];
assert(isnumeric(values)&&isreal(values)&&all(isfinite(values))&& ...
    all([tracker.consecutive,tracker.completedWindows,tracker.observations]>=0)&& ...
    all([tracker.consecutive,tracker.completedWindows,tracker.observations]== ...
    fix([tracker.consecutive,tracker.completedWindows,tracker.observations])), ...
    'ipm_cr:SteadyTracker','Tracker counters and times must remain finite.');
if s.step==tracker.lastStep
    assert(tau==tracker.lastTau&&isequaln(data,tracker.last),'ipm_cr:SteadySequence', ...
        'A repeated accepted step must contain identical time and numerical fields.');
    record=tracker.lastRecord;record.duplicate=true;return;
end
assert(s.step>tracker.lastStep&&tau>tracker.lastTau&&tau>=tracker.anchorTau, ...
    'ipm_cr:SteadySequence','Observe strictly increasing accepted steps and times.');
tracker.windowResidualPassed=tracker.windowResidualPassed&&current;
tracker.residualMaxima=merge_maxima(tracker.residualMaxima,nowMax);
if ~current,tracker.consecutive=0;end
elapsed=tau-tracker.anchorTau;
if elapsed>=policy.windowTau
    changes=struct();maximum=0;
    for name=fields()
        f=name{1};a=data.(f);b=tracker.anchor.(f);
        d=struct('core',a.core-b.core,'polar',a.polar-b.polar);
        changes.(f)=metric(d,{a,b},g,masks.global);
        maximum=max(maximum,changes.(f).mixedInfinity);
    end
    a=data.rates;b=tracker.anchor.rates;
    difference=a-b;denominator=max([ones(size(a));abs(a);abs(b)],[],1);
    rateChange=abs(difference)./denominator;
    assert(all(isfinite(difference))&&all(isfinite(denominator))&&all(isfinite(rateChange)), ...
        'ipm_cr:SteadyFinite','Physical rate differences must remain finite.');
    changes.rates=struct('signedDifference',difference,'denominator',denominator, ...
        'componentMixedDifference',rateChange,'mixedInfinity',max(rateChange));
    maximum=max(maximum,changes.rates.mixedInfinity);
    passed=tracker.windowResidualPassed&&current&&maximum<=policy.changeTol;
    if passed,tracker.consecutive=tracker.consecutive+1;else,tracker.consecutive=0;end
    tracker.completedWindows=tracker.completedWindows+1;
    record.window=struct('performed',true,'fromStep',tracker.anchorStep, ...
        'fromTau',tracker.anchorTau,'toTau',tau,'deltaTau',elapsed, ...
        'changes',changes,'maximumMixedChange',maximum, ...
        'residualMaxima',tracker.residualMaxima,'residualPassed',tracker.windowResidualPassed, ...
        'passed',passed,'scope','All active fixed-grid endpoints; residual maxima include every unique observation since the anchor.');
    record.candidate=passed&&current&&tau>=policy.minTau&&tracker.consecutive>=policy.windows;
    tracker.anchor=data;tracker.anchorStep=s.step;tracker.anchorTau=tau;
    tracker.windowResidualPassed=current;tracker.residualMaxima=nowMax;
end
record.consecutive=tracker.consecutive;record.completedWindows=tracker.completedWindows;
tracker.observations=tracker.observations+1;
tracker.last=data;tracker.lastStep=s.step;tracker.lastTau=tau;tracker.lastRecord=record;
end

function p=resolve(overrides)
p=struct('residualTol',1e-6,'changeTol',1e-5,'windowTau',10,'windows',3, ...
    'minTau',0,'nearRadius',3,'bodyRadius',200);
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:SteadyPolicy','Use a scalar candidate policy.');
names=fieldnames(overrides);
assert(all(ismember(names,fieldnames(p))),'ipm_cr:SteadyPolicy', ...
    'Only residualTol, changeTol, windowTau, windows, minTau, nearRadius and bodyRadius belong to this kernel.');
for k=1:numel(names),p.(names{k})=overrides.(names{k});end
names=fieldnames(p);
for k=1:numel(names)
    v=p.(names{k});assert(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v), ...
        'ipm_cr:SteadyPolicy','Every candidate parameter must be a finite real scalar.');
end
assert(p.residualTol>=0&&p.changeTol>=0&&p.windowTau>0&&p.windows>=1&& ...
    p.windows==fix(p.windows)&&p.minTau>=0&&p.nearRadius>0&&p.bodyRadius>=p.nearRadius, ...
    'ipm_cr:SteadyPolicy','Use nonnegative tolerances/minTau, positive windows/radii and bodyRadius >= nearRadius.');
end

function [masks,radius]=regions(g,p)
assert(isstruct(g)&&isscalar(g)&&all(isfield(g,{'core','polar'})), ...
    'ipm_cr:SteadyMesh','Use the native core and polar value mesh.');
masks=struct('near',struct(),'body',struct(),'global',struct());outer=0;
for patch={'core','polar'}
    key=patch{1};v=g.(key);
    assert(isstruct(v)&&isscalar(v)&&all(isfield(v,{'X','Y','active'}))&&isnumeric(v.X)&&isreal(v.X)&& ...
        isnumeric(v.Y)&&isreal(v.Y)&&isequal(size(v.X),size(v.Y),size(v.active))&& ...
        islogical(v.active)&&all(isfinite(v.X),'all')&&all(isfinite(v.Y),'all'), ...
        'ipm_cr:SteadyMesh','Native coordinates and logical active masks must agree and be finite.');
    r=hypot(v.X,v.Y);
    assert(all(isfinite(r),'all'),'ipm_cr:SteadyMesh','Native radii must remain finite.');
    masks.near.(key)=v.active&(r<=p.nearRadius);
    masks.body.(key)=v.active&(r<=p.bodyRadius);masks.global.(key)=v.active;
    if any(v.active,'all'),outer=max(outer,max(r(v.active)));end
end
for name={'near','body','global'}
    m=masks.(name{1});assert(any(m.core,'all')||any(m.polar,'all'), ...
        'ipm_cr:SteadyMesh','Each candidate region needs at least one active native node.');
end
radius=struct('near',p.nearRadius,'body',p.bodyRadius,'global',outer);
end

function data=snapshot(s,g)
f=s.rhs;
assert(isstruct(f)&&isscalar(f)&&all(isfield(f, ...
    {'q','F','Q','physical','c_l','c_omega','c_r','delta','Uref'})), ...
    'ipm_cr:SteadyState','Use the complete physical RHS, not eta_rhs.');
if isfield(f,'diagnosticsComplete')
    assert(isequal(f.diagnosticsComplete,true),'ipm_cr:SteadyState','Observe a complete accepted endpoint.');
end
rates=[f.c_l,f.c_omega,f.c_r,f.delta,f.Uref];
assert(isnumeric(rates)&&isreal(rates)&&isequal(size(rates),[1,5])&&all(isfinite(rates)), ...
    'ipm_cr:SteadyFinite','Every physical rate must be a finite real scalar.');
assert(isstruct(f.physical)&&all(isfield(f.physical,{'core','polar'})), ...
    'ipm_cr:SteadyField','Physical velocity must cover both patches.');
u1=struct();u2=struct();centered=struct();
for patch={'core','polar'}
    key=patch{1};v=f.physical.(key);
    assert(isstruct(v)&&all(isfield(v,{'u1','u2'})),'ipm_cr:SteadyField','Use both physical velocity components.');
    u1.(key)=v.u1;u2.(key)=v.u2;centered.(key)=v.u1-f.Uref;
end
data=struct('W',s.W,'q',f.q,'F',f.F,'Q',f.Q,'u1',u1,'u1Centered',centered,'u2',u2,'rates',rates);
check_snapshot(data,g);
end

function check_snapshot(data,g)
for name=fields()
    f=name{1};assert(isfield(data,f)&&isstruct(data.(f))&&all(isfield(data.(f),{'core','polar'})), ...
        'ipm_cr:SteadyField','Every observed field must cover both native patches.');
    for patch={'core','polar'}
        key=patch{1};v=data.(f).(key);
        assert(isnumeric(v)&&isreal(v)&&isequal(size(v),size(g.(key).X)), ...
            'ipm_cr:SteadyField','Observed arrays must have native geometry dimensions.');
        assert(all(isfinite(v),'all'),'ipm_cr:SteadyFinite','Every observed native field must be finite.');
    end
end
assert(isfield(data,'rates')&&isnumeric(data.rates)&&isreal(data.rates)&& ...
    isequal(size(data.rates),[1,5])&&all(isfinite(data.rates)), ...
    'ipm_cr:SteadyFinite','Anchor rates must remain finite.');
end

function names=fields()
names={'W','q','F','Q','u1','u1Centered','u2'};
end

function m=metric(field,references,g,masks)
peak=-1;scale=0;location=struct();
for patch={'core','polar'}
    key=patch{1};indices=find(masks.(key));if isempty(indices),continue;end
    values=field.(key)(indices);
    assert(all(isfinite(values)),'ipm_cr:SteadyFinite','Window differences must remain finite.');
    for j=1:numel(references),scale=max(scale,max(abs(references{j}.(key)(indices))));end
    [value,k]=max(abs(values));
    if value>peak
        index=indices(k);p=g.(key);peak=value;
        location=struct('patch',key,'linearIndex',index,'X',p.X(index),'Y',p.Y(index), ...
            'r',hypot(p.X(index),p.Y(index)),'signedValue',field.(key)(index));
    end
end
denominator=max(1,scale);
m=struct('absoluteInfinity',peak,'denominator',denominator, ...
    'mixedInfinity',peak/denominator,'maximumLocation',location);
end

function out=residuals(norms)
out=struct();
for name={'near','body','global'}
    key=name{1};out.(key)=struct('F',norms.(key).F.mixedInfinity,'Q',norms.(key).Q.mixedInfinity);
end
end

function out=merge_maxima(a,b)
out=a;
for name={'near','body','global'}
    key=name{1};
    assert(isfinite(a.(key).F)&&isfinite(a.(key).Q),'ipm_cr:SteadyTracker','Stored residual maxima must remain finite.');
    out.(key).F=max(a.(key).F,b.(key).F);out.(key).Q=max(a.(key).Q,b.(key).Q);
end
end

function w=empty_window()
w=struct('performed',false,'fromStep',[],'fromTau',[],'toTau',[], ...
    'deltaTau',[],'changes',struct(),'maximumMixedChange',[], ...
    'residualMaxima',struct(),'residualPassed',false,'passed',false, ...
    'scope','No new complete nonoverlapping window.');
end
