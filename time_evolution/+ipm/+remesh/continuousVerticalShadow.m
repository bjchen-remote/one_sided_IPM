function shadow=continuousVerticalShadow(shadow,state,policy,baseline)
%IPM.REMESH.CONTINUOUSVERTICALSHADOW Restartable, bounded Y observer audit.
% The native controller and all safety/transaction gates remain unchanged.
if isempty(shadow)
    shadow=struct('version',1,'originStep',state.step, ...
        'observations',0,'differentRequests',0, ...
        'baselineOnlyRequests',0,'continuousOnlyRequests',0, ...
        'maximumRelativeCellGap',0, ...
        'window',struct('time',[],'step',[],'epoch',[], ...
        'nodalCells',[],'continuousCells',[]), ...
        'latest',struct());
end
t=state.scale.canonicalTime;n=state.step;epoch=state.ops.remeshCount;
q=ipm.remesh.continuousVerticalCore(state.flow,state.ops.x, ...
    state.ops.y,state.config.scaling.adaptiveLevels);
w=shadow.window;
if ~isempty(w.time)
    assert(t>=w.time(end) && n>=w.step(end) && epoch>=w.epoch(end), ...
        'ipm:ContinuousVerticalShadowClock','The accepted state moved backwards.');
    if epoch~=w.epoch(end)
        w.time=[];w.step=[];w.epoch=[];w.nodalCells=[];w.continuousCells=[];
    end
end
w.time(end+1,1)=t;w.step(end+1,1)=n;w.epoch(end+1,1)=epoch;
w.nodalCells(end+1,1)=state.flow.verticalCoreGridPoints;
w.continuousCells(end+1,1)=q.coreCells;
keep=w.time>=t-policy.trendWindow;
names=fieldnames(w);
for k=1:numel(names),w.(names{k})=w.(names{k})(keep,:);end
shadow.window=w;
trend=numel(w.time)>=4;decay=0;
if trend
    fit=polyfit(w.time-w.time(end),log(w.continuousCells),1);
    secants=-diff(log(w.continuousCells))./diff(w.time);
    decay=max([0;-fit(1);secants]);
end
predictedY=q.coreCells*exp(-policy.maximumReviewInterval*decay);
core=[baseline.coreCells(1),q.coreCells];
predicted=[baseline.predictedCoreCells(1),predictedY];
reason='resolved';requested=false;
if any(core<policy.regridCoreTrigger)
    reason='current_core_trigger';requested=true;
elseif trend && any(predicted<policy.predictedCoreBuffer & ...
        core<policy.targetCoreCells)
    reason='forecast_core_trigger';requested=true;
elseif state.flow.safetyFactor>=policy.maximumSafety
    reason='safety_buffer';requested=true;
end
shadow.observations=shadow.observations+1;
different=baseline.requested~=requested;
shadow.differentRequests=shadow.differentRequests+double(different);
shadow.baselineOnlyRequests=shadow.baselineOnlyRequests+ ...
    double(baseline.requested && ~requested);
shadow.continuousOnlyRequests=shadow.continuousOnlyRequests+ ...
    double(requested && ~baseline.requested);
gap=abs(q.coreCells-w.nodalCells(end))/q.coreCells;
shadow.maximumRelativeCellGap=max(shadow.maximumRelativeCellGap,gap);
shadow.latest=struct('step',n,'canonicalTime',t,'remeshCount',epoch, ...
    'nodalCells',w.nodalCells(end),'continuousCells',q.coreCells, ...
    'continuousWidth',q.coreWidth,'continuousPeakX',q.peakX, ...
    'nodalPeakX',state.flow.omegaGaugePeakX, ...
    'baselineRequested',baseline.requested,'baselineReason',baseline.reason, ...
    'continuousRequested',requested,'continuousReason',reason, ...
    'continuousPredictedCells',predicted,'continuousDecayEstimate',decay, ...
    'continuousTrendAvailable',trend, ...
    'baselineSafetyShared',true,'decisionDifferent',different);
end
