function varargout=history(action,varargin)
%HISTORY Accepted-step columns, observations and durable run summaries.
switch action
    case 'columns',varargout{1}=columns();
    case 'row',varargout{1}=history_row(varargin{:});
    case 'summary',varargout{1}=finish_summary(varargin{:});
    case 'write',write_summary(varargin{:});
    otherwise,error('ipm_cr:HistoryAction','Unknown history action %s.',action);
end
end

function names=columns()
names={'step','tau','physicalTime','physicalTimeLow','dt','maximumStageRate','maximumStageCfl', ...
    'configuredCfl','retries','stepSeconds','accumulatedWallSeconds','gaugeResidualInfinity', ...
    'poissonBackwardError','poissonStrongResidual','poissonStrongAvailable','overlapValueMismatch', ...
    'overlapVelocityMismatch','normalFluxMismatch','traceResidual','cl','cw','cr', ...
    'wallI','wallP','wallC','centeredDensityMin','centeredDensityMax','densityOrigin', ...
    'totalDensityMin','totalDensityMax','QInfinity','centeredFInfinity', ...
    'rhsEvaluations','peakEventLimited','peakEventLocated','peakEventProbeCount','peakEventExtraRhs', ...
    'peakTieStages','minimumPeakCurvatureMagnitude'};
end

function row=history_row(s,info,wall)
p=s.rhs.poisson;event=get(info,'peakEvent',struct());
peakTies=0;peakCurvature=NaN;
if isfield(info,'stagePeakDiagnostics')&&~isempty(info.stagePeakDiagnostics)
    peaks=info.stagePeakDiagnostics;peakTies=nnz([peaks.tieNonunique]);
    peakCurvature=min(abs([peaks.curvature]));
end
centeredMin=min([min(s.W.core,[],'all'),min(s.W.polar,[],'all')]);
centeredMax=max([max(s.W.core,[],'all'),max(s.W.polar,[],'all')]);
row=[s.step,s.scale.canonicalTime,s.scale.physicalTime,s.scale.physicalTimeLow, ...
    get(info,'dt',NaN),get(info,'maximumStageRate',NaN),get(info,'realizedMaximumStageCfl',NaN), ...
    get(info,'configuredCfl',NaN),get(info,'retries',0),get(info,'seconds',0),wall, ...
    max(max(abs(s.rhs.gaugeResidual)),get(info,'maximumStageGaugeResidual',0)),get(p,'linearSolveBackwardError',NaN), ...
    get(p,'strongResidual',NaN),get(p,'strongResidualAvailable',false), ...
    get(p,'overlapValueMismatch',NaN),get(p,'overlapVelocityMismatch',NaN), ...
    get(p,'normalFluxMismatch',NaN),get(p,'traceResidual',NaN), ...
    s.rhs.c_l,s.rhs.c_omega,s.rhs.c_r,s.rhs.measurement.I,s.rhs.measurement.P,s.rhs.measurement.C, ...
    centeredMin,centeredMax,s.scale.densityOrigin,centeredMin+s.scale.densityOrigin, ...
    centeredMax+s.scale.densityOrigin, ...
    max([max(abs(s.rhs.Q.core),[],'all'),max(abs(s.rhs.Q.polar),[],'all')]), ...
    max([max(abs(s.rhs.F.core),[],'all'),max(abs(s.rhs.F.polar),[],'all')]), ...
    get(info,'rhsEvaluations',0),get(event,'limited',false),get(event,'located',false), ...
    get(event,'probeCount',0),get(event,'extraRhsEvaluations',0),peakTies,peakCurvature];
end

function s=finish_summary(s,status,state,g,method,history,names,pointer,priorWall,elapsed)
s.status=status;s.endedUtc=ipm_cr.io.identity('utc');s.tau=state.scale.canonicalTime;s.steps=state.step;
s.invocationWallSeconds=elapsed;s.accumulatedWallSeconds=priorWall+elapsed;
s.physicalTime=state.scale.physicalTime;s.physicalTimeLow=state.scale.physicalTimeLow;
s.rates=[state.rhs.c_l,state.rhs.c_omega,state.rhs.c_r];s.explicitOptions=method;
if isfield(method,'timeIntegrator'),s.timeIntegrator=method.timeIntegrator;else,s.timeIntegrator='ssprk3';end
s.coreSize=size(state.W.core);s.polarSize=size(state.W.polar);s.activeDofCount=g.activeDofCount;
s.acceptedCheckpoint=pointer;s.historyRows=size(history,1);s.historyColumns=names;
s.maximumAcceptedStageCfl=finite_max(history(:,strcmp(names,'maximumStageCfl')));
s.maximumRecordedGaugeResidual=finite_max(history(:,strcmp(names,'gaugeResidualInfinity')));
s.maximumRecordedPoissonBackwardError=finite_max(history(:,strcmp(names,'poissonBackwardError')));
s.maximumRecordedStrongResidual=finite_max(history(:,strcmp(names,'poissonStrongResidual')));
s.strongResidualAvailableRows=nnz(history(:,strcmp(names,'poissonStrongAvailable')));
s.totalStageCflRetries=sum(history(:,strcmp(names,'retries')));
s.totalRhsEvaluations=sum(history(:,strcmp(names,'rhsEvaluations')));
s.originalPeakEventLimitedSteps=sum(history(:,strcmp(names,'peakEventLimited')));
s.originalPeakEventsLocated=sum(history(:,strcmp(names,'peakEventLocated')));
s.originalPeakEventProbes=sum(history(:,strcmp(names,'peakEventProbeCount')));
s.originalPeakEventExtraRhsEvaluations=sum(history(:,strcmp(names,'peakEventExtraRhs')));
s.nonuniquePeakStageCount=sum(history(:,strcmp(names,'peakTieStages')));
curvature=history(:,strcmp(names,'minimumPeakCurvatureMagnitude'));
curvature=curvature(isfinite(curvature));s.minimumRecordedPeakCurvatureMagnitude=min(curvature);
s.finalDensityRange=history(end,ismember(names,{'totalDensityMin','totalDensityMax'}));
s.finalQInfinity=max([max(abs(state.rhs.Q.core),[],'all'),max(abs(state.rhs.Q.polar),[],'all')]);
s.targetReached=state.scale.canonicalTime>=s.targetTau;
s.finalPoissonDiagnostics=state.rhs.poisson;
end

function write_summary(folder,summary,history,historyColumns)
temporary=fullfile(folder,'summary.writing.mat');save(temporary,'summary','history','historyColumns','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'summary.mat'),'f');assert(ok,'ipm_cr:ExplicitSummary','%s',message);
ipm_cr.io.json(fullfile(folder,'summary.json'),summary);
end

function value=get(s,name,default)
if isfield(s,name),value=s.(name);else,value=default;end
end

function value=finite_max(values)
values=values(isfinite(values));if isempty(values),value=NaN;else,value=max(values);end
end
