function report=plot_ipm_cr_tau10_radial_diagnostic(studyDir,radialDir,outputDir)
%PLOT_IPM_CR_TAU10_RADIAL_DIAGNOSTIC Frozen-state dense radial display.
%   report=...(studyDir,radialDir,outputDir) uses three absolute directories.
%   STUDYDIR contains the original ivp_spatial.mat and level_2/level_3.
%   RADIALDIR contains the completed radial-only run and its accepted CP.
%   OUTPUTDIR must be fresh. No time advance, Poisson solve, field alignment
%   or numerical-package modification is performed.
%
% W = H(scale,X,Y) + I(eta); q = G_X(scale,X,Y) + I(sync(D_X eta)).
% Each trajectory retains its own saved background, complete scale/clock
% and native derivative/ghost operators. The denser display queries are
% not the original convergence-test query set or a new precision norm.
assert(nargin==3,'ipm_cr:RadialPlotArguments', ...
    'Supply explicit studyDir, radialDir and fresh outputDir.');
studyDir=char(studyDir);radialDir=char(radialDir);outputDir=char(outputDir);
assert(startsWith(studyDir,filesep)&&startsWith(radialDir,filesep)&& ...
    startsWith(outputDir,filesep)&&isfolder(studyDir)&&isfolder(radialDir), ...
    'ipm_cr:RadialPlotPaths','Use absolute paths and existing input directories.');
assert(~isfolder(outputDir)&&~isfile(outputDir)&& ...
    ~startsWith(outputDir,[studyDir,filesep])&&~startsWith(outputDir,[radialDir,filesep]), ...
    'ipm_cr:RadialPlotFresh','Use a fresh output directory outside both input directories.');
experimentDir=fileparts(studyDir);runtime=fullfile(experimentDir,'runtime');
previousPath=path;pathCleanup=onCleanup(@()path(previousPath)); %#ok<NASGU>
addpath(fullfile(runtime,'test'),fullfile(runtime,'time_evolution'),'-begin');
assert(strcmp(canonical(which('ipm_c.output.sha256')), ...
    canonical(fullfile(runtime,'time_evolution','+ipm_c','+output','sha256.m'))), ...
    'ipm_cr:RadialPlotFrozenPath','SHA and all numerical functions must use the frozen runtime.');
timer=tic;
manifestPath=fullfile(experimentDir,'source_inventory.json');
manifest=jsondecode(fileread(manifestPath));verify_manifest(manifest,experimentDir);
resolved=verify_resolution(manifest,experimentDir,runtime);
source=mfilename('fullpath');if ~endsWith(source,'.m'),source=[source,'.m'];end
sourceIdentity=file_identity(source);
inputFiles=[file_identity(fullfile(studyDir,'ivp_spatial.mat')), ...
    file_identity(fullfile(studyDir,'ivp_spatial.json')),file_identity(manifestPath), ...
    file_identity(fullfile(radialDir,'radial_refinement.mat')), ...
    file_identity(fullfile(radialDir,'radial_refinement.json'))];
loaded=load(inputFiles(1).path,'report');base=loaded.report;clear loaded
loaded=load(inputFiles(4).path,'report');radial=loaded.report;clear loaded
assert(numel(base.records)==3&&numel(base.sampledStates)==3&& ...
    base.records(2).run.tau==10&&base.records(3).run.tau==10&& ...
    radial.completed&&radial.tau==10&&radial.numericalSourcesUnchanged&& ...
    radial.originalReportUnchanged, ...
    'ipm_cr:RadialPlotEndpoint','Both original levels and the radial-only trajectory must be completed at tau=10.');
assert(strcmp(radial.originalReportMatSHA256,inputFiles(1).sha256)&& ...
    strcmp(radial.originalReportJsonSHA256,inputFiles(2).sha256)&& ...
    strcmp(canonical(radial.baseStudy),canonical(studyDir)), ...
    'ipm_cr:RadialPlotStudyIdentity','Radial trial must use these unchanged original study files.');
verify_source_set(base.sourceInventory,manifest,experimentDir);
verify_source_set(radial.sourceInventory,manifest,experimentDir);

% Use the exact physical angle of the original coarse-minus-middle q peak.
% The original query set and its numerical precision report remain intact.
cloud=base.queryCloud;
assert(isequaln(radial.queryCloud,cloud),'ipm_cr:RadialPlotQueryIdentity', ...
    'The radial replay must use the unchanged original physical query cloud.');
[~,peakIndex]=max(abs(base.sampledStates{1}.q-base.sampledStates{2}.q));
theta=[0,atan2(cloud.y(peakIndex),cloud.x(peakIndex))];
assert(theta(2)>0&&theta(2)<pi/2,'ipm_cr:RadialPlotAngle', ...
    'The original q-peak ray must be an interior positive-X ray.');
r=linspace(20,180,1601)';
X=r*cos(theta);Y=r*sin(theta);
folders={fullfile(studyDir,'level_2'),fullfile(studyDir,'level_3'),radialDir};
labels={'mid','fine','radial'};
expected={base.records(2).run.acceptedCheckpoint,base.records(3).run.acceptedCheckpoint, ...
    radial.run.acceptedCheckpoint};
references={base.sampledStates{2},base.sampledStates{3},radial.sampledState};
samples=cell(1,3);levels=struct([]);replay=struct([]);
commonAuthority=[];commonOriginalAuthority=[];
for j=1:3
    [cp,pointer,pointerFile]=load_accepted(folders{j});
    assert_same_pointer(pointer,expected{j});
    verify_source_set(cp.sourceInventory,manifest,experimentDir);
    state=cp.state;g=cp.geometry;method=cp.explicitOptions;
    assert(isfield(state,'residual')&&isfield(method,'background'), ...
        'ipm_cr:RadialPlotRepresentation','Use saved eta plus the saved analytic background.');
    if j==1
        commonAuthority=cp.authority;commonOriginalAuthority=state.originalAuthority;
        midMethod=method;
    else
        assert(isequaln(cp.authority,commonAuthority)&& ...
            isequaln(state.originalAuthority,commonOriginalAuthority), ...
            'ipm_cr:RadialPlotAuthority','All three trajectories must retain the same original wall authority.');
    end
    if j==3
        assert(isequaln(method,radial.method)&&isequaln(cp.authority,radial.authority), ...
            'ipm_cr:RadialPlotMethod','Saved radial checkpoint must match its completed report.');
        unchanged=method;
        unchanged.geometry.grid.nr=midMethod.geometry.grid.nr;
        unchanged.geometry.grid.firstSpacing=midMethod.geometry.grid.firstSpacing;
        assert(isequaln(unchanged,midMethod), ...
            'ipm_cr:RadialPlotOnlyRadialChange','Only nr and firstSpacing may differ from the original middle method.');
    end
    % Preserve the original full arrays and their ghost collars. No fresh
    % evaluate/evaluateResidual, synchronization of eta, or solve is needed.
    d=ipm_cr.explicit.differentiate(state.residual,g);
    qResidual=ipm_cr.explicit.synchronize(struct('core',d.core.x,'polar',d.polar.x),g);
    samples{j}=sample_state(state,g,method,qResidual,X,Y);
    original=sample_state(state,g,method,qResidual,cloud.x,cloud.y);
    check=struct('trajectory',labels{j},'originalQueryCount',numel(cloud.x), ...
        'W',replay_check(original.W,references{j}.W), ...
        'q',replay_check(original.q,references{j}.q));
    replay=append(replay,check);
    entry=struct('trajectory',labels{j},'checkpoint',pointer, ...
        'pointerFile',pointerFile,'tau',state.scale.canonicalTime,'step',state.step, ...
        'sourceInventory',cp.sourceInventory,'backgroundModel',method.background, ...
        'scale',state.scale,'explicitOptions',method,'authority',cp.authority, ...
        'originalAuthority',state.originalAuthority, ...
        'coreSpacing',[g.core.dx,g.core.dy], ...
        'polarNodeCount',numel(g.polar.r),'angularNodeCount',numel(g.polar.theta), ...
        'radialNodesInDisplay',g.polar.r(g.polar.r>=r(1)&g.polar.r<=r(end)), ...
        'interpolationWidth',samples{j}.interpolationWidth, ...
        'interpolationOperatorInfinityNorm',samples{j}.interpolationOperatorInfinityNorm);
    levels=append(levels,entry);
    clear cp state g method d qResidual original
end
differences=struct('W_midMinusFine',samples{1}.W-samples{2}.W, ...
    'W_radialMinusFine',samples{3}.W-samples{2}.W, ...
    'q_midMinusFine',samples{1}.q-samples{2}.q, ...
    'q_radialMinusFine',samples{3}.q-samples{2}.q);
assert(all(isfinite([differences.W_midMinusFine(:);differences.W_radialMinusFine(:); ...
    differences.q_midMinusFine(:);differences.q_radialMinusFine(:)])), ...
    'ipm_cr:RadialPlotFinite','Every dense display query must remain finite; no samples are deleted.');
displayMetrics=struct([]);
for ray=1:2
    item=struct('thetaRadians',theta(ray),'thetaDegrees',theta(ray)*180/pi);
    for name=fieldnames(differences)'
        key=name{1};value=differences.(key);item.(key)=sample_metric(value(:,ray),r);
    end
    displayMetrics=append(displayMetrics,item);
end

mkdir(outputDir);
csvPath=fullfile(outputDir,'radial_samples.csv');
nodeCsvPath=fullfile(outputDir,'native_radial_nodes.csv');
write_samples(csvPath,r,theta,X,Y,samples,differences);
write_nodes(nodeCsvPath,levels);
[plotFiles,rugGroups]=make_plot(outputDir,r,theta,differences,levels);
caption=['At canonical tau=10, signed differences from the original full fine trajectory ', ...
    'on two physical rays, sampled at 1601 equally spaced radii r=20..180. ', ...
    'Mid-minus-fine changes core, radial and angular resolution; radial-minus-fine ', ...
    'retains the middle core/angular grid but uses the fine radial grid. ', ...
    'Each saved scale/background is evaluated directly; W=H+I(eta), ', ...
    'q=G_X+I(sync(D_X eta)), using the frozen production seven-point interpolation ', ...
    'and native constant-preserving derivative. No field fitting, shifting or alignment ', ...
    'is applied. Short gray marks at the top are actual native radial nodes; ', ...
    'their vertical positions are only a visual reference. These are dense display ', ...
    'samples, not the original convergence-test norms or an exact-solution error. ', ...
    'Physical r alone does not establish that the transported body is in its asymptotic tail. ', ...
    'All plotted quantities and coordinates are dimensionless rescaled variables.'];
captionPath=fullfile(outputDir,'radial_diagnostic_caption.txt');
write_text(captionPath,caption);
report=struct('kind','frozen_tau10_dense_radial_display_v1','status','completed', ...
    'studyDir',studyDir,'radialDir',radialDir,'outputDir',outputDir,'tau',10, ...
    'script',sourceIdentity,'inputFiles',inputFiles,'levels',levels, ...
    'frozenRuntime',runtime,'frozenManifestFile',manifestPath, ...
    'verifiedFrozenFileCount',numel(manifest),'resolvedNumericalFunctions',resolved, ...
    'originalQueryReplay',replay,'displaySampleCountPerRay',numel(r), ...
    'displaySampleCountTotal',numel(X),'radiusInterval',[r(1),r(end)], ...
    'thetaRadians',theta,'thetaDegrees',theta*180/pi, ...
    'secondRayDefinition','Exact angle of the original coarse-minus-middle q maximum on its original query cloud.', ...
    'secondRayOriginalQuery',struct('oneBasedIndex',peakIndex, ...
        'X',cloud.x(peakIndex),'Y',cloud.y(peakIndex)), ...
    'fieldDefinition','W=background.H+interpolate(eta); q=background.x+interpolate(synchronize(native Dx eta)).', ...
    'units','Dimensionless rescaled coordinates, centered density W and rescaled density derivative q.', ...
    'reference','Original finite-resolution full fine trajectory, not an exact solution.', ...
    'displayMetricsOnly',displayMetrics,'nativeNodeRugGroups',rugGroups, ...
    'allDisplayQueriesRetained',true,'fieldAlignmentOrFitApplied',false, ...
    'originalNormsReplaced',false,'convergenceOrderEstimated',false, ...
    'timeAdvanced',0,'poissonRecomputed',false,'numericalStateChanged',false, ...
    'style','journal-minimal; 7.0 x 4.9 inch, 9 point serif, paired color/line styles, zero lines, no overall title', ...
    'visualInspectionCompleted',false, ...
    'visualInspectionNote','The helper exports PNG/PDF; inspect those files at final size after execution.', ...
    'files',struct('samplesCsv',csvPath,'nativeNodesCsv',nodeCsvPath, ...
        'samplesMat',fullfile(outputDir,'radial_diagnostic.mat'), ...
        'png',plotFiles.png,'pdf',plotFiles.pdf,'caption',captionPath, ...
        'json',fullfile(outputDir,'radial_diagnostic.json')), ...
    'caption',caption, ...
    'limits',{{'Dense one-dimensional display sampling is not the original precision norm or a continuum accuracy certificate.', ...
        'Radial-only sensitivity versus the full fine grid is not an isotropic refinement order.', ...
        'Both rays retain the physical coordinates and the separate saved scale/background of each trajectory.', ...
        'A visual oscillation does not by itself identify a transport, ghost, interpolation or Poisson cause.', ...
        'This tau10 display does not qualify a tau200 trajectory.'}});
verify_manifest(manifest,experimentDir);
verify_file_identities(inputFiles);verify_file_identities(sourceIdentity);
for j=1:numel(levels)
    assert(strcmp(ipm_c.output.sha256(levels(j).checkpoint.path),levels(j).checkpoint.sha256), ...
        'ipm_cr:RadialPlotInputChanged','Accepted checkpoint changed during plotting.');
    verify_file_identities(levels(j).pointerFile);
end
assert(isequaln(resolved,verify_resolution(manifest,experimentDir,runtime)), ...
    'ipm_cr:RadialPlotSourceChanged','Actual numerical resolutions changed during plotting.');
report.inputsAndSourcesUnchanged=true;report.elapsedSeconds=toc(timer);
save(report.files.samplesMat,'report','r','theta','X','Y','samples','differences','-v7.3');
write_text(report.files.json,jsonencode(report,'PrettyPrint',true));
fprintf('TAU10 RADIAL DISPLAY tau=10 rays=2 points/ray=%d elapsed=%.3fs %s\n', ...
    numel(r),report.elapsedSeconds,outputDir);
end

function values=sample_state(state,g,method,qResidual,X,Y)
background=ipm_cr.explicit.background(method.background,state.scale,X,Y);
[eta,map]=ipm_cr.explicit.interpolate(state.residual,g,X,Y);
assert(map.allValid,'ipm_cr:RadialPlotQueries','All physical display/replay queries must be valid.');
dxEta=reshape(map.core*qResidual.core(:)+map.polar*qResidual.polar(:),size(X));
values=struct('W',background.H+eta,'q',background.x+dxEta, ...
    'backgroundH',background.H,'eta',eta,'backgroundX',background.x,'dxEta',dxEta, ...
    'interpolationWidth',map.interpolationWidth, ...
    'interpolationOperatorInfinityNorm',map.operatorInfinityNorm);
for name={'W','q','backgroundH','eta','backgroundX','dxEta'}
    assert(all(isfinite(values.(name{1})),'all'),'ipm_cr:RadialPlotFinite', ...
        'Nonfinite reconstructed value at a retained query.');
end
end

function info=replay_check(actual,recorded)
assert(isequal(size(actual),size(recorded)),'ipm_cr:RadialPlotReplay','Original query shapes must agree.');
d=max(abs(actual-recorded),[],'all');
amplitude=max([1;abs(actual(:));abs(recorded(:))]);bound=64*eps(amplitude);
info=struct('absoluteInfinity',d,'bitwiseIdentical',isequaln(actual,recorded), ...
    'roundoffReplayBound',bound,'boundDefinition','64*eps(max(1,maximum magnitude of actual and recorded field))', ...
    'withinReplayBound',d<=bound);
assert(info.withinReplayBound,'ipm_cr:RadialPlotReplay', ...
    'Frozen production reconstruction disagrees with the saved original sample (%.17g > %.17g).',d,bound);
end

function [cp,pointer,pointerFile]=load_accepted(folder)
pointerFile=file_identity(fullfile(folder,'latest.json'));
pointer=jsondecode(fileread(pointerFile.path));
assert(all(isfield(pointer,{'kind','path','sha256','sequence','tau'}))&& ...
    strcmp(pointer.kind,'ipm_cr_explicit_accepted_pointer_v1')&&pointer.tau==10&& ...
    any(strcmp(pointer.path,{fullfile(folder,'accepted_0.mat'),fullfile(folder,'accepted_1.mat')}))&& ...
    strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256), ...
    'ipm_cr:RadialPlotCheckpoint','Require an intact accepted tau10 pointer in its own two-slot folder.');
loaded=load(pointer.path,'checkpoint');cp=loaded.checkpoint;
assert(strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256)&& ...
    strcmp(cp.kind,'ipm_cr_explicit_accepted_checkpoint_v1')&&cp.acceptedOnly&& ...
    cp.sequence==pointer.sequence&&cp.state.step==pointer.sequence&& ...
    cp.state.scale.canonicalTime==10&&isequaln(cp.rho,cp.state.rho), ...
    'ipm_cr:RadialPlotCheckpoint','Accepted checkpoint and pointer must agree before and after loading.');
verify_file_identities(pointerFile);
end

function assert_same_pointer(actual,expected)
for name={'kind','path','sha256','sequence','tau'}
    key=name{1};
    assert(isfield(expected,key)&&isequaln(actual.(key),expected.(key)), ...
        'ipm_cr:RadialPlotRecordedCheckpoint','The accepted checkpoint must equal the endpoint recorded by its completed experiment.');
end
end

function verify_manifest(manifest,experimentDir)
assert(isstruct(manifest)&&~isempty(manifest)&&all(isfield(manifest,{'path','sha256'})), ...
    'ipm_cr:RadialPlotManifest','Require a complete frozen source/input manifest.');
for j=1:numel(manifest)
    p=manifest_path(manifest(j).path,experimentDir);
    assert(isfile(p)&&strcmp(ipm_c.output.sha256(p),manifest(j).sha256), ...
        'ipm_cr:RadialPlotFrozenSource','Frozen input/source changed: %s',p);
end
end

function verify_source_set(inventory,manifest,experimentDir)
paths=arrayfun(@(s)canonical(manifest_path(s.path,experimentDir)),manifest,'UniformOutput',false);
assert(isstruct(inventory)&&~isempty(inventory)&&all(isfield(inventory,{'path','sha256'})), ...
    'ipm_cr:RadialPlotSourceInventory','Each accepted checkpoint/report must include source SHA identities.');
for j=1:numel(inventory)
    p=canonical(inventory(j).path);index=find(strcmp(paths,p));
    assert(isscalar(index)&&isfile(p)&&strcmp(inventory(j).sha256,manifest(index).sha256)&& ...
        strcmp(ipm_c.output.sha256(p),inventory(j).sha256), ...
        'ipm_cr:RadialPlotSourceInventory','Recorded source is changed or absent from the frozen manifest: %s',p);
end
end

function resolved=verify_resolution(manifest,experimentDir,runtime)
symbols={'ipm_c.output.sha256','ipm_cr.explicit.background', ...
    'ipm_cr.explicit.gaussianBackground','ipm_cr.explicit.dipoleTail', ...
    'ipm_cr.explicit.interpolate','ipm_cr.explicit.differentiate', ...
    'ipm_cr.explicit.synchronize','ipm_cr.carrier.tensorMap','ipm_cr.carrier.localMap'};
paths=arrayfun(@(s)canonical(manifest_path(s.path,experimentDir)),manifest,'UniformOutput',false);
resolved=struct([]);
for j=1:numel(symbols)
    p=canonical(which(symbols{j}));index=find(strcmp(paths,p));
    assert(isscalar(index)&&startsWith(p,[canonical(runtime),filesep])&& ...
        strcmp(ipm_c.output.sha256(p),manifest(index).sha256), ...
        'ipm_cr:RadialPlotResolution','Actual function %s must resolve to its SHA-verified frozen runtime implementation.',symbols{j});
    resolved=append(resolved,struct('symbol',symbols{j},'path',p,'sha256',manifest(index).sha256));
end
end

function p=manifest_path(p,folder)
p=char(p);if ~startsWith(p,filesep),p=fullfile(folder,p);end
end

function p=canonical(p)
assert(~isempty(p),'ipm_cr:RadialPlotPathResolution','A required path cannot be empty.');
p=char(java.io.File(char(p)).getCanonicalPath());
end

function identity=file_identity(p)
assert(isfile(p),'ipm_cr:RadialPlotInput','Required input does not exist: %s',p);
identity=struct('path',canonical(p),'sha256',ipm_c.output.sha256(p));
end

function verify_file_identities(identities)
for j=1:numel(identities)
    assert(isfile(identities(j).path)&&strcmp(ipm_c.output.sha256(identities(j).path),identities(j).sha256), ...
        'ipm_cr:RadialPlotInputChanged','Input/helper changed during plotting: %s',identities(j).path);
end
end

function values=append(values,item)
if isempty(values),values=item;else,values(end+1)=item;end
end

function item=sample_metric(d,r)
[value,k]=max(abs(d));
item=struct('sampledAbsoluteInfinity',value,'sampledRms',sqrt(mean(d.^2)), ...
    'radiusAtMaximum',r(k),'signedDifferenceAtMaximum',d(k),'sampleCount',numel(d));
end

function write_samples(file,r,theta,X,Y,samples,d)
% One row for every retained physical query; no downsampling or smoothing.
n=numel(X);ray=repelem((1:2)',numel(r));rColumn=repmat(r,2,1);
T=table(ray,rColumn,repelem(theta(:),numel(r)),X(:),Y(:), ...
    'VariableNames',{'ray','r','thetaRadians','X','Y'});
names={'mid','fine','radial'};
for j=1:3
    for field={'W','q','backgroundH','eta','backgroundX','dxEta'}
        key=field{1};T.([key,'_',names{j}])=reshape(samples{j}.(key),n,1);
    end
end
for field=fieldnames(d)'
    key=field{1};T.(key)=reshape(d.(key),n,1);
end
writetable(T,file);
end

function write_nodes(file,levels)
trajectory=strings(0,1);radius=zeros(0,1);
for j=1:numel(levels)
    v=levels(j).radialNodesInDisplay(:);
    trajectory=[trajectory;repmat(string(levels(j).trajectory),numel(v),1)]; %#ok<AGROW>
    radius=[radius;v]; %#ok<AGROW>
end
writetable(table(trajectory,radius),file);
end

function [files,groups]=make_plot(folder,r,theta,d,levels)
fig=figure('Visible','off','Color','white','Units','inches', ...
    'Position',[1,1,7,4.9],'InvertHardcopy','off');
figureCleanup=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors=[31,78,121;163,74,40]/255;styles={'-','--'};
groups=struct([]);
for j=1:numel(levels)
    matched=[];
    for k=1:numel(groups)
        if isequaln(groups(k).radialNodes,levels(j).radialNodesInDisplay),matched=k;break;end
    end
    if isempty(matched)
        groups=append(groups,struct('trajectories',{{levels(j).trajectory}}, ...
            'radialNodes',levels(j).radialNodesInDisplay));
    else
        groups(matched).trajectories{end+1}=levels(j).trajectory;
    end
end
axesList=gobjects(1,4);lineHandles=gobjects(1,2);
for ray=1:2
    for quantity=1:2
        panel=(ray-1)*2+quantity;ax=nexttile(layout);axesList(panel)=ax;hold(ax,'on');
        if quantity==1
            a=d.W_midMinusFine(:,ray);b=d.W_radialMinusFine(:,ray);
            allValues=[d.W_midMinusFine(:);d.W_radialMinusFine(:)];label='\Delta W';
        else
            a=d.q_midMinusFine(:,ray);b=d.q_radialMinusFine(:,ray);
            allValues=[d.q_midMinusFine(:);d.q_radialMinusFine(:)];label='\Delta q';
        end
        yline(ax,0,'-','Color',[.55,.55,.55],'LineWidth',.65,'HandleVisibility','off');
        h1=plot(ax,r,a,'Color',colors(1,:),'LineStyle',styles{1},'LineWidth',1.3);
        h2=plot(ax,r,b,'Color',colors(2,:),'LineStyle',styles{2},'LineWidth',1.3);
        if panel==1,lineHandles=[h1,h2];end
        lo=min([0;allValues]);hi=max([0;allValues]);span=hi-lo;
        if span==0,span=eps;lo=-span/2;hi=span/2;end
        bounds=[lo-.06*span,hi+.18*span];ylim(ax,bounds);
        % Node rugs occupy reserved top padding and are not sampled values.
        for k=1:numel(groups)
            nodes=groups(k).radialNodes(:)';
            high=bounds(2)-(.025+.038*(k-1))*diff(bounds);
            low=high-.022*diff(bounds);
            plot(ax,reshape([nodes;nodes;NaN(size(nodes))],[],1), ...
                reshape([low+zeros(size(nodes));high+zeros(size(nodes));NaN(size(nodes))],[],1), ...
                'Color',[.54,.54,.54],'LineWidth',.55,'HandleVisibility','off');
        end
        xlim(ax,[r(1),r(end)]);xticks(ax,20:40:180);
        xlabel(ax,'r (rescaled)');ylabel(ax,label);
        title(ax,sprintf('(%c) \\theta = %.5g^\\circ','a'+panel-1,theta(ray)*180/pi), ...
            'FontWeight','normal','Interpreter','tex');
        set(ax,'FontName','Times New Roman','FontSize',9,'LineWidth',.7, ...
            'TickDir','out','Box','off','XColor',[.2,.2,.2], ...
            'YColor',[.2,.2,.2],'Color','white');
        ax.XLabel.FontSize=10;ax.YLabel.FontSize=10;ax.Title.FontSize=10;
        ax.XAxis.Exponent=0;
    end
end
legendHandle=legend(axesList(1),lineHandles,{'mid - fine','radial - fine'}, ...
    'Orientation','horizontal','Box','off','FontName','Times New Roman','FontSize',9);
legendHandle.Layout.Tile='north';
drawnow;
files=struct('png',fullfile(folder,'radial_diagnostic.png'), ...
    'pdf',fullfile(folder,'radial_diagnostic.pdf'));
fig.PaperUnits='inches';fig.PaperSize=fig.Position(3:4);
fig.PaperPosition=[0,0,fig.Position(3:4)];fig.PaperPositionMode='manual';
print(fig,files.pdf,'-dpdf','-painters');
print(fig,files.png,'-dpng','-r300');
for k=1:numel(groups)
    groups(k).displayOnlyRugRowFromTop=k;
    groups(k).note='Rug height is a visual reference, not a field value.';
end
end

function write_text(file,value)
fid=fopen(file,'w');assert(fid>=0,'ipm_cr:RadialPlotOutput','Cannot write %s',file);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,value,'char');
end
