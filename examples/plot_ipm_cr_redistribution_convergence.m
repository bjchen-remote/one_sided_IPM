function report=plot_ipm_cr_redistribution_convergence(oldFile,newFile,validationFile,outputDir)
%PLOT_IPM_CR_REDISTRIBUTION_CONVERGENCE Completed tau10 JSON reports only.
% Optional arguments: absolute old/new IVP JSON, validation JSON, fresh output
% directory. With no arguments, use the named experiments below. No numerical
% package, checkpoint loading, field interpolation or Poisson solve is used.
root=fileparts(fileparts(mfilename('fullpath')));results=fullfile(root,'test','results');
if nargin<1||isempty(oldFile),oldFile=fullfile(results,'tau10_order_20261022_1323','study','ivp_spatial.json');end
if nargin<2||isempty(newFile),newFile=fullfile(results,'radial_redistribution_tau10_20261022','study','ivp_spatial.json');end
if nargin<3||isempty(validationFile),validationFile=fullfile(results,'radial_redistribution_tau10_20261022','validation','radial_redistribution.json');end
if nargin<4||isempty(outputDir),outputDir=fullfile(results,'radial_redistribution_tau10_20261022','figures','redistribution_convergence');end
oldFile=char(oldFile);newFile=char(newFile);validationFile=char(validationFile);outputDir=char(outputDir);
assert(all(cellfun(@(p)java.io.File(p).isAbsolute(),{oldFile,newFile,validationFile,outputDir})), ...
    'ipm_cr:RedistributionPlotPaths','Use absolute input and output paths.');
timer=tic;source=mfilename('fullpath');if ~endsWith(source,'.m'),source=[source,'.m'];end
scriptIdentity=file_identity(source);
inputFiles=[file_identity(oldFile),file_identity(newFile),file_identity(validationFile)];
old=jsondecode(fileread(oldFile));new=jsondecode(fileread(newFile));validation=jsondecode(fileread(validationFile));
assert(isfield(new,'sampledStates')&&numel(new.records)==3&&numel(new.sampledStates)==3&&numel(new.comparisons)==2, ...
    'ipm_cr:RedistributionPlotIncomplete','New three-level tau10 study is incomplete; wait for records, sampledStates and comparisons.');
studies={read_study(old,'original'),read_study(new,'redistributed')};base=studies{1};refined=studies{2};
assert(isequaln(base.levels,refined.levels)&&isequaln(old.options.solverOptions,new.options.solverOptions), ...
    'ipm_cr:RedistributionPlotConfiguration','Requested levels and Poisson options must agree.');
for key={'tauEnd','maxDt','cfl'}
    assert(isequaln(old.options.(key{1}),new.options.(key{1})), ...
        'ipm_cr:RedistributionPlotConfiguration','Time option differs: %s.',key{1});
end
assert(isequaln(base.x,refined.x)&&isequaln(base.y,refined.y), ...
    'ipm_cr:RedistributionPlotCloud','Original and new query coordinates must agree exactly.');
for key={'near','wall','hotspot','interfaceBand','remoteWallStrip'}
    assert(isequaln(base.masks.(key{1}),refined.masks.(key{1})), ...
        'ipm_cr:RedistributionPlotCloud','Query membership differs: %s.',key{1});
end
[convergence,orders,plotDifferences]=convergence_data(studies);
[axesData,gaps]=radial_axes(validation,base.levels);
for j=1:3
    if isfield(refined.records(j),'actualGeometry')
        assert(isequal(refined.records(j).actualGeometry.polarSize(:)', ...
            [axesData(j).geometry.nr,axesData(j).geometry.ntheta]), ...
            'ipm_cr:RedistributionPlotGeometry','Study and validation annulus sizes differ.');
    end
end
qReference=base.states(3).q;difference=[base.states(2).q-qReference,refined.states(2).q-qReference];
assert(all(isfinite(difference),'all'),'ipm_cr:RedistributionPlotFinite','Every common-cloud difference must be finite.');
qSamples=table((1:base.count)',base.x,base.y,base.r,atan2(base.y,base.x), ...
    'VariableNames',{'sample_index_one_based','X','Y','r','theta'});
for key=fieldnames(base.masks)',qSamples.(key{1})=base.masks.(key{1});end
qSamples.q_original_middle=base.states(2).q;qSamples.q_redistributed_middle=refined.states(2).q;
qSamples.q_original_fine_reference=qReference;
qSamples.original_signed_difference=difference(:,1);qSamples.redistributed_signed_difference=difference(:,2);
qSamples.original_absolute_difference=abs(difference(:,1));qSamples.redistributed_absolute_difference=abs(difference(:,2));
maxima=struct();names={'original','redistributed'};
for j=1:2
    [value,k]=max(abs(difference(:,j)));
    maxima.(names{j})=struct('sampledAbsoluteInfinity',value,'sampleIndexOneBased',k, ...
        'X',base.x(k),'Y',base.y(k),'r',base.r(k));
end
[performance,timingFiles]=performance_data(studies,axesData);inputFiles=[inputFiles,timingFiles];
assert(~isfolder(outputDir)&&~isfile(outputDir),'ipm_cr:RedistributionPlotFresh','Use a fresh output directory: %s',outputDir);
mkdir(outputDir);
write_csv(fullfile(outputDir,'convergence.csv'),struct2table(convergence));
write_csv(fullfile(outputDir,'radial_gaps.csv'),gaps);
write_csv(fullfile(outputDir,'q_common_cloud.csv'),qSamples);
write_csv(fullfile(outputDir,'performance.csv'),struct2table(performance));
[figureFiles,settings]=make_plot(outputDir,studies,axesData,plotDifferences,difference,performance);
caption=make_caption(studies,validation,orders,settings,axesData);
captionFile=fullfile(outputDir,'caption.txt');write_text(captionFile,caption);
for j=1:numel(inputFiles),verify_identity(inputFiles(j));end
verify_identity(scriptIdentity);
qualification=struct('original',qualification_flags(old),'redistributed',qualification_flags(new));
validationFlags=struct();
for key={'status','passed','structuralChecksPassed','spatialQualificationPassed', ...
        'allResolvedOrdersAtLeastFour','fourthOrderEvidenceResolvedForAllPairs'}
    validationFlags.(key{1})=get_option(validation,key{1},[]);
end
inventories=struct();
for j=1:2
    inventories.(names{j})=struct('count',numel(studies{j}.raw.sourceInventory), ...
        'savedInventoryJsonSHA256',hash_text(jsonencode(studies{j}.raw.sourceInventory)));
end
report=struct('kind','ipm_cr_radial_redistribution_figure_v2','status','completed', ...
    'script',scriptIdentity,'inputFiles',inputFiles,'outputDir',outputDir,'tau',10, ...
    'queryCount',base.count,'queryRadiusRange',[min(base.r),max(base.r)], ...
    'allQueriesRetained',true,'exteriorMaskApplied',false,'errorInterpolationPerformed',false, ...
    'continuousGlobalInfinityBoundClaimed',false,'observedSelfDifferenceOrders',orders, ...
    'qMiddleVersusOriginalFineMaxima',maxima,'qualification',qualification, ...
    'validationQualification',validationFlags,'sourceInventories',inventories, ...
    'style','journal-minimal; Times New Roman; recorded points only; no overall title', ...
    'plotSettings',settings,'files',figureFiles,'captionFile',captionFile, ...
    'visualInspectionCompleted',false,'matlabVersion',version,'timeAdvanced',0, ...
    'poissonRecomputed',false,'inputsAndScriptUnchanged',true,'elapsedSeconds',toc(timer));
write_text(fullfile(outputDir,'provenance.json'),jsonencode(report,'PrettyPrint',true));
files=dir(outputDir);manifest=struct([]);
for j=1:numel(files)
    if ~files(j).isdir
        entry=file_identity(fullfile(outputDir,files(j).name));entry.path=files(j).name;manifest=append(manifest,entry);
    end
end
write_csv(fullfile(outputDir,'sha256.csv'),struct2table(manifest));
fprintf('REDISTRIBUTION FIGURE tau=10 queries=%d elapsed=%.3fs %s\n',base.count,report.elapsedSeconds,outputDir);
end

function s=read_study(raw,label)
assert(all(isfield(raw,{'levels','records','sampledStates','comparisons','queryCloud','sourcesUnchanged'}))&& ...
    numel(raw.levels)==3&&numel(raw.records)==3&&numel(raw.sampledStates)==3&&numel(raw.comparisons)==2, ...
    'ipm_cr:RedistributionPlotIncomplete','%s study must contain three complete spatial levels and two comparisons.',label);
assert(~any(strcmp(raw.status,{'running','initializing'}))&&raw.sourcesUnchanged, ...
    'ipm_cr:RedistributionPlotSources','%s report must be finished with a passing source-inventory check.',label);
cloud=raw.queryCloud;x=finite_vector(cloud.x,numel(cloud.x));y=finite_vector(cloud.y,numel(x));
assert(cloud.allQueriesRetained&&~cloud.exteriorMaskApplied,'ipm_cr:RedistributionPlotCloud', ...
    'The full original cloud must be retained without an exterior mask.');
masks=struct();
for key={'near','wall','hotspot','interfaceBand','remoteWallStrip'}
    mask=cloud.(key{1})(:);
    assert(numel(mask)==numel(x)&&all(mask==0|mask==1),'ipm_cr:RedistributionPlotCloud','Invalid query mask.');
    masks.(key{1})=logical(mask);
end
levels=struct_array(raw.levels);records=struct_array(raw.records);samples=struct_array(raw.sampledStates);states=struct([]);
for j=1:3
    run=records(j).run;
    assert(isequaln(records(j).level,levels(j))&&run.completed&&abs(run.tau-10)<=1e-12&& ...
        run.seconds>0&&isfinite(run.seconds)&&run.steps>0&&run.fullExplicitCflRequired&& ...
        isfinite(run.maximumStageCfl)&&run.maximumStageCfl<=run.configuredCfl*(1+64*eps), ...
        'ipm_cr:RedistributionPlotEndpoint','%s level %d lacks a valid tau10 endpoint/full-stage CFL.',label,j);
    item=struct();
    for key={'W','q','F','Q'},item.(key{1})=finite_vector(samples(j).(key{1}),numel(x));end
    states=append(states,item);
end
h=[levels.h];
assert(all(h>0)&&all(diff(h)<0)&&abs(h(1)/h(2)-h(2)/h(3))<1e-12, ...
    'ipm_cr:RedistributionPlotLevels','Use three positive h values with equal refinement ratios.');
s=struct('raw',raw,'label',label,'levels',levels,'records',records,'states',states, ...
    'comparisons',struct_array(raw.comparisons),'x',x,'y',y,'r',hypot(x,y),'masks',masks,'h',h,'count',numel(x));
end

function [rows,orders,plotData]=convergence_data(studies)
rows=struct([]);orders=struct();plotData=struct();
for method=1:2
    s=studies{method};values=struct();p=struct();
    for name={'W','q','F','Q'}
        key=name{1};d=zeros(1,2);
        for j=1:2
            actual=max(abs(s.states(j).(key)-s.states(j+1).(key)));saved=s.comparisons(j).fields.(key).global;
            assert(saved.sampleCount==s.count&&isfinite(saved.absoluteInfinity)&&saved.absoluteInfinity>=0&& ...
                abs(saved.absoluteInfinity-actual)<=64*eps(max([saved.absoluteInfinity,actual,realmin])), ...
                'ipm_cr:RedistributionPlotMetric','%s/%s saved metric disagrees with complete sampled arrays.',s.label,key);
            d(j)=saved.absoluteInfinity;loc=saved.maximumLocation;
            rows=append(rows,struct('method',s.label,'field',key,'pair',j,'h_coarser',s.h(j), ...
                'h_finer',s.h(j+1),'query_count',s.count,'recorded_absolute_infinity',d(j), ...
                'recomputed_absolute_infinity',actual,'recorded_rms',saved.rms, ...
                'maximum_X',loc.X,'maximum_Y',loc.Y,'maximum_r',loc.r));
        end
        values.(key)=d;
        if all(d>0),p.(key)=log(d(1)/d(2))/log(s.h(1)/s.h(2));else,p.(key)=[];end
    end
    orders.(s.label)=p;plotData.(s.label)=values;
end
end

function [data,T]=radial_axes(validation,levels)
records=struct_array(validation.records);
assert(numel(records)==3&&validation.structuralChecksPassed,'ipm_cr:RedistributionPlotValidation', ...
    'Require three validation grids with passing structural checks.');
data=struct([]);method=strings(0,1);level=[];spacing=[];index=[];left=[];right=[];retained=false(0,1);
for j=1:3
    record=records(j);g=record.geometry;info=g.radialRedistribution;
    assert(record.h==levels(j).h&&record.requestedNr==levels(j).nr&& ...
        record.requestedTheta==levels(j).ntheta&&info.enabled,'ipm_cr:RedistributionPlotGeometry','Validation levels differ.');
    base=finite_vector(info.baseR,info.baseNr);u=finite_vector(info.oldIndexAtNewNodes,info.actualNr);
    map=info.analyticMap;s=u/info.baseIntervals;
    if map.alpha==0
        tail=s.^map.power;
    else
        tail=zeros(size(s));positive=s>0;tail(positive)=exp(map.power*(log_sinh(map.alpha*s(positive))-log_sinh(map.alpha)));
    end
    r=map.a*s+map.tailCoefficient*tail;
    r(1:info.preservedPrefixCount)=base(1:info.preservedPrefixCount);r(1)=0;r(end)=base(end);
    assert(all(diff(base)>0)&&all(diff(r)>0),'ipm_cr:RedistributionPlotGeometry','Saved radial axes must increase.');
    [~,k]=min(abs(r-97));probe=g.spacingAt97;
    assert(k<numel(r)&&abs(r(k)-probe.r)<=2e-12*max(abs(probe.r),1)&& ...
        abs(r(k+1)-r(k)-probe.rightGap)<=2e-12*max(abs(probe.rightGap),1), ...
        'ipm_cr:RedistributionPlotGeometry','Recovered saved nodes disagree with spacingAt97.');
    data=append(data,struct('base',base,'new',r,'geometry',g,'oldGeometry',record.oldGeometry));
    for m=1:2
        if m==1,v=base;name='original';halo=record.oldGeometry.innerHaloRadius;
        else,v=r;name='redistributed';halo=g.innerHaloRadius;end
        n=numel(v)-1;method=[method;repmat(string(name),n,1)]; %#ok<AGROW>
        level=[level;repmat(j,n,1)];spacing=[spacing;repmat(levels(j).h,n,1)]; %#ok<AGROW>
        index=[index;(1:n)'];left=[left;v(1:end-1)];right=[right;v(2:end)];retained=[retained;v(1:end-1)>=halo]; %#ok<AGROW>
    end
end
T=table(method,level,spacing,index,left,right,(left+right)/2,right-left,retained, ...
    'VariableNames',{'method','level','h','left_full_index_one_based','r_left','r_right','r_midpoint','physical_neighbor_gap','inside_retained_annulus'});
end

function [rows,identities]=performance_data(studies,axesData)
rows=struct([]);identities=struct([]);
for method=1:2
    s=studies{method};
    for j=1:3
        run=s.records(j).run;identity=file_identity(fullfile(run.outputDir,'summary.json'));
        summary=jsondecode(fileread(identity.path));identities=append(identities,identity);
        assert(summary.tau==10&&summary.steps==run.steps&&isfinite(summary.accumulatedWallSeconds)&& ...
            summary.accumulatedWallSeconds>0&&abs(summary.invocationWallSeconds-run.seconds)<=64*eps(max(run.seconds,1)), ...
            'ipm_cr:RedistributionPlotTiming','Summary and recorded endpoint/timing disagree.');
        if method==1,g=axesData(j).oldGeometry;else,g=axesData(j).geometry;end
        rows=append(rows,struct('method',s.label,'level',j,'h',s.h(j),'requested_nr',s.levels(j).nr, ...
            'actual_full_axis_nr',g.actualFullAxisNr,'actual_annulus_nr',g.nr,'ntheta',s.levels(j).ntheta, ...
            'tau',run.tau,'elapsed_total_seconds',summary.accumulatedWallSeconds, ...
            'invocation_seconds',run.seconds,'resumed',logical(summary.resumed),'steps',run.steps, ...
            'total_seconds_per_step',summary.accumulatedWallSeconds/run.steps, ...
            'maximum_stage_cfl',run.maximumStageCfl,'configured_cfl',run.configuredCfl, ...
            'retries',run.retries,'checkpoint_sha256',run.acceptedCheckpoint.sha256));
    end
end
end

function [files,settings]=make_plot(folder,studies,axesData,D,qDifference,performance)
fig=figure('Visible','off','Color','white','Units','inches','Position',[1,1,7.2,6.15],'InvertHardcopy','off');
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors=[0,114,178;213,94,0;0,158,115;204,121,167]/255;
methodColors=[.45,.45,.45;0,114/255,178/255];styles={'--','-'};fields={'W','q','F','Q'};
ax=gobjects(1,4);zeroValues={};ax(1)=nexttile(layout);hold(ax(1),'on');handles=gobjects(1,4);
for m=1:2
    s=studies{m};
    for f=1:4
        v=D.(s.label).(fields{f});
        for k=find(v==0),zeroValues{end+1}=sprintf('%s/%s/pair%d',s.label,fields{f},k);end %#ok<AGROW>
        v(v==0)=NaN;face=colors(f,:);if m==1,face='white';end
        h=plot(ax(1),s.h(1:2),v,'Color',colors(f,:),'LineStyle',styles{m}, ...
            'Marker','o','MarkerSize',4,'MarkerFaceColor',face,'LineWidth',1.2);
        if m==2,handles(f)=h;end
    end
end
assert(any(arrayfun(@(h)any(isfinite(h.YData)),handles)), ...
    'ipm_cr:RedistributionPlotZero','All new self-differences are zero; a log convergence panel is undefined.');
set(ax(1),'XScale','log','YScale','log','XTick',sort(studies{1}.h(1:2)));
xlim(ax(1),[min(studies{1}.h(1:2))/1.05,max(studies{1}.h(1:2))*1.05]);
xlabel(ax(1),'Coarser h in each pair');ylabel(ax(1),'max_{cloud} |f_h - f_{h/2}|');
legend(ax(1),handles,fields,'Location','northwest','NumColumns',2,'Box','off','FontSize',8);
text(ax(1),.97,.04,'dashed: original; solid: redistributed','Units','normalized','HorizontalAlignment','right','FontSize',6.8);
if ~isempty(zeroValues),text(ax(1),.03,.14,sprintf('%d zero values: see provenance',numel(zeroValues)),'Units','normalized','FontSize',7);end
ax(2)=nexttile(layout);hold(ax(2),'on');levelColors=[230,159,0;0,158,115;0,114,178]/255;lh=gobjects(1,3);
for j=1:3
    for m=1:2
        if m==1,r=axesData(j).base;else,r=axesData(j).new;end
        h=plot(ax(2),(r(1:end-1)+r(2:end))/2,diff(r),'Color',levelColors(j,:), ...
            'LineStyle',styles{m},'LineWidth',1.2);
        if m==2,lh(j)=h;end
    end
end
info=axesData(1).geometry.radialRedistribution;
xline(ax(2),info.preserveRadius,':','Color',[.6,.6,.6],'HandleVisibility','off');
xline(ax(2),info.fullRefinementRadius,':','Color',[.6,.6,.6],'HandleVisibility','off');
set(ax(2),'XScale','log','YScale','log');xlim(ax(2),[.6,300]);
xlabel(ax(2),'Neighbor midpoint radius r');ylabel(ax(2),'Physical gap \Delta r');
legend(ax(2),lh,compose('h = %g',studies{1}.h),'Location','northwest','Box','off','FontSize',8);
text(ax(2),.97,.04,'dashed: original; solid: redistributed','Units','normalized','HorizontalAlignment','right','FontSize',6.8);
ax(3)=nexttile(layout);hold(ax(3),'on');r=studies{1}.r;E=abs(qDifference);positive=E(E>0);x0=.1;
if isempty(positive),e0=1e-12;else,e0=max(realmin*eps,min(1e-12,max(positive)/1000));end
transform=@(v,s)log1p(v/s)/log(10);markers={'.','+'};sh=gobjects(1,2);
for m=1:2
    sh(m)=scatter(ax(3),transform(r,x0),transform(E(:,m),e0),5,methodColors(m,:),markers{m}, ...
        'MarkerEdgeAlpha',.35,'LineWidth',.35);
    [v,k]=max(E(:,m));plot(ax(3),transform(r(k),x0),transform(v,e0),'o', ...
        'Color',methodColors(m,:),'MarkerFaceColor','none','MarkerSize',5,'HandleVisibility','off');
end
xt=[0,1,1e2,1e4,1e6,1e8,1e10];xt=xt(xt<=max(r)*(1+64*eps));
largest=max(E,[],'all');yt=0;
if largest>0,yt=[0,10.^(floor(log10(e0))+2:2:floor(log10(largest)))];end
set(ax(3),'XTick',transform(xt,x0),'XTickLabel',numeric_labels(xt), ...
    'YTick',transform(yt,e0),'YTickLabel',numeric_labels(yt));xlim(ax(3),[0,transform(max(r)*1.2,x0)]);
xlabel(ax(3),'Original query radius r (all samples)');ylabel(ax(3),'|q_{mid} - q_{old fine}|');
legend(ax(3),sh,{'Original','Redistributed'},'Location','northeast','Box','off','FontSize',7.5);
ax(4)=nexttile(layout);hold(ax(4),'on');minutes=zeros(3,2);counts=zeros(3,2);
for m=1:2
    for j=1:3,k=(m-1)*3+j;minutes(j,m)=performance(k).elapsed_total_seconds/60;counts(j,m)=performance(k).steps;end
end
bars=bar(ax(4),1:3,minutes,'grouped');
for m=1:2
    bars(m).FaceColor=methodColors(m,:);bars(m).EdgeColor='none';
    px=bars(m).XEndPoints;py=bars(m).YEndPoints;
    text(ax(4),px(:),py(:)+.012*max(minutes,[],'all'), ...
        compose('%d',counts(:,m)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7);
end
set(ax(4),'XTick',1:3,'XTickLabel',compose('%g',studies{1}.h));
ylim(ax(4),[0,1.2*max(minutes,[],'all')]);xlabel(ax(4),'Core spacing h');ylabel(ax(4),'Measured wall time (min)');
legend(ax(4),bars,{'Original','Redistributed'},'Location','northwest','Box','off','FontSize',7.5);
text(ax(4),.97,.97,'labels: accepted steps','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
for j=1:4
    heading=title(ax(j),sprintf('(%c)','a'+j-1),'FontWeight','normal','HorizontalAlignment','left');
    heading.Units='normalized';position=heading.Position;position(1)=0;heading.Position=position;
    set(ax(j),'FontName','Times New Roman','FontSize',9,'LineWidth',.7,'TickDir','out', ...
        'Box','off','XColor',[.2,.2,.2],'YColor',[.2,.2,.2],'Color','white');
    ax(j).XLabel.FontSize=9;ax(j).YLabel.FontSize=9;ax(j).Title.FontSize=10;
end
drawnow;files=struct('png',fullfile(folder,'redistribution_convergence.png'),'pdf',fullfile(folder,'redistribution_convergence.pdf'));
fig.PaperUnits='inches';fig.PaperSize=fig.Position(3:4);fig.PaperPosition=[0,0,fig.Position(3:4)];fig.PaperPositionMode='manual';
exportgraphics(fig,files.png,'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,files.pdf,'ContentType','vector','BackgroundColor','white');
settings=struct('figureSizeInches',[7.2,6.15],'pngDpi',300,'spacingDisplayRadiusWindow',[.6,300], ...
    'qPanelDisplayTransform','log10(1+v/scale), ticks in original recorded units; no error floor', ...
    'qPanelRadiusScale',x0,'qPanelErrorScale',e0,'exactZeroLogConvergenceValues',{zeroValues});
end

function text=make_caption(studies,validation,orders,settings,axesData)
parts={};
for m=1:2
    terms={};p=orders.(studies{m}.label);
    for key={'W','q','F','Q'}
        if isempty(p.(key{1})),value='undefined (zero difference)';else,value=sprintf('%.5g',p.(key{1}));end
        terms{end+1}=sprintf('%s=%s',key{1},value); %#ok<AGROW>
    end
    parts{end+1}=[studies{m}.label,': ',strjoin(terms,', ')]; %#ok<AGROW>
end
info=axesData(1).geometry.radialRedistribution;
text=sprintf(['Radial redistribution at canonical tau=10, using the original common physical query cloud.\n\n', ...
    '(a) Recorded sampled infinity self-differences of W, q, F and Q over ALL %d original queries, ', ...
    'including the exterior through r=1e10. Each x coordinate is the coarser h in its adjacent-level pair. ', ...
    'Lines only join two recorded differences; no smoothing, fitting or error interpolation. ', ...
    'Three-level p=log(Dcoarse/Dfine)/log(2): %s. These are sampled differences, not bounds for a continuous ', ...
    'global infinity norm. Exact zero values are recorded in provenance, never replaced by a positive log floor.\n\n', ...
    '(b) Physical neighbor gaps Delta r=r[j+1]-r[j], at neighbor midpoints, displayed for r=0.6..300. ', ...
    'CSV retains EVERY full-axis interval. Nodes are recovered from saved oldIndexAtNewNodes/analyticMap ', ...
    'metadata, with baseR copied on the protected prefix and at the outer endpoint; the independent spacingAt97 ', ...
    'record is checked. No field/error metric is reconstructed. Dotted lines mark preserveRadius=%g ', ...
    'and fullRefinementRadius=%g. Physical gaps are not the mapped CFL derivative metric.\n\n', ...
    '(c) All original query samples, including exact zeros and r=0, show each middle-grid q minus the ORIGINAL ', ...
    'fine-grid q. That fine grid is a finite-resolution reference, NOT an exact solution. This difference ', ...
    'alone does not establish accuracy. Open circles mark sampled maxima. No angular/radial selection, ', ...
    'exterior deletion, alignment, fitting or reinterpolation. Display axes use log10(1+r/%.17g) and ', ...
    'log10(1+absoluteDifference/%.17g), with ticks in original units; these monotone display transforms ', ...
    'retain exact zeros without altering CSV values or norms.\n\n', ...
    '(d) Measured summary.accumulatedWallSeconds/60 across runner invocations; labels give accepted steps. ', ...
    'CSV also preserves run.seconds=summary.invocationWallSeconds and resumed. Runner timing includes ', ...
    'initialization, advancement, rejection and checkpoint/output work, but excludes MATLAB startup and outer ', ...
    'harness grid preparation/sampling. These are recorded runs, not repeated controlled timings or extrapolated ', ...
    'costs. Both reports use the same saved Poisson options, full-stage CFL 0.3 and maxDt 0.2. Timing does not ', ...
    'establish equal accuracy or guarantee a tau200 runtime.\n\n', ...
    'Original status: %s. New status: %s. The old time-half check was incomplete; the new time-half check ', ...
    'is configured %d. Qualification flags are copied unchanged to provenance.json. Validation ', ...
    'structuralChecksPassed=%d, spatialQualificationPassed=%d, passed=%d; unresolved static-order pairs ', ...
    'are not promoted to passing. This figure replaces none of the temporal, Poisson, spatial or long-trajectory ', ...
    'accuracy tests. All values are dimensionless rescaled fields/coordinates.\n\n', ...
    'Style: journal-minimal, Times New Roman, 7.2 x 6.15 inches; no overall title. No new PDE evolution ', ...
    'or Poisson solve. Inspect the exported PNG/PDF separately.\n'], ...
    studies{1}.count,strjoin(parts,'; '),info.preserveRadius,info.fullRefinementRadius, ...
    settings.qPanelRadiusScale,settings.qPanelErrorScale,studies{1}.raw.status,studies{2}.raw.status, ...
    studies{2}.raw.options.timeHalfCheck,validation.structuralChecksPassed,validation.spatialQualificationPassed,validation.passed);
end

function flags=qualification_flags(report)
flags=struct('status',report.status);
for key={'globalFourthOrderObserved','timeContributionScreenPassed','fullIVPShortWindowQualified', ...
        'rateAndScalarFourthOrderObserved','tau200Qualified','sourcesUnchanged'}
    flags.(key{1})=get_option(report,key{1},[]);
end
end

function a=struct_array(a)
if iscell(a),a=[a{:}];end
assert(isstruct(a),'ipm_cr:RedistributionPlotSchema','Expected JSON struct array.');
end

function v=finite_vector(v,n)
v=v(:);assert(isnumeric(v)&&numel(v)==n&&all(isfinite(v)), ...
    'ipm_cr:RedistributionPlotFinite','Every retained vector needs its complete finite query count.');
end

function value=log_sinh(t)
value=t+log(-expm1(-2*t))-log(2);
end

function values=numeric_labels(v)
values=string(compose('%.0e',v));values(v==0)="0";
for j=1:numel(v),if v(j)>=.1&&v(j)<=100,values(j)=string(v(j));end;end
end

function value=get_option(s,key,fallback)
if isfield(s,key),value=s.(key);else,value=fallback;end
end

function s=append(s,item)
if isempty(s),s=item;else,s(end+1)=item;end
end

function identity=file_identity(p)
assert(isfile(p),'ipm_cr:RedistributionPlotInput','Missing input: %s',p);
p=char(java.io.File(char(p)).getCanonicalPath());info=dir(p);
identity=struct('path',p,'bytes',info.bytes,'sha256',hash_file(p));
end

function verify_identity(identity)
assert(isfile(identity.path)&&strcmp(hash_file(identity.path),identity.sha256), ...
    'ipm_cr:RedistributionPlotInputChanged','Input/helper changed during plotting: %s',identity.path);
end

function value=hash_file(p)
d=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(p,'rb');
assert(fid>=0,'ipm_cr:RedistributionPlotInput','Cannot read %s',p);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
while true,b=fread(fid,1048576,'*uint8');if isempty(b),break;end;d.update(b);end
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2)',1,[]));
end

function value=hash_text(text)
d=java.security.MessageDigest.getInstance('SHA-256');d.update(unicode2native(text,'UTF-8'));
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2)',1,[]));
end

function write_csv(file,T)
% Explicit 17 significant digits preserve raw doubles independently of format.
values=strings(height(T),width(T));
for j=1:width(T)
    column=T{:,j};
    if isnumeric(column),values(:,j)=compose('%.17g',column);
    elseif islogical(column),values(:,j)=compose('%d',column);
    else
        quote=string(char(34));column=string(column);
        values(:,j)=quote+replace(column,quote,quote+quote)+quote;
    end
end
lines=[string(strjoin(T.Properties.VariableNames,','));join(values,',',2)];
write_text(file,char(join(lines,newline)+newline));
end

function write_text(file,text)
fid=fopen(file,'w');assert(fid>=0,'ipm_cr:RedistributionPlotOutput','Cannot write %s',file);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,unicode2native(text,'UTF-8'),'uint8');
end
