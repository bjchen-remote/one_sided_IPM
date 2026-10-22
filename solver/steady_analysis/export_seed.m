function export_seed(resultFile,checkpointFile,outputDir)
%EXPORT_SEED Read-only analysis of a saved run and centered steady seeds.
% No time evolution, symmetry projection, or changes to the saved run.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,'-begin');
maxNumCompThreads(2);
assert(~isfolder(outputDir),'Use a new output directory.');
loaded = load(resultFile,'result');
r = loaded.result;
[state,~,~] = ipm.output.restoreCheckpoint(checkpointFile, ...
    struct('makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false));
ops = state.ops;
assert(isequal(ops.x,r.grid.x) && isequal(ops.y,r.grid.y));
assert(isequal(state.rho,r.state.rho));
mkdir(outputDir);
xFull = ops.x; yFull = ops.y;
tau = r.snapshots.canonicalTime(:);
n = numel(tau);
wallQ = zeros(n,numel(xFull));
geometry = zeros(n,4);
for k=1:n
    assert(isequal(r.snapshots.x{k},xFull) && isequal(r.snapshots.y{k},yFull), ...
        'Snapshot meshes differ; use per-frame operators.');
    wallQ(k,:) = r.snapshots.rho{k}(1,:)*ops.Dx';
    geometry(k,:) = peak_geometry(wallQ(k,:),xFull);
end
rho = state.rho; flow = state.flow;
[rhoRate,~] = ipm.evolve.flow(rho,ops,state.scale);
rhoX = rho*ops.Dx'; rhoXRate = rhoRate*ops.Dx';
psi = flow.psi; u1 = flow.u1; u2 = flow.u2;
a = geometry(end,1); P = geometry(end,2);
C = geometry(end,3); ell = geometry(end,4);
dt = tau(end)-tau(end-1);
beta = log(ell/geometry(end-1,4))/dt;
piRate = log(P/geometry(end-1,2))/dt;
aDot = (a-geometry(end-1,1))/dt;
sourceTau = tau(end); sourceSnapshotIndex = n;
sourceCL = flow.c_l; sourceCOmega = flow.c_omega; sourceCR = flow.c_r;
uAtPeak = interp1(xFull,u1(1,:),a,'spline');
v0 = (uAtPeak+sourceCL*a+sourceCR-aDot)/(P*ell);
rates = [(sourceCL-beta)/P,(sourceCOmega-beta-piRate)/P,v0,v0];
common = r.history.common;
physicalTime = r.state.physicalTime;
stopReason = r.state.stopReason;
poissonResidual = flow.poissonResidual;
metadata = struct('sourceResultFile',resultFile,'sourceCheckpoint',checkpointFile, ...
    'stopReason',stopReason,'tau',sourceTau,'physicalTime',physicalTime, ...
    'c_l',sourceCL,'c_omega',sourceCOmega,'c_r',sourceCR, ...
    'a',a,'P',P,'C',C,'ell',ell,'beta',beta,'piRate',piRate,'aDot',aDot, ...
    'uAtPeak',uAtPeak,'innerRates',rates,'poissonResidual',poissonResidual, ...
    'rhoRateInf',max(abs(rhoRate),[],'all'), ...
    'rhoXRateInf',max(abs(rhoXRate),[],'all'),'remeshCount',r.grid.remeshCount);
save(fullfile(outputDir,'run_analysis.mat'),'xFull','yFull','tau','wallQ', ...
    'geometry','rho','rhoX','rhoRate','rhoXRate','psi','u1','u2','common','metadata','-v7');
fid = fopen(fullfile(outputDir,'source_summary.json'),'w');
fprintf(fid,'%s\n',jsonencode(metadata,PrettyPrint=true)); fclose(fid);
rhoInterp = griddedInterpolant({yFull,xFull},rho,'spline','none');
% Remove a harmonic linear function to reduce cancellation of bulk drift.
psiInterp = griddedInterpolant({yFull,xFull},psi+uAtPeak*yFull(:),'spline','none');
centerRho = rhoInterp(0,a);
specs = [4,65,33;4,129,65;6,193,97];
for j=1:size(specs,1)
    L = specs(j,1); nx = specs(j,2); ny = specs(j,3);
    x = linspace(-L,L,nx); y = linspace(0,L,ny);
    [X,Y] = meshgrid(x,y);
    Q = (rhoInterp(ell*Y,a+ell*X)-centerRho)/(P*ell);
    psiQ = psiInterp(ell*Y,a+ell*X)/(P*ell^2);
    psi_boundary = zeros(ny,nx);
    psi_boundary(end,:) = psiQ(end,:);
    psi_boundary(:,[1,end]) = psiQ(:,[1,end]);
    psi_boundary(1,:) = 0;
    sourceResultFile = resultFile;
    geometryConvention = 'a,P: three-node quadratic of discrete R_X; C: seven-node degree-six fit; ell=sqrt(-2P/C)';
    boundaryConvention = 'Full saved half-plane psi plus harmonic uAtPeak*Y, normalized by P*ell^2; frozen exterior response';
    rateConvention = 'Centered density Q, ds/dtau=P; rates=(cl,cw,v0,b), b=v0; frame rates from last two outputs';
    assert(all(isfinite(Q),'all') && all(isfinite(psi_boundary),'all'));
    save(fullfile(outputDir,sprintf('seed_L%d_%dx%d.mat',L,nx,ny)), ...
        'x','y','Q','psi_boundary','rates','a','P','C','ell','centerRho', ...
        'sourceTau','sourceSnapshotIndex','sourceCL','sourceCOmega','sourceCR', ...
        'uAtPeak','beta','piRate','aDot','poissonResidual', ...
        'sourceResultFile','geometryConvention','boundaryConvention','rateConvention','-v7');
end
fprintf('tau=%.9g t=%.9g norm=%.6e C=(%.9g %.9g %.9g) ell=%.9g\n', ...
    sourceTau,physicalTime,metadata.rhoXRateInf,sourceCL,sourceCOmega,sourceCR,ell);
fprintf('inner rates=(%.9g %.9g %.9g %.9g), beta=%.3e\n',rates,beta);
end

function g = peak_geometry(q,x)
[~,i] = max(q);
p = ipm.evolve.quadraticPeakFunctional(q,x,i);
idx = i-3:i+3;
[c,~,scale] = polyfit(x(idx),q(idx),6);
C = polyval(polyder(polyder(c)),(p.x-scale(1))/scale(2))/scale(2)^2;
assert(C<0 && p.value>0);
g = [p.x,p.value,C,sqrt(-2*p.value/C)];
end
