function export_geometry_frames()
% Read saved checkpoints only; export small MAT-v7 fields for Python analysis.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'solver'),'-begin');
maxNumCompThreads(2);
folder = fullfile(root,'runs','run_20260927_133548_651','checkpoints');
out = fullfile(root,'runs','geometry_loop_20260928');
lines = splitlines(string(fileread(fullfile(folder,'checkpoint_manifest.jsonl'))));
records = cellfun(@jsondecode,cellstr(lines(strlength(lines)>0)),'UniformOutput',false);
times = cellfun(@(r)r.canonicalTime,records);
for target = [10,15,20]
    [~,i] = min(abs(times-target));
    sourceCheckpoint = records{i}.checkpointFile;
    loaded = load(sourceCheckpoint,'checkpoint');
    state = loaded.checkpoint.payload.state;
    xFull = state.x; yFull = state.y; rho = state.rho;
    scheme = state.config.transport.spatialDiscretization;
    if strcmp(scheme,'high_order')
        width = 7;
    elseif strcmp(scheme,'sixth_order')
        width = 9;
    else
        error('Unimplemented derivative export for %s',scheme);
    end
    Dx = ipm.mesh.fdMatrix(xFull,1,width);
    rhoX = rho*Dx';
    [~,im] = max(rhoX(1,:));
    peak = ipm.evolve.quadraticPeakFunctional(rhoX(1,:),xFull,im);
    idx = im-3:im+3;
    [coef,~,scale] = polyfit(xFull(idx),rhoX(1,idx),6);
    C = polyval(polyder(polyder(coef)),(peak.x-scale(1))/scale(2))/scale(2)^2;
    a = peak.x; P = peak.value; ell = sqrt(-2*P/C);
    assert(isreal(ell) && ell>0);
    rinterp = griddedInterpolant({yFull,xFull},rho,'spline','none');
    centerRho = rinterp(0,a);
    x = (xFull-a)/ell; y = yFull/ell;
    Q = (rho-centerRho)/(P*ell); q = rhoX/P;
    tau = state.scale.canonicalTime;
    coordinateNote = 'x=(X-a)/ell,y=Y/ell,Q=(rho-rho(a,0))/(P*ell),q=rho_X/P; 3-node peak, 7-node degree-six curvature';
    save(fullfile(out,sprintf('frame_tau%d.mat',target)), ...
        'x','y','Q','q','tau','a','P','C','ell','sourceCheckpoint', ...
        'coordinateNote','scheme','-v7');
    fprintf('Export tau %.9f, %d x %d, ell %.9f\n',tau,numel(x),numel(y),ell);
    clear loaded state rho rhoX;
end
end
