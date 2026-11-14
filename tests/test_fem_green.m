function report=test_fem_green(outputDir)
%TEST_FEM_GREEN Check full native sources against independent reflected Green data.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=tempname;end
file=fullfile(root,'tests','reference','fem_green_reference.json');
fixture=jsondecode(fileread(file));data=fixture.green;
assert(strcmp(fixture.format,'fem_independent_green_reference_v1')&&data.radius==1);
x=reshape(data.x,1,[]);y=reshape(data.y,1,[]);
gates=struct('weak',1e-10,'trace',1e-9,'wall',5e-12);
levels=[24,48,80];errors=zeros(numel(levels),3,2);
for k=1:numel(levels)
    n=levels(k);h=.4/n;r=exp(linspace(log(.12),0,3*n+1)).';
    axes=struct('x',(-n:n)*h,'y',(0:n)'*h,'r',r, ...
        'theta',linspace(0,pi,2*n+1),'dxMetric',h*ones(1,2*n+1), ...
        'dyMetric',h*ones(n+1,1),'drMetric',r*log(1/.12)/(3*n));
    g=ipm_cr.mesh.build(struct('axes',axes,'staircase',true,'interpolationWidth',7));
    solver=ipm_cr.space.poisson(g);
    core=hypot(x,y)<g.polar.innerActiveRadius;
    for j=1:2
        power=j-1;index=find(all(data.powers==[power,0],2),1);
        q=struct('core',g.core.X.^power,'polar',g.polar.X.^power);
        [flow,info]=solver.solve(q);
        assert(info.sourceQuadratureCount==0&&info.weakBackwardError<gates.weak&& ...
            info.traceAbsolute<gates.trace,'ipm_cr:GreenWeakRows');
        assert(max(abs(flow.core.psi(1,:)))<gates.wall&& ...
            max(abs(flow.polar.u2(:,[1,end])),[],'all')<gates.wall, ...
            'ipm_cr:GreenWall');
        assert(all(q.polar(end,[1,end])~=0),'ipm_cr:GreenCornerSource');
        vc=solver.evaluate(flow,x(core),y(core),'core');
        vp=solver.evaluate(flow,x(~core),y(~core),'polar');
        for i=1:3
            names={'psi','u1','u2'};name=names{i};val=zeros(size(x));
            val(core)=vc.(name);val(~core)=vp.(name);
            errors(k,i,j)=max(abs(val-data.(name)(index,:)));
        end
        if j==1
            [direct,~]=solver.solve(@(X,Y)ones(size(X)));
            assert(max(abs(direct.core.psi-flow.core.psi),[],'all')<2e-12&& ...
                max(abs(direct.polar.psi-flow.polar.psi),[],'all')<2e-12, ...
                'ipm_cr:GreenCachedLoad');
        end
    end
    fprintf('Green n=%d %s\n',n,sprintf(' %.4g',errors(k,:,:)));
end
% The velocity at the source cutoff corners need not have smooth-data order.
limits=[1e-9,2e-5,1e-8];fine=reshape(errors(end,:,:),3,[]);
pass=all(fine<limits','all')&&all(fine<reshape(errors(1,:,:),3,[]),'all');
report=struct('passed',pass,'intervals',levels,'errors',errors, ...
    'columns',{{'psi','u1','u2'}},'sources',{{'1','x'}},'finalLimits',limits, ...
    'fixtureSha256',ipm_cr.io.identity('sha',file),'exteriorModel','zero_source_outside_unit_disk', ...
    'wallCornerSourceNonzero',true,'timeAdvanced',0);
if ~isfolder(outputDir),mkdir(outputDir);end
ipm_cr.io.identity('json',fullfile(outputDir,'fem_green.json'),report);
assert(report.passed,'ipm_cr:GreenReference','Independent full-source Green check failed.');
end
