function report=diagnose_ipm_cr_poisson_backward(outputDir,options)
%DIAGNOSE_IPM_CR_POISSON_BACKWARD Read-only original-system residual audit.
% requireAccuracy=false is used ONLY to retain the failed solution for this
% diagnostic. No solver, gate, RHS, boundary condition or production test is
% changed. The analytic source is sampled exactly as in the fourth-order
% Poisson adapter's zero-residual analytic-offset path.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','poisson_backward_diagnostic');end
o=resolve(options);outputDir=char(outputDir);assert(startsWith(outputDir,filesep),'ipm_cr:BackwardOutput');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'poisson_backward.json')),'ipm_cr:BackwardPreserve','Choose a fresh output directory.');
timer=tic;report=struct('kind','original_discrete_Poisson_componentwise_backward_diagnostic_v1', ...
    'status','running','options',o,'requireAccuracy',false,'diagnosticOnly',true, ...
    'originalAccuracyThresholdChanged',false,'qualificationPassedClaimed',false, ...
    'source','q=4*x*y*(3-x^2-y^2)*exp(-x^2-y^2)', ...
    'sourceInventory',source_inventory(root),'sourceHashesUnchanged',false,'elapsedSeconds',0, ...
    'interpretation',{{'Backward errors are recomputed from ORIGINAL unscaled operators, stored modal solutions and original RHS.', ...
        'A small normwise residual need not imply a small componentwise residual at an exponentially tiny row.', ...
        'Condition estimates concern forward sensitivity; ill conditioning alone is not proof of a large backward error.', ...
        'Power-scaled/Kahan dot diagnostics use double mantissa products; they avoid subnormal products but are not arbitrary-precision arithmetic.'}});
persist(outputDir,report);
try
    solver=ipm_cr.psolver(o.solverOptions);c=solver.Core;p=solver.Polar;
    qc=source(c.X,c.Y);qp=source(p.X,p.Y);
    [s,info]=solver.solve(struct('core',qc,'polar',qp));
    report.solverAudit=info;report.config=solver.Config;
    report.geometry=struct('coreShape',size(c.X),'polarShape',size(p.X), ...
        'nModes',numel(p.mu),'nPanels',p.nPanels,'rOuter',p.rOuter, ...
        'coreEigenvectorCondition',c.eigenvectorCondition,'schurCondition',solver.Coupling.condition);
    modes=struct([]);perRow=struct('ratio',zeros(p.nr,1),'mode',zeros(p.nr,1), ...
        'residual',zeros(p.nr,1),'denominator',zeros(p.nr,1),'rhs',zeros(p.nr,1),'solution',zeros(p.nr,1));
    allRatios=zeros(p.nr,numel(p.mu));
    for k=1:numel(p.mu)
        A=p.operators{k};u=s.polar.modal(:,k);b=s.polar.rhs(:,k);
        residual=A*u-b;productMagnitude=abs(A)*abs(u);
        denominator=productMagnitude+abs(b);ratios=abs(residual)./max(denominator,realmin);
        allRatios(:,k)=ratios;[maximum,row]=max(ratios);
        rowDetail=examine_row(A,u,b,p.rowScale(:,k),row,p,k);
        [~,order]=sort(ratios,'descend');top=order(1:min(o.topRows,numel(order)));
        worstRows=struct([]);
        for j=1:numel(top)
            i=top(j);worstRows=append(worstRows,summary_row(i,k,residual(i),denominator(i),u(i),b(i),p));
        end
        item=struct('mode',k,'mu',p.mu(k),'maximumComponentwiseBackward',maximum, ...
            'worstRow',row,'worstRowDetail',rowDetail,'worstRows',worstRows, ...
            'rhsInfinity',norm(b,inf),'solutionInfinity',norm(u,inf),'residualInfinity',norm(residual,inf), ...
            'operatorInfinity',norm(A,inf),'rowScaledOperatorInfinity',norm(spdiags(p.rowScale(:,k),0,p.nr,p.nr)*A,inf), ...
            'normwiseBackward',norm(residual,inf)/max(norm(A,inf)*norm(u,inf)+norm(b,inf),realmin), ...
            'subnormalSolutionCount',nnz(abs(u)>0&abs(u)<realmin), ...
            'subnormalRhsCount',nnz(abs(b)>0&abs(b)<realmin), ...
            'denominatorBelowRealminCount',nnz(denominator>0&denominator<realmin), ...
            'scaledNonzeroRhsLostToZeroCount',nnz(b~=0&p.rowScale(:,k).*b==0));
        modes=append(modes,item);
        use=ratios>perRow.ratio;
        perRow.ratio(use)=ratios(use);perRow.mode(use)=k;perRow.residual(use)=residual(use);
        perRow.denominator(use)=denominator(use);perRow.rhs(use)=b(use);perRow.solution(use)=u(use);
    end
    [~,modeOrder]=sort([modes.maximumComponentwiseBackward],'descend');
    conditionRows=struct([]);
    for j=1:min(o.conditionModes,numel(modeOrder))
        k=modeOrder(j);A=p.operators{k};scaled=spdiags(p.rowScale(:,k),0,p.nr,p.nr)*A;
        conditionRows=append(conditionRows,struct('mode',k,'mu',p.mu(k), ...
            'unscaledCondestOneNorm',condest(A),'rowScaledCondestOneNorm',condest(scaled)));
    end
    [maximum,index]=max(allRatios(:));[worstRow,worstMode]=ind2sub(size(allRatios),index);
    report.polar=struct('perMode',modes,'perRadialRowWorstMode',perRow, ...
        'maximumComponentwiseBackward',maximum,'worstMode',worstMode,'worstRow',worstRow, ...
        'worstDetail',modes(worstMode).worstRowDetail,'conditionEstimates',conditionRows);
    report.core=core_audit(c,s.core.psi,qc,o.topRows);
    report.recomputedGlobalBackward=max(maximum,report.core.maximumComponentwiseBackward);
    report.recomputedMinusSolverReported=report.recomputedGlobalBackward-info.linearBackwardError;
    report.independentManufactureProbes=manufacture_probe(solver,s);
    diagnosticData=struct('polarOperators',{p.operators},'polarModal',s.polar.modal, ...
        'polarRhs',s.polar.rhs,'polarRowScale',p.rowScale,'polarRatios',allRatios, ...
        'polarR',p.r,'polarMu',p.mu,'polarEquationRows',p.equationRows, ...
        'polarInterfaceRows',p.interfaceRows,'corePsi',s.core.psi,'coreSource',qc, ...
        'coreDxx',c.Dxx,'coreDyy',c.Dyy,'polarSource',qp); %#ok<NASGU>
    save(fullfile(outputDir,'original_system_data.mat'),'diagnosticData','-v7.3');
    report.originalSystemsMat=fullfile(outputDir,'original_system_data.mat');
    report.sourceHashesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourceHashesUnchanged,'ipm_cr:BackwardSources','Numerical sources changed during diagnosis.');
    report.status='completed_failed_or_accepted_solution_diagnostic';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
fprintf('POISSON BACKWARD diagnostic %.9g mode=%d row=%d core=%.3e solver=%.9g %.3fs\n', ...
    report.polar.maximumComponentwiseBackward,report.polar.worstMode,report.polar.worstRow, ...
    report.core.maximumComponentwiseBackward,report.solverAudit.linearBackwardError,report.elapsedSeconds);
end

function o=resolve(options)
o=struct('solverOptions',struct('coreDegree',[64,64],'angularDegree',80, ...
    'polarDegree',20,'rOuter',1e10,'requireAccuracy',false),'topRows',8,'conditionModes',3);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(o))),'ipm_cr:BackwardOptions');
names=fieldnames(options);for j=1:numel(names),o.(names{j})=options.(names{j});end
o.solverOptions.requireAccuracy=false;
validateattributes(o.topRows,{'numeric'},{'scalar','integer','>=',1});
validateattributes(o.conditionModes,{'numeric'},{'scalar','integer','>=',0});
end

function info=examine_row(A,u,b,scale,row,p,mode)
[~,columns,coefficients]=find(A(row,:));columns=columns(:);coefficients=coefficients(:);
values=u(columns);products=coefficients.*values;Au=A*u;absAu=abs(A)*abs(u);
r=Au(row)-b(row);denominator=absAu(row)+abs(b(row));
count=numel(columns);gamma=(count+2)*eps/(1-(count+2)*eps);
quantizationBound=abs(coefficients)'*eps(abs(values));
scaledA=spdiags(scale,0,p.nr,p.nr)*A;scaledB=scale.*b;
scaledR=scaledA*u-scaledB;scaledDen=abs(scaledA)*abs(u)+abs(scaledB);
power=power_dot(coefficients,values,b(row));
info=summary_row(row,mode,r,denominator,u(row),b(row),p);
info.rowScale=scale(row);info.coefficientInfinity=max(abs(coefficients));
info.rowNonzeros=count;info.referencedSolutionInfinity=max(abs(values));
nz=abs(values(values~=0));info.minimumReferencedNonzeroSolution=0;if ~isempty(nz),info.minimumReferencedNonzeroSolution=min(nz);end
info.referencedSubnormalSolutionCount=nnz(abs(values)>0&abs(values)<realmin);
info.subnormalProductsCount=nnz(abs(products)>0&abs(products)<realmin);
info.nonzeroFactorsWithUnderflowedZeroProductCount=nnz(coefficients~=0&values~=0&products==0);
info.ordinaryDotRoundoffGamma=gamma;
info.normalArithmeticRoundoffMagnitude=gamma*denominator;
info.normalArithmeticBoundApplicable=all(products==0|abs(products)>=realmin)&& ...
    all(values==0|abs(values)>=realmin)&&denominator>=realmin;
info.residualOverNormalArithmeticBound=abs(r)/max(gamma*denominator,realmin);
info.oneUlpStoredSolutionQuantizationMagnitude=quantizationBound;
info.quantizationMagnitudeOverAuditDenominator=quantizationBound/max(denominator,realmin);
info.scaledRhs=scaledB(row);info.scaledResidual=scaledR(row);info.scaledDenominator=scaledDen(row);
info.scaledComponentwiseBackward=abs(scaledR(row))/max(scaledDen(row),realmin);
info.scaledNonzeroRhsLostToZero=b(row)~=0&&scaledB(row)==0;
info.powerNormalizedProductSum=power;
info.stencilColumns=columns;info.stencilCoefficients=coefficients;
info.stencilSolution=values;info.stencilProducts=products;
info.containsSubnormalOrUnderflow=info.referencedSubnormalSolutionCount>0|| ...
    info.subnormalProductsCount>0||info.nonzeroFactorsWithUnderflowedZeroProductCount>0|| ...
    (denominator>0&&denominator<realmin)||info.scaledNonzeroRhsLostToZero;
end

function info=summary_row(row,mode,residual,denominator,solution,rhs,p)
type='PDE';if row==1,type='inner_trace';elseif row==p.nr,type='outer_boundary';
elseif ismember(row,p.interfaceRows),type='panel_flux';end
info=struct('mode',mode,'row',row,'rowType',type,'r',p.r(row),'logR',p.t(row), ...
    'residual',residual,'denominator',denominator,'auditClampedDenominator',max(denominator,realmin), ...
    'componentwiseBackward',abs(residual)/max(denominator,realmin), ...
    'solutionAtRow',solution,'rhsAtRow',rhs, ...
    'solutionSubnormal',abs(solution)>0&&abs(solution)<realmin, ...
    'rhsSubnormal',abs(rhs)>0&&abs(rhs)<realmin, ...
    'residualSubnormal',abs(residual)>0&&abs(residual)<realmin, ...
    'denominatorBelowRealmin',denominator>0&&denominator<realmin);
end

function info=power_dot(a,u,b)
% Form each product from normal mantissas and combine at a common power.
% This retains products whose direct a*u would round into subnormals.
nonzero=a~=0&u~=0;a=a(nonzero);u=u(nonzero);
[fa,ea]=log2(abs(a));[fu,eu]=log2(abs(u));
mantissa=sign(a).*sign(u).*fa.*fu;exponent=ea+eu;
if b~=0,[fb,eb]=log2(abs(b));mantissa=[mantissa;-sign(b)*fb];exponent=[exponent;eb];end
if isempty(mantissa)
    info=struct('commonPowerOfTwo',0,'normalizedResidual',0,'normalizedAbsoluteSum',0, ...
        'residualRoundedToDouble',0,'unclampedBackward',0,'auditClampedBackward',0);return;
end
power=max(exponent);terms=pow2(mantissa,exponent-power);total=0;compensation=0;
for j=1:numel(terms)
    y=terms(j)-compensation;t=total+y;compensation=(t-total)-y;total=t;
end
absoluteSum=sum(abs(terms));floorScaled=pow2(realmin,-power);
info=struct('commonPowerOfTwo',power,'normalizedResidual',total, ...
    'normalizedAbsoluteSum',absoluteSum,'residualRoundedToDouble',pow2(total,power), ...
    'unclampedBackward',abs(total)/absoluteSum, ...
    'auditClampedBackward',abs(total)/max(absoluteSum,floorScaled));
end

function info=core_audit(c,psi,q,topCount)
residual=-c.Dyy*psi-psi*c.Dxx'-q;
denominator=abs(c.Dyy)*abs(psi)+abs(psi)*abs(c.Dxx')+abs(q);
mask=false(size(psi));mask(c.iy,c.ix)=true;ratios=abs(residual)./max(denominator,realmin);
indices=find(mask);[maximum,k]=max(ratios(indices));worst=indices(k);
[~,order]=sort(ratios(indices),'descend');rows=struct([]);
for j=1:min(topCount,numel(order))
    index=indices(order(j));[iy,ix]=ind2sub(size(psi),index);
    rows=append(rows,struct('row',iy,'column',ix,'X',c.X(index),'Y',c.Y(index), ...
        'residual',residual(index),'denominator',denominator(index), ...
        'componentwiseBackward',ratios(index),'solution',psi(index),'rhsSource',q(index)));
end
info=struct('maximumComponentwiseBackward',maximum,'worstLinearIndex',worst, ...
    'worstRows',rows,'allInteriorRows',struct('linearIndices',indices, ...
        'residual',residual(indices),'denominator',denominator(indices), ...
        'ratio',ratios(indices),'solution',psi(indices),'rhsSource',q(indices)));
end

function info=manufacture_probe(solver,s)
theta=pi*((0:40)+.371)/41;r=[.1;.4;.8;1.1;1.6;2;4;10;1e3;1e10];
[R,A]=ndgrid(r,theta);X=R.*cos(A);Y=R.*sin(A);v=solver.evaluate(s,X,Y);
e=exp(-X.^2-Y.^2);exactPsi=X.*Y.*e;
u1=-X.*(1-2*Y.^2).*e;u2=Y.*(1-2*X.^2).*e;
info=struct('psiInfinity',max(abs(v.psi-exactPsi),[],'all'), ...
    'velocityInfinity',max(hypot(v.u1-u1,v.u2-u2),[],'all'), ...
    'definition','Independent exact streamfunction and Cartesian velocity at fractionally offset radial rays; finite probes only.');
end

function q=source(X,Y)
r2=X.^2+Y.^2;q=4*X.*Y.*(3-r2).*exp(-r2);
end

function inventory=source_inventory(root)
files=dir(fullfile(root,'src','+ipm_cr','**','*.m'));paths={mfilename('fullpath')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for j=1:numel(files),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths)
    f=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',f.bytes,'sha256',ipm_cr.io.sha(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'poisson_backward.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'poisson_backward.mat'),'f');assert(ok,'ipm_cr:BackwardSave','%s',message);
temporary=fullfile(folder,'poisson_backward.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:BackwardSave');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'poisson_backward.json'),'f');assert(ok,'ipm_cr:BackwardSave','%s',message);
end

function a=append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
