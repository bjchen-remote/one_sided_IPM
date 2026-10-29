function value = evaluate(p, solution, x, y, requested, map)
%EVALUATE Interpolate psi and its analytic spectral derivatives.
%   Query arrays have matching sizes. Values outside the half-annulus are
%   NaN. Panel interfaces use the panel on their right, except at rOuter.
%   Derivatives are obtained from the same interpolating polynomials as psi.

if ~isequal(size(x),size(y))
    error('ipm_cr:polar:QuerySize','x and y must have matching sizes.');
end
if nargin<5||isempty(requested),requested={'psi','u1','u2','laplacian','dt','dtheta'};end
if ischar(requested),requested={requested};end
assert(iscellstr(requested)&&all(ismember(requested, ...
    {'psi','u1','u2','laplacian','dt','dtheta'})), ...
    'ipm_cr:polar:Fields','Unknown requested evaluation field.');
if nargin<6||isempty(map),map=ipm_cr.poisson.polar.prepare(p,x,y,requested);end
shape = size(x); x = x(:); y = y(:);
assert(isequaln(map.x,x)&&isequaln(map.y,y)&&isequal(map.shape,shape)&& ...
    isequal(map.thetaNodes,p.theta)&&isequal(map.logEdges,p.logEdges)&&map.radialDegree==p.degree, ...
    'ipm_cr:polar:EvaluationMap','The evaluation map must match these queries and geometry.');
if isfield(solution,'psi')
    psi = solution.psi;
elseif isfield(solution,'modal')
    psi = zeros(p.nr,p.nt);
    psi(:,p.interior) = solution.modal*p.V.';
else
    error('ipm_cr:polar:Solution','solution requires psi or modal data.');
end
if ~isequal(size(psi),[p.nr,p.nt])
    error('ipm_cr:polar:SolutionSize','solution.psi has the wrong size.');
end
radius=map.radius;theta=map.theta;valid=map.valid;
needVelocity=any(ismember(requested,{'u1','u2'}));
needDt=needVelocity||ismember('dt',requested);
needDa=needVelocity||ismember('dtheta',requested);
needLap=ismember('laplacian',requested);
stream = nan(size(x)); dt = stream; dtheta = stream;
dtt = stream; dtheta2 = stream;
for j = 1:p.nPanels
    entry=map.panels{j};selected=entry.indices;
    if isempty(selected), continue; end
    Br=entry.Br;Ba=entry.Ba;
    local = psi(p.panelIndices{j},:);
    interpolated = Br*local;
    if ismember('psi',requested),stream(selected)=sum(interpolated.*Ba,2);end
    if needDt||needLap,localDt=p.panelD{j}*local;end
    if needDt,dt(selected)=sum((Br*localDt).*Ba,2);end
    if needDa
        if isfield(entry,'angularFirst'),Da=entry.angularFirst;else,Da=Ba*p.Dtheta;end
        dtheta(selected)=sum(interpolated.*Da,2);
    end
    if needLap
        if isfield(entry,'angularSecond'),Daa=entry.angularSecond;else,Daa=Ba*p.Dtheta2;end
        dtt(selected)=sum((Br*(p.panelD{j}*localDt)).*Ba,2);
        dtheta2(selected)=sum(interpolated.*Daa,2);
    end
end
value=struct();
if ismember('psi',requested),value.psi=reshape(stream,shape);end
if ismember('u1',requested)
    value.u1=reshape(-(sin(theta).*dt+cos(theta).*dtheta)./radius,shape);
end
if ismember('u2',requested)
    psiX=(cos(theta).*dt-sin(theta).*dtheta)./radius;
    % Homogeneous wall data imply an exactly zero tangential derivative.
    psiX(valid&y==0)=0;value.u2=reshape(psiX,shape);
end
if needLap,value.laplacian=reshape((dtt+dtheta2)./radius.^2,shape);end
if ismember('dt',requested),value.dt=reshape(dt,shape);end
if ismember('dtheta',requested),value.dtheta=reshape(dtheta,shape);end
end
