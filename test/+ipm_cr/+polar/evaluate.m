function value = evaluate(p, solution, x, y)
%EVALUATE Interpolate psi and its analytic spectral derivatives.
%   Query arrays have matching sizes. Values outside the half-annulus are
%   NaN. Panel interfaces use the panel on their right, except at rOuter.
%   Derivatives are obtained from the same interpolating polynomials as psi.

if ~isequal(size(x),size(y))
    error('ipm_cr:polar:QuerySize','x and y must have matching sizes.');
end
shape = size(x); x = x(:); y = y(:);
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
radius = hypot(x,y);
wallTolerance = 128*eps(max(1,radius));
innerTolerance = 128*eps(p.rInner);
outerTolerance = 128*eps(p.rOuter);
valid = isfinite(x) & isfinite(y) & y >= -wallTolerance ...
    & radius >= p.rInner-innerTolerance ...
    & radius <= p.rOuter+outerTolerance;
theta = atan2(max(y,0),x);
theta = min(pi,max(0,theta));
t = log(max(p.rInner,min(p.rOuter,radius)));
stream = nan(size(x)); dt = stream; dtheta = stream;
dtt = stream; dtheta2 = stream;
for j = 1:p.nPanels
    if j == p.nPanels
        selected = valid & t >= p.logEdges(j) & t <= p.logEdges(j+1);
    else
        selected = valid & t >= p.logEdges(j) & t < p.logEdges(j+1);
    end
    if ~any(selected), continue; end
    Br = ipm_cr.internal.bary(p.panelNodes{j},t(selected));
    Ba = ipm_cr.internal.bary(p.theta,theta(selected));
    local = psi(p.panelIndices{j},:);
    interpolated = Br*local;
    localDt = p.panelD{j}*local;
    stream(selected) = sum(interpolated.*Ba,2);
    dt(selected) = sum((Br*localDt).*Ba,2);
    dtheta(selected) = sum(interpolated.*(Ba*p.Dtheta),2);
    dtt(selected) = sum((Br*(p.panelD{j}*localDt)).*Ba,2);
    dtheta2(selected) = sum(interpolated.*(Ba*p.Dtheta2),2);
end
psiX = (cos(theta).*dt-sin(theta).*dtheta)./radius;
psiY = (sin(theta).*dt+cos(theta).*dtheta)./radius;
% Homogeneous wall data imply an exactly zero tangential derivative.
psiX(valid & y == 0) = 0;
value = struct('psi',reshape(stream,shape), ...
    'u1',reshape(-psiY,shape),'u2',reshape(psiX,shape), ...
    'laplacian',reshape((dtt+dtheta2)./radius.^2,shape), ...
    'dt',reshape(dt,shape),'dtheta',reshape(dtheta,shape));
end
