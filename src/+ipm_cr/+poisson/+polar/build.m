function p = build(cfg)
%BUILD Factor a log-polar spectral Poisson subdomain once.
%   The radial panels share endpoint values. At each shared endpoint a
%   derivative-continuity row replaces the differential equation. The walls
%   theta=0,pi have homogeneous Dirichlet data; their source values need not
%   vanish. Angular eigenmodes are those of the collocation operator itself.

if nargin == 0, cfg = struct(); end
degree = option(cfg, 'polarDegree', 16);
angularDegree = option(cfg, 'angularDegree', 32);
rInner = option(cfg, 'rInner', 0.7);
rOuter = option(cfg, 'rOuter', 16);
panelWidth = option(cfg, 'logPanelWidth', log(2));
outerBoundary = char(option(cfg, 'outerBoundary', 'decay'));
validateattributes(degree, {'numeric'}, {'scalar','integer','>=',3});
validateattributes(angularDegree, {'numeric'}, {'scalar','integer','>=',3});
validateattributes(rInner, {'numeric'}, {'scalar','real','positive'});
validateattributes(rOuter, {'numeric'}, {'scalar','real','>',rInner});
validateattributes(panelWidth, {'numeric'}, {'scalar','real','positive'});
if ~ismember(outerBoundary, {'decay','dirichlet'})
    error('ipm_cr:polar:OuterBoundary', ...
        'outerBoundary must be ''decay'' or ''dirichlet''.');
end

nPanels = max(1, ceil(log(rOuter/rInner)/panelWidth));
edges = linspace(log(rInner), log(rOuter), nPanels+1);
nr = nPanels*degree+1;
t = zeros(nr,1);
Dt = zeros(nr);
Dtt = zeros(nr);
panelNodes = cell(nPanels,1);
panelD = cell(nPanels,1);
panelIndices = cell(nPanels,1);
for j = 1:nPanels
    indices = (j-1)*degree+(1:degree+1);
    [tj,Dj] = ipm_cr.internal.cheb(degree, edges(j), edges(j+1));
    tj = tj(:);
    t(indices) = tj;
    panelNodes{j} = tj;
    panelD{j} = Dj;
    panelIndices{j} = indices;
    Dt(indices,indices) = Dt(indices,indices) + Dj;
    Dtt(indices,indices) = Dtt(indices,indices) + Dj*Dj;
end
interfaceRows = (degree+1:degree:nr-degree).';
% The averaged derivative is useful for diagnostics and nodal evaluation.
Dt(interfaceRows,:) = 0.5*Dt(interfaceRows,:);
Dtt(interfaceRows,:) = 0.5*Dtt(interfaceRows,:);
flux = zeros(nr);
for j = 1:nPanels-1
    row = j*degree+1;
    il = panelIndices{j}; ir = panelIndices{j+1};
    flux(row,il) = flux(row,il) + panelD{j}(end,:);
    flux(row,ir) = flux(row,ir) - panelD{j+1}(1,:);
end
equationRows = setdiff((2:nr-1).', interfaceRows);

[theta,Dtheta] = ipm_cr.internal.cheb(angularDegree, 0, pi);
theta = theta(:).';
interior = 2:angularDegree;
Dtheta2 = Dtheta*Dtheta;
[V,lambda] = eig(-Dtheta2(interior,interior));
mu = diag(lambda);
if max(abs(imag(mu))) > 1e-9*max(1,max(abs(mu)))
    error('ipm_cr:polar:AngularSpectrum', ...
        'The angular Dirichlet spectrum is not real.');
end
[mu,order] = sort(real(mu));
V = real(V(:,order));
if any(mu <= 0)
    error('ipm_cr:polar:AngularSpectrum', ...
        'The angular Dirichlet spectrum must be positive.');
end
% Fix arbitrary eigenvector signs for reproducible cached operators.
for k = 1:size(V,2)
    [~,pivot] = max(abs(V(:,k)));
    if V(pivot,k) < 0, V(:,k) = -V(:,k); end
end
Vinv = V\eye(size(V));
nModes = numel(mu);
factors = cell(nModes,1);
operators = cell(nModes,1);
absOperators = cell(nModes,1);
rowScale = zeros(nr,nModes);
harmonicModes = zeros(nr,nModes);
harmonicResidualInfinity = zeros(1,nModes);
harmonicBackwardError = zeros(1,nModes);
harmonicRefinementUpdates = zeros(1,nModes);
roundoffGuard = 128*eps;
for k = 1:nModes
    A = -Dtt + mu(k)*eye(nr);
    A(interfaceRows,:) = flux(interfaceRows,:);
    A(1,:) = 0; A(1,1) = 1;
    if strcmp(outerBoundary, 'dirichlet')
        A(end,:) = 0; A(end,end) = 1;
    else
        A(end,:) = Dt(end,:);
        A(end,end) = A(end,end) + sqrt(mu(k));
    end
    scale = 1./max(abs(A),[],2);
    operators{k} = sparse(A);
    absOperators{k} = abs(operators{k});
    factors{k} = decomposition(spdiags(scale,0,nr,nr)*operators{k}, 'lu');
    rowScale(:,k) = scale;
    rhs = zeros(nr,1); rhs(1) = 1;
    harmonic = factors{k}\(scale.*rhs);
    harmonic(1) = 1;
    % Refine the cached unit-inner response against the ORIGINAL rows,
    % exactly as solve does for its source/outer-boundary right hand side.
    % Subsequent trace updates are linear combinations of this response.
    residual = rhs-operators{k}*harmonic;
    backward = max(abs(residual)./max(absOperators{k}*abs(harmonic)+abs(rhs),realmin));
    for attempt = 1:3
        if backward <= roundoffGuard, break; end
        candidate = harmonic+factors{k}\(scale.*residual);
        candidate(1) = 1;
        candidateResidual = rhs-operators{k}*candidate;
        candidateBackward = max(abs(candidateResidual)./ ...
            max(absOperators{k}*abs(candidate)+abs(rhs),realmin));
        if ~isfinite(candidateBackward) || candidateBackward >= backward, break; end
        harmonic = candidate;
        residual = candidateResidual;
        backward = candidateBackward;
        harmonicRefinementUpdates(k) = harmonicRefinementUpdates(k)+1;
    end
    harmonicModes(:,k) = harmonic;
    harmonicResidualInfinity(k) = norm(residual,inf);
    harmonicBackwardError(k) = backward;
end
r = exp(t);
r(1) = rInner; r(end) = rOuter;
p = struct('degree',degree,'angularDegree',angularDegree, ...
    'rInner',rInner,'rOuter',rOuter,'outerBoundary',outerBoundary, ...
    'nPanels',nPanels,'nr',nr,'nt',angularDegree+1, ...
    't',t,'r',r,'theta',theta,'logEdges',edges, ...
    'panelNodes',{panelNodes},'panelD',{panelD}, ...
    'panelIndices',{panelIndices},'Dt',Dt,'Dtt',Dtt, ...
    'Dtheta',Dtheta,'Dtheta2',Dtheta2,'interior',interior, ...
    'V',V,'Vinv',Vinv,'mu',mu,'factors',{factors}, ...
    'operators',{operators},'absOperators',{absOperators},'flux',flux, ...
    'rowScale',rowScale,'equationRows',equationRows, ...
    'interfaceRows',interfaceRows,'harmonicModes',harmonicModes, ...
    'harmonicResidualInfinity',harmonicResidualInfinity, ...
    'harmonicBackwardError',harmonicBackwardError, ...
    'harmonicRefinementRoundoffGuard',roundoffGuard, ...
    'harmonicRefinementUpdates',harmonicRefinementUpdates);
p.X = r*cos(theta);
p.Y = r*sin(theta);
p.Y(:,[1,end]) = 0;
end

function value = option(cfg, name, fallback)
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = fallback;
end
end
