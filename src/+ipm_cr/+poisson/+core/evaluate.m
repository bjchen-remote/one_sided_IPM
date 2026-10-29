function v=evaluate(c,s,x,y,requested,map)
%EVALUATE Values and derivatives of the same tensor interpolant.
% REQUESTED optionally selects fields; MAP is a geometry-only prepare map.
if nargin<5||isempty(requested),requested={'psi','u1','u2','laplacian'};end
if ischar(requested),requested={requested};end
assert(iscellstr(requested)&&all(ismember(requested,{'psi','u1','u2','laplacian'})), ...
    'ipm_cr:core:Fields','Unknown requested evaluation field.');
if nargin<6||isempty(map),map=ipm_cr.poisson.core.prepare(c,x,y);end
assert(isequaln(map.x,x(:))&&isequaln(map.y,y(:))&& ...
    isequal(map.xNodes,c.x)&&isequal(map.yNodes,c.y), ...
    'ipm_cr:core:EvaluationMap','The evaluation map must match these queries and geometry.');
Bx=map.Bx;By=map.By;
contract=@(A)sum((By*A).*Bx,2);
v=struct();
if ismember('psi',requested),v.psi=contract(s.psi);end
if ismember('u1',requested),v.u1=-contract(c.Dy*s.psi);end
if ismember('u2',requested),v.u2=contract(s.psi*c.Dx');end
if ismember('laplacian',requested)
    v.laplacian=contract(c.Dyy*s.psi+s.psi*c.Dxx');
end
end
