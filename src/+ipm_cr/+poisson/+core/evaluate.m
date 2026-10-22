function v=evaluate(c,s,x,y)
%EVALUATE Values and derivatives of the same tensor interpolant.
Bx=ipm_cr.internal.bary(c.x,x);By=ipm_cr.internal.bary(c.y,y);
contract=@(A)sum((By*A).*Bx,2);
v=struct('psi',contract(s.psi),'u1',-contract(c.Dy*s.psi), ...
    'u2',contract(s.psi*c.Dx'), ...
    'laplacian',contract(c.Dyy*s.psi+s.psi*c.Dxx'));
end
