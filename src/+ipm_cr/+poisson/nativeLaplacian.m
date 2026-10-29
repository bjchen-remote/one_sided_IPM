function [core,polar]=nativeLaplacian(c,p,solution,plan)
%NATIVELAPLACIAN Separable fixed-query reconstruction audit, not a PDE gate.
% The same polynomial is evaluated on the audit tensor coordinates. Matrix
% regrouping, and avoiding Cartesian->polar coordinate round trips, change
% floating-point rounding relative to generic pointwise evaluation. The
% continuous-callback strong-residual acceptance path does not use this.
m=plan.nativeCoreTensor;
lap=c.Dyy*solution.core.psi+solution.core.psi*c.Dxx';
core=m.By*lap*m.Bx';core=core(:);
m=plan.nativePolarTensor;polar=nan(m.shape);
for j=1:p.nPanels
    panel=m.panels{j};U=solution.polar.psi(p.panelIndices{j},:);
    sampled=panel.Br*U;
    radialSecond=panel.Br*(p.panelD{j}*(p.panelD{j}*U));
    polar(panel.rows,:)=(radialSecond*m.angular'+sampled*m.angularSecond')./panel.r.^2;
end
end
