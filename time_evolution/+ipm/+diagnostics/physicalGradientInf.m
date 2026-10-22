function value = physicalGradientInf(measurement,scale,rho,ops)
%IPM.DIAGNOSTICS.PHYSICALGRADIENTINF Physical gradient on the stored grid.

if isfield(scale,'logC_y')
    rhoX = rho*ops.Dx';
    rhoY = ops.Dy*rho;
    physicalRhoX = exp(scale.logC_l-scale.logC_omega)*rhoX;
    physicalRhoY = exp(scale.logC_y-scale.logC_omega)*rhoY;
    value = max(hypot(physicalRhoX,physicalRhoY),[],'all');
else
    value = exp(scale.logC_l-scale.logC_omega)*measurement.gradInf;
end
end
