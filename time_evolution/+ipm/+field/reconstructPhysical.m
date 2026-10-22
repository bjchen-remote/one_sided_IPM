function physical = reconstructPhysical(rho,flow,scale,ops)
%IPM.FIELD.RECONSTRUCTPHYSICAL Convert one rescaled state to physical variables.

isAnisotropic = strcmp(ops.dynamicScaleGeometry,'anisotropic');
Cx = exp(scale.logC_l);
if isAnisotropic
    Cy = exp(scale.logC_y);
    aspect = exp(scale.logC_y-scale.logC_l);
else
    Cy = Cx;
    aspect = 1;
end

amplitudeInverse = exp(-scale.logC_omega);
rhoY = ops.Dy*rho;
if isAnisotropic
    omega = Cx*amplitudeInverse*flow.source;
    physicalRhoY = Cy*amplitudeInverse*rhoY;
    u2 = amplitudeInverse*flow.u2/aspect;
    psi = amplitudeInverse*flow.psi/Cy;
else
    omega = exp(scale.logC_l-scale.logC_omega)*flow.source;
    physicalRhoY = exp(scale.logC_l-scale.logC_omega)*rhoY;
    u2 = amplitudeInverse*flow.u2;
    psi = exp(-scale.logC_omega-scale.logC_l)*flow.psi;
end

physical = struct( ...
    'rho',amplitudeInverse*rho, ...
    'omega',omega,'rhoX',omega,'rhoY',physicalRhoY, ...
    'psi',psi,'u1',amplitudeInverse*flow.u1,'u2',u2, ...
    'x',(ops.x-scale.X_shift)/Cx,'y',ops.y/Cy, ...
    'Cx',Cx,'Cy',Cy,'Lx',1/Cx,'Ly',1/Cy,'aspect',aspect);
end
