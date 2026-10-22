function scale = initialScale(ops)
%IPM.EVOLVE.INITIALSCALE Identity scale state for the selected geometry.

if strcmp(ops.dynamicScaleGeometry,'anisotropic')
    scale = struct('logC_l',0,'logC_x',0,'logC_y',0, ...
        'logC_omega',0,'physicalTime',0,'X_shift',0,'canonicalTime',0);
else
    scale = struct('logC_l',0,'logC_omega',0,'physicalTime',0, ...
        'X_shift',0,'canonicalTime',0);
end
end
