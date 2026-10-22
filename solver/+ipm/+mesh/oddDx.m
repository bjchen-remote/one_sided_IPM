function D = oddDx(ops)
%IPM.MESH.ODDDX X derivative for a field odd across the omitted x=0 side.
if ipm.mesh.isQuadrant(ops)
    D = ops.DxOdd;
else
    D = ops.Dx;
end
end
