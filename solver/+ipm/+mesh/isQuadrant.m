function yes = isQuadrant(ops)
%IPM.MESH.ISQUADRANT Explicit first-quadrant runtime geometry marker.
yes = isfield(ops,'quadrantOnly') && ops.quadrantOnly;
end
