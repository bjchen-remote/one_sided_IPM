function geometry = weno5Geometry(axis)
%IPM.FIELD.WENO5GEOMETRY Precompute fifth-order nonuniform face interpolation.
%   Values retain the solver-wide nodal-value meaning. The grid-dependent
%   Lagrange and divided-difference coefficients are rebuilt only when an
%   adaptive axis changes.

axis = axis(:)';
numberOfNodes = numel(axis);
if numberOfNodes < 7
    error('ipm:WenoGridSize', ...
        'weno5_nonuniform requires at least seven nodes per axis.');
end
faces = (3:numberOfNodes-3)';
geometry = struct('faces',faces, ...
    'left',side_geometry(axis,faces,-2:2), ...
    'right',side_geometry(axis,faces,-1:3));
end

function side = side_geometry(axis,faces,offsets)
numberOfFaces = numel(faces);
indices = faces+offsets;
fullWeights = zeros(numberOfFaces,5);
candidateWeights = zeros(numberOfFaces,3,5);
differenceWeights = zeros(numberOfFaces,5);
deltaCoordinates = zeros(numberOfFaces,4);
for faceNumber = 1:numberOfFaces
    face = faces(faceNumber);
    target = 0.5*(axis(face)+axis(face+1));
    nodes = axis(indices(faceNumber,:));
    scale = max(abs(nodes-target));
    coordinates = (nodes-target)/scale;
    fullWeights(faceNumber,:) = lagrange_weights(coordinates,0);
    for candidate = 1:3
        local = candidate:candidate+2;
        candidateWeights(faceNumber,candidate,local) = ...
            lagrange_weights(coordinates(local),0);
    end
    for node = 1:5
        others = [1:node-1,node+1:5];
        differenceWeights(faceNumber,node) = ...
            factorial(4)/prod(coordinates(node)-coordinates(others));
    end
    deltaCoordinates(faceNumber,:) = diff(coordinates);
end
side = struct('indices',indices,'fullWeights',fullWeights, ...
    'candidateWeights',candidateWeights, ...
    'differenceWeights',differenceWeights, ...
    'deltaCoordinates',deltaCoordinates);
end

function weights = lagrange_weights(nodes,target)
numberOfNodes = numel(nodes);
weights = ones(1,numberOfNodes);
for node = 1:numberOfNodes
    others = [1:node-1,node+1:numberOfNodes];
    weights(node) = prod((target-nodes(others)) ./ ...
        (nodes(node)-nodes(others)));
end
end
