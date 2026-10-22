function weights = quadrature6(axis)
%IPM.MESH.QUADRATURE6 Positive mapped degree-seven quadrature norm.
%   An arbitrary smooth monotone nodal axis is regarded as the image of a
%   uniform logical grid. Degree-seven quadrature integrates f(x(xi))*J on
%   that logical grid, with J reconstructed by the paired nine-point
%   derivative. The final normalization integrates constants exactly.

validateattributes(axis,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'axis');
if numel(axis) < 9
    error('ipm:SixthOrderQuadratureGrid', ...
        'Sixth-order mapped quadrature requires at least nine nodes.');
end
wasRow = isrow(axis);
axis = axis(:);
reference = linspace(0,1,numel(axis))';
metric = ipm.mesh.fdMatrix(reference,1,9)*axis;
if any(metric <= 0)
    error('ipm:SixthOrderQuadratureMetric', ...
        'The reconstructed logical-grid metric must remain positive.');
end
weights = ipm.mesh.quadrature(reference,8).*metric;
if any(weights <= 0)
    error('ipm:SixthOrderQuadratureWeights', ...
        'The mapped degree-seven quadrature norm must remain positive.');
end
weights = weights*(axis(end)-axis(1))/sum(weights);
if wasRow
    weights = weights';
end
end
