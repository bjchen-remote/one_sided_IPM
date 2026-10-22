function reconstructed = weno5Reconstruct(values,geometry,epsilon)
%IPM.FIELD.WENO5RECONSTRUCT Fifth-order WENO reconstruction from point values.
%   VALUES contains independent nodal signals by row. GEOMETRY describes
%   one side of every interior face on an arbitrary monotone axis.

numberOfLines = size(values,1);
numberOfFaces = size(geometry.indices,1);
samples = zeros(numberOfLines,numberOfFaces,5,'like',values);
for node = 1:5
    samples(:,:,node) = values(:,geometry.indices(:,node));
end
localScale = max(abs(samples),[],3);
safeScale = max(localScale,sqrt(realmin(class(values))));
normalized = samples./safeScale;

fullPolynomial = sum(normalized.*reshape( ...
    geometry.fullWeights,1,numberOfFaces,5),3);
candidatePolynomial = zeros(numberOfLines,numberOfFaces,3,'like',values);
for candidate = 1:3
    weights = reshape(geometry.candidateWeights(:,candidate,:), ...
        1,numberOfFaces,5);
    candidatePolynomial(:,:,candidate) = sum(normalized.*weights,3);
end

slopes = diff(normalized,1,3)./reshape( ...
    geometry.deltaCoordinates,1,numberOfFaces,4);
indicators = cat(3, ...
    slopes(:,:,1).^2+slopes(:,:,2).^2, ...
    slopes(:,:,2).^2+slopes(:,:,3).^2, ...
    slopes(:,:,3).^2+slopes(:,:,4).^2);
globalDifference = sum(normalized.*reshape( ...
    geometry.differenceWeights,1,numberOfFaces,5),3);
globalIndicator = globalDifference.^2;
indicatorPower = indicators.^2;
globalPower = globalIndicator.^2;
alpha = (1/3)*(1+reshape(globalPower,numberOfLines,numberOfFaces,1) ./ ...
    (indicatorPower+epsilon));
nonlinearWeights = alpha./sum(alpha,3);
lowOrder = sum(nonlinearWeights.*candidatePolynomial,3);
globalAverageWeight = 1./(1+globalPower.*sum( ...
    1./(indicatorPower+epsilon),3));
normalizedResult = globalAverageWeight.*fullPolynomial + ...
    (1-globalAverageWeight).*lowOrder;

reconstructed = normalizedResult.*safeScale;
end
