function prefix = continuousTrustedPrefix(mask)
%IPM.OUTPUT.CONTINUOUSTRUSTEDPREFIX Keep trust only before the first failure.

if ~islogical(mask) || ~isvector(mask)
    error('ipm:TrustedPrefixMask', ...
        'The trusted mask must be a logical vector.');
end
prefix = mask;
firstFailure = find(~mask,1,'first');
if ~isempty(firstFailure)
    prefix(firstFailure:end) = false;
end
end
