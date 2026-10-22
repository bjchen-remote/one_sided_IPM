function s=state(s,form)
%STATE Lossless field-name adapter at the public/legacy numerical boundary.
% No field arithmetic, interpolation, reordering or cache serialization.
if strcmp(form,'new'),from={'rho','residual','flow'};to={'W','eta','rhs'};
elseif strcmp(form,'old'),from={'W','eta','rhs'};to={'rho','residual','flow'};
else,error('ipm_cr:StateForm','Use new or old state field names.');end
for k=1:numel(from)
    if isfield(s,from{k})
        assert(~isfield(s,to{k}),'ipm_cr:StateAlias','State contains duplicate representations.');
        s.(to{k})=s.(from{k});s=rmfield(s,from{k});
    end
end
end
