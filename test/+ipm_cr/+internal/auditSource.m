function f=auditSource(callback,x,y)
%AUDITSOURCE A missing audit and an invalid audit are distinct conditions.
f=callback(x,y);
if isscalar(f),f=repmat(f,size(x));end
assert(isnumeric(f)&&isreal(f)&&isequal(size(f),size(x))&&all(isfinite(f),'all'), ...
    'ipm_cr:SourceAudit','The off-grid source must be finite, real and match the query shape.');
end
