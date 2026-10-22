function [result,rate,position,partials]=transport(rho,varargin)
%TRANSPORT Compatibility adapter; new callers use op.prepare/op.apply.
% Only this adapter preserves the historical output-count convention.
if numel(varargin)==1||numel(varargin)==2
    g=varargin{1};options=struct();
    if numel(varargin)==2,options=varargin{2};end
    prepared=ipm_cr.op.prepare(rho,g,options);
    result=@legacy_apply;rate=[];position=[];partials=[];
elseif numel(varargin)==3||numel(varargin)==4
    vx=varargin{1};vy=varargin{2};g=varargin{3};options=struct();
    if numel(varargin)==4,options=varargin{4};end
    prepared=ipm_cr.op.prepare(rho,g,options);
    request=struct('cfl',nargout>1,'partials',nargout>3);
    value=ipm_cr.op.apply(prepared,vx,vy,request);
    result=value.B;rate=value.rate;position=value.position;partials=value.partials;
else
    error('ipm_cr:ExplicitTransportArguments', ...
        'Use prepare(rho,g,options) and apply(prepared,vx,vy,request).');
end
    function [B,rate,position,partials]=legacy_apply(vx,vy,computeCfl)
        if nargin<3,computeCfl=true;end
        validateattributes(computeCfl,{'logical'},{'scalar'});
        request=struct('cfl',computeCfl&&nargout>1,'partials',nargout>3);
        value=ipm_cr.op.apply(prepared,vx,vy,request);
        B=value.B;rate=value.rate;position=value.position;partials=value.partials;
    end
end
