function [options,boundary]=options(input)
%OPTIONS Explicit finite-source-disk exterior contract for time evolution.
% The direct FEM solver also accepts independently supplied Green boundary
% data. The time-dependent adapter presently supports zero exterior source
% only; this declares a model, not a bound on an unknown physical tail.
% These options select mathematics/discretization, not run-time accuracy
% gates. Independent tests qualify errors before experimental evolution.
if nargin<1,input=struct();end
options=struct('degree',5,'quadOrder',8,'stride',[],'nodalLayout','native', ...
    'outerBoundary','forced_dtn','exteriorModel','finite_source_disk_zero_exterior', ...
    'cornerCorrection',true,'linearSolver','sparse_direct');
assert(isstruct(input)&&isscalar(input)&&all(ismember(fieldnames(input),fieldnames(options))), ...
    'ipm_cr:FEMOptions','Unknown FEM runtime option.');
for name=fieldnames(input)',options.(name{1})=input.(name{1});end
validateattributes(options.degree,{'numeric'},{'scalar','real','finite','integer','>=',2,'<=',6});
validateattributes(options.quadOrder,{'numeric'},{'scalar','real','finite','integer','>=',options.degree+1});
assert(text_equal(options.nodalLayout,'native'),'ipm_cr:FEMNativeGrid', ...
    'Production FEM unknowns must use native transport nodes; gll_reference is a direct-solver comparison only.');
options.nodalLayout='native';
if ~isempty(options.stride)
    validateattributes(options.stride,{'numeric'},{'scalar','real','finite','integer','positive'});
    assert(options.stride==options.degree,'ipm_cr:FEMNativeGrid', ...
        'Native FEM stride must be empty or equal the target element degree; it cannot create another grid.');
end
validateattributes(options.cornerCorrection,{'logical'},{'scalar'});
assert(text_equal(options.linearSolver,'sparse_direct')||text_equal(options.linearSolver,'modal'), ...
    'ipm_cr:FEMOptions','Use linearSolver sparse_direct or modal.');
options.linearSolver=char(options.linearSolver);
assert(text_equal(options.outerBoundary,'forced_dtn')&& ...
    text_equal(options.exteriorModel,'finite_source_disk_zero_exterior'), ...
    'ipm_cr:FEMExterior','The FEM runtime requires the explicit finite-source-disk zero-exterior model.');
options.outerBoundary=char(options.outerBoundary);options.exteriorModel=char(options.exteriorModel);
if nargout>1
    boundary=struct('kind','forced_dtn','exteriorStatus','declared_compact_source', ...
        'value',@(X,Y)zeros(size(X)),'closureModel',options.exteriorModel);
end
end

function yes=text_equal(value,target)
yes=((ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)))&&strcmp(value,target);
end
