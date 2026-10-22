function derivatives = derivatives(values,grid,options)
%POLAR_DERIVATIVES Cartesian and polar derivatives on a half-disk tensor.
% d.rx=d/dX, d.ry=d/dY, d.rr=d/dr (FIRST radial derivative),
% d.theta=d/dtheta.  d.radialSecond and d.angularSecond are second
% derivatives.  options.computeLaplacian=true additionally returns d.laplace.
% All origin samples must be the identical single scalar.  Cartesian first
% and second derivatives there use its Cartesian polynomial originFit,
% reading the unchanged first five coefficients, not a zero radial slope.
if nargin<3,options=struct();end
validateattributes(options,{'struct'},{'scalar'},mfilename,'options');
assert(isstruct(grid)&&isscalar(grid)&&isfield(grid,'originFit')&& ...
    isfield(grid,'kind')&&strcmp(grid.kind,'independent_polar_half_disk_tensor_v1'), ...
    'polar:Geometry','Use a geometry produced by ipm_c.mesh.build.');
assert(isnumeric(values)&&isreal(values)&& ...
    isequal(size(values),[grid.nr,grid.ntheta])&&all(isfinite(values),'all'), ...
    'polar:Field','The field must be a finite real nr-by-ntheta array.');
originValue=values(1,1);
assert(all(values(1,:)==originValue),'polar:OriginValueMismatch', ...
    'The complete origin row represents one scalar and must agree exactly.');
computeLaplacian=false;
if isfield(options,'computeLaplacian'),computeLaplacian=options.computeLaplacian;end
assert(isscalar(computeLaplacian)&& ...
    (islogical(computeLaplacian)||(isnumeric(computeLaplacian)&& ...
    isreal(computeLaplacian)&&isfinite(computeLaplacian)&& ...
    any(computeLaplacian==[0,1]))), ...
    'polar:LaplacianOption','computeLaplacian must be a scalar logical or zero/one.');
assert(~isfield(options,'poleRegularization')||isequal(options.poleRegularization,false), ...
    'polar:RetiredPoleRegularization','The experimental pole partition has been retired.');
% Apply the same local-polynomial stencil to differences about its own
% output node. Every local constant plateau is annihilated before the sum.
resFr=difference_apply(grid.Dr,values);resFrr=difference_apply(grid.Drr,values);
resFt=difference_apply(grid.Dtheta,values.').';
resFtt=difference_apply(grid.Dthetatheta,values.').';
fr=resFr;frr=resFrr;ft=resFt;ftt=resFtt;
positive=2:grid.nr;inverseR=1./grid.r(positive);
rx=zeros(size(values));ry=rx;
rx(positive,:)=resFr(positive,:).*grid.cosTheta- ...
    (resFt(positive,:).*inverseR).*grid.sinTheta;
ry(positive,:)=resFr(positive,:).*grid.sinTheta+ ...
    (resFt(positive,:).*inverseR).*grid.cosTheta;
fitValues=values(grid.originFit.rows,:)-originValue;
scaledCoefficients=grid.originFit.map*fitValues(:);
fitScale=grid.originFit.scale;
gx=scaledCoefficients(1)/fitScale;gy=scaledCoefficients(2)/fitScale;
hxx=2*scaledCoefficients(3)/fitScale^2;
hxy=scaledCoefficients(4)/fitScale^2;
hyy=2*scaledCoefficients(5)/fitScale^2;
rx(1,:)=gx;ry(1,:)=gy;
fr(1,:)=gx*grid.cosTheta+gy*grid.sinTheta;
frr(1,:)=hxx*grid.cosTheta.^2+ ...
    2*hxy*(grid.cosTheta.*grid.sinTheta)+hyy*grid.sinTheta.^2;
ft(1,:)=0;ftt(1,:)=0;
derivatives=struct('rx',rx,'ry',ry,'rr',fr,'theta',ft, ...
    'radialSecond',frr,'angularSecond',ftt, ...
    'origin',struct('value',originValue,'rx',gx,'ry',gy,'xx',hxx, ...
    'xy',hxy,'yy',hyy,'laplace',hxx+hyy, ...
    'closure',sprintf('shared-value Cartesian degree-%d fit; first-five-coefficient derivative jet',grid.originFit.degree)), ...
    'poleRegularization',false);
if computeLaplacian
    laplace=zeros(size(values));
    laplace(positive,:)=resFrr(positive,:)+resFr(positive,:).*inverseR+ ...
        resFtt(positive,:).*inverseR.^2;
    laplace(1,:)=hxx+hyy;
    derivatives.laplace=laplace;
end
fields={'rx','ry','rr','theta','radialSecond','angularSecond'};
if computeLaplacian,fields{end+1}='laplace';end
for k=1:numel(fields)
    assert(all(isfinite(derivatives.(fields{k})),'all'), ...
        'polar:NonfiniteDerivative','A polar derivative is not representable.');
end
end

function result=difference_apply(D,values)
[rows,columns,weights]=find(D);
differences=values(columns,:)-values(rows,:);
accumulator=sparse(rows,(1:numel(rows)).',weights,size(D,1),numel(rows));
result=accumulator*differences;
end
