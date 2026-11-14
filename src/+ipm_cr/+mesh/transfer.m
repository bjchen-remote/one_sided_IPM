function [out,audit]=transfer(s,old,new,allowContinuation)
% Transfer the active field without changing the state clock or references.
if nargin<4,allowContinuation=true;end
width=6;
out=s;key='W';if isfield(s,'eta'),key='eta';end
[out.(key),audit]=transfer_field(s.(key),old,new,allowContinuation,width);
if strcmp(key,'eta'),out.W=struct('core',zeros(size(new.core.X)),'polar',zeros(size(new.polar.X)));end
out.ctxId='';out.ctxKey=[];
out.rhs=struct('c_l',s.rhs.c_l,'c_r',s.rhs.c_r);
audit.integratedField=key;audit.timeAdvanced=0;audit.originalReferencesPreserved=true;
end

function [out,audit]=transfer_field(input,old,new,allowContinuation,width)
out=struct();counts=struct();maxNorm=0;
staircase=isfield(old.core,'stair');
for patch={'core','polar'}
    key=patch{1};target=new.(key);indices=target.activeIndices;
    x=target.X(indices);y=target.Y(indices);
    if strcmp(key,'polar')
        % Native axes avoid hypot roundoff outside the outer circle.
        [ir,it]=ind2sub(size(target.X),indices);
        r=target.r(ir);theta=target.theta(it);theta=theta(:);
    else
        r=hypot(x,y);theta=atan2(y,x);
    end
    coreCovered=covered(x,old.core.activeXBounds)&covered(y,old.core.activeYBounds);
    if staircase,coreCovered=ipm_cr.mesh.staircase('contains',old.core.stair,x,y);end
    polarCovered=covered(r,[old.polar.innerActiveRadius,old.polar.r(end)])& ...
        covered(theta,[old.polar.theta(1),old.polar.theta(end)]);
    fixedFrame=same_frame(old.(key),target,key);
    continueLocal=allowContinuation&&fixedFrame&&(~staircase||strcmp(key,'polar'));
    if strcmp(key,'core'),coveredBySame=coreCovered;else,coveredBySame=polarCovered;end
    same=coveredBySame;
    if continueLocal,same(:)=true;end
    if strcmp(key,'core'),coreOwner=same;else,coreOwner=~same;end
    % Standalone transfers retain fixed-frame continuation when requested.
    % Mapped remeshing uses the same patch only within its active coverage.
    values=zeros(size(x));
    for owner={'core','polar'}
        name=owner{1};take=coreOwner;if strcmp(name,'polar'),take=~take;end
        if ~any(take),continue;end
        p=old.(name);rows=p.activeRows;columns=p.activeColumns;
        if strcmp(name,'core')&&staircase
            M=ipm_cr.mesh.staircase('donors',p,x(take),y(take),width);
            donor=input.core(p.activeIndices);
        elseif strcmp(name,'core')
            M=ipm_cr.mesh.interpolate.tensor(p.y(rows),p.x(columns),y(take),x(take),width,continueLocal);
            donor=input.core(rows,columns);
        else
            M=ipm_cr.mesh.interpolate.tensor(p.r(rows),p.theta(columns),r(take),theta(take),width,continueLocal);
            donor=input.polar(rows,columns);
        end
        values(take)=difference_apply(M,donor(:));maxNorm=max(maxNorm,norm(M,inf));
    end
    out.(key)=zeros(size(target.X));out.(key)(indices)=values;
    counts.(key)=struct('samePatch',nnz(same),'crossPatch',nnz(~same), ...
        'samePatchBoundaryExtension',nnz(same&~coveredBySame),'fixedPatchFrame',fixedFrame);
end
out=ipm_cr.mesh.exchange(out,new);
audit=struct('kind','independent_active_patch_tensor_same_time_transfer_v2', ...
    'routes',counts,'mapInfinityNorm',maxNorm,'ghostValuesRead',false, ...
    'activeOverlapProjected',false,'polynomialDegree',width-1,'interpolationWidth',width, ...
    'fixedFrameContinuationAllowed',allowContinuation, ...
    'continuumAccuracyCertified',false);
end

function yes=same_frame(a,b,patch)
% Complete axes keep the patch frame; active endpoints may move with their
% third interior nodes. Counts and full endpoints identify redistribution,
% as distinct from cross-method import or patch-domain replacement.
if strcmp(patch,'core'),first='y';second='x';else,first='r';second='theta';end
yes=isequal(size(a.X),size(b.X))&& ...
    isequal(a.(first)([1,end]),b.(first)([1,end]))&& ...
    isequal(a.(second)([1,end]),b.(second)([1,end]));
end

function yes=covered(value,bounds)
% Use each endpoint's local ulp, as in localMap's roundoff repair.
yes=value>=bounds(1)-64*eps(max(1,abs(bounds(1))))& ...
    value<=bounds(2)+64*eps(max(1,abs(bounds(2))));
end

function v=difference_apply(M,donor)
[row,column,weight]=find(M);[~,anchor]=max(abs(M),[],2);anchor=full(anchor);
base=donor(anchor);v=base+accumarray(row,weight.*(donor(column)-base(row)),[size(M,1),1]);
end
