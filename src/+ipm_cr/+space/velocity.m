classdef velocity
    %FIELD Evaluation maps and velocity recovery for the same FEM field.
    methods (Static)
        function plan=native(a,b,g,chart)
            % Evaluate with sparse 1D maps on the original sites.
            [Ea,Da]=ipm_cr.space.basis.evaluate(a,a.nodes);
            [Eb,Db]=ipm_cr.space.basis.evaluate(b,b.nodes);
            plan=struct('native',true,'chart',chart,'Da',Da,'Db',Db,'shape',size(g.X));
            plan.Ea=Ea;plan.Eb=Eb;
            if strcmp(chart,'polar')
                plan.radius=g.r(:);plan.cosTheta=g.cosTheta;plan.sinTheta=g.sinTheta;
            end
        end

        function plan=plan(a,b,x,y,chart)
            % Build sparse tensor evaluation maps for arbitrary physical queries.
            assert(isequal(size(x),size(y)),'ipm_cr:FEMQueryShape');shape=size(x);x=x(:);y=y(:);
            if strcmp(chart,'core'),qa=x;qb=y;else,qa=log(hypot(x,y));qb=atan2(y,x);end
            qa=roundoff_clamp(qa,a.breaks);qb=roundoff_clamp(qb,b.breaks);
            [Ea,Da]=ipm_cr.space.basis.evaluate(a,qa);[Eb,Db]=ipm_cr.space.basis.evaluate(b,qb);
            plan=struct('E',row_tensor(Ea,Eb,chart),'Da',row_tensor(Da,Eb,chart), ...
                'Db',row_tensor(Ea,Db,chart),'shape',shape,'x',x,'y',y,'chart',chart);
        end

        function v=evaluate(plan,U)
            % Recover the stream function and physical velocity from FEM coefficients.
            if isfield(plan,'native')&&plan.native
                    if strcmp(plan.chart,'core')
                        psi=plan.Eb*U*plan.Ea.';u1=-plan.Db*U*plan.Ea.';u2=plan.Eb*U*plan.Da.';
                    else
                        psi=plan.Ea*U*plan.Eb.';dt=plan.Da*U*plan.Eb.';da=plan.Ea*U*plan.Db.';
                        u1=-(plan.sinTheta.*dt+plan.cosTheta.*da)./plan.radius;
                        u2=(plan.cosTheta.*dt-plan.sinTheta.*da)./plan.radius;
                    end
                    v=struct('psi',psi,'u1',u1,'u2',u2);return
                end
            psi=plan.E*U(:);a=plan.Da*U(:);b=plan.Db*U(:);
            if strcmp(plan.chart,'core'),u1=-b;u2=a;
            else
                r=hypot(plan.x,plan.y);co=plan.x./r;si=plan.y./r;
                u1=-(si.*a+co.*b)./r;u2=(co.*a-si.*b)./r;
            end
            v=struct('psi',reshape(psi,plan.shape),'u1',reshape(u1,plan.shape),'u2',reshape(u2,plan.shape));
        end
    end
end

function q=roundoff_clamp(q,interval)
lo=interval(1);hi=interval(end);tol=64*eps(max([1,abs(lo),abs(hi)]));
assert(all(q>=lo-tol&q<=hi+tol),'ipm_cr:FEMQueryOutside','Query is outside its FEM patch.');
q=max(lo,min(hi,q));
% Polar physical points make an r -> (x,y) -> hypot -> log round trip.
% Snap only a few coordinate ulps to an existing element face so its
% declared averaged derivative does not become an accidental one-sided one.
[distance,index]=min(abs(q-interval(:).'),[],2);
onFace=distance<=8*eps(max(1,abs(q)));q(onFace)=interval(index(onFace));
end
function M=row_tensor(A,B,chart)
if size(A,1)==0,M=sparse(0,size(A,2)*size(B,2));return;end
[ca,va]=row_lists(A);[cb,vb]=row_lists(B);n=size(A,1);na=size(A,2);nb=size(B,2);
ii={};jj={};vv={};count=0;
for a=1:size(ca,2),for b=1:size(cb,2)
    use=ca(:,a)>0&cb(:,b)>0&va(:,a)~=0&vb(:,b)~=0;if ~any(use),continue;end
    count=count+1;ii{count}=find(use); %#ok<AGROW>
    if strcmp(chart,'core'),jj{count}=cb(use,b)+(ca(use,a)-1)*nb;
    else,jj{count}=ca(use,a)+(cb(use,b)-1)*na;end %#ok<AGROW>
    vv{count}=va(use,a).*vb(use,b); %#ok<AGROW>
end,end
if count==0,M=sparse(n,na*nb);else,M=sparse(vertcat(ii{:}),vertcat(jj{:}),vertcat(vv{:}),n,na*nb);end
end
function [columns,values]=row_lists(A)
[c,r,v]=find(A');counts=accumarray(r,1,[size(A,1),1]);w=max([0;counts]);
columns=zeros(size(A,1),w);values=columns;counts(:)=0;
for j=1:numel(r),i=r(j);counts(i)=counts(i)+1;columns(i,counts(i))=c(j);values(i,counts(i))=v(j);end
end
