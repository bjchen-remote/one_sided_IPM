function result=compactBackground(model,scale,X,Y,second,sourceOnly)
%COMPACTBACKGROUND Keep the original stable Gaussian arithmetic inside.
base=model;base.kind='gaussian_primitive';
ell=scale.logC_l;omega=scale.logC_omega;inverse=exp(-ell);
x=(X-scale.X_shift)*inverse;y=Y*inverse;
outside=abs(x)>model.cutoff(1)|abs(y)>model.cutoff(1);
if sourceOnly
    if ~any(outside,'all')
        result=ipm_cr.model.source(base,scale,X,Y);return
    end
    result=zeros(size(X));
    inside=~outside;
    result(inside)=ipm_cr.model.source(base,scale,X(inside),Y(inside));
    taper=outside&abs(x)<model.cutoff(2)&abs(y)<model.cutoff(2);
    if any(taper,'all')
        value=ipm_cr.model.compactGaussian(model,x(taper),y(taper),false);
        result(taper)=exp(omega-ell)*value.x;
    end
    return
end
result=ipm_cr.model.gaussian(base,scale,X,Y,second);
origin=-scale.X_shift*inverse;
if any(outside,'all')
    value=ipm_cr.model.compactGaussian(model,x(outside),y(outside),second);
    result.H(outside)=result.H(outside)+exp(omega)*value.correction;
    result.x(outside)=exp(omega-ell)*value.x;
    result.y(outside)=exp(omega-ell)*value.y;
    if second
        result.xx(outside)=exp(omega-2*ell)*value.xx;
        result.xy(outside)=exp(omega-2*ell)*value.xy;
    end
end
if abs(origin)>model.cutoff(1)
    v=ipm_cr.model.compactGaussian(model,origin,0,false);
    result.H=result.H-exp(omega)*v.correction;
    result.originValue=exp(omega)*v.rho;
    result.originX=exp(omega-ell)*v.x;
end
result.H(X==0&Y==0)=0;
end
