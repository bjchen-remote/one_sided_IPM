function q=backgroundSource(model,scale,X,Y)
%BACKGROUNDSOURCE Analytic G_x without forming unused density/derivatives.
% This is the same Gaussian source as background().x. Extreme scales use
% that routine's combined-log arithmetic instead of ordinary products.
assert(isequal(size(X),size(Y))&&isreal(X)&&isreal(Y)&& ...
    all(isfinite(X),'all')&&all(isfinite(Y),'all'), ...
    'ipm_cr:BackgroundQueries','Use matching finite real queries.');
if isfield(model,'tailCorrection')&&model.tailCorrection
    gaussianModel=model;gaussianModel.tailCorrection=false;
    q=ipm_cr.explicit.backgroundSource(gaussianModel,scale,X,Y);
    assert(isfield(scale,'physicalTime'),'ipm_cr:TailBackgroundClock');
    tail=ipm_cr.explicit.dipoleTail(model,scale,X,Y,false,true);
    low=0;if isfield(scale,'physicalTimeLow'),low=scale.physicalTimeLow;end
    q=q+scale.physicalTime*tail.x+low*tail.x;
    assert(all(isfinite(q),'all'),'ipm_cr:TailBackgroundRange');return
end
if strcmp(model.kind,'none'),q=zeros(size(X));return;end
ell=scale.logC_l;omega=scale.logC_omega;
if abs(ell)<=80&&abs(omega)<=80
    z=((X-scale.X_shift)*exp(-ell)-model.centerX)/model.widthX;
    y=Y*(exp(-ell)/model.widthY);
    coefficient=model.amplitude*exp(omega-ell);
    if isfinite(coefficient)&&all(isfinite(z),'all')&&all(isfinite(y),'all')
        w=model.skewWeight;
        q=coefficient*(((1-w)*exp(-z.^2)+(w/1.25)*exp(-((z-.8)/1.25).^2)).*exp(-y.^2));
        return;
    end
end
b=ipm_cr.explicit.background(model,scale,X,Y);q=b.x;
end
