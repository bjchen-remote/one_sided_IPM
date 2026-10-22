function report=test_ipm_c_transport_manufactured(options)
%TEST_IPM_C_TRANSPORT_MANUFACTURED Actual production WENO path, zero Poisson.
% The optional pole-partition experiment is retired. These records retain
% the complete-domain source, transport and Q errors of the direct method.
if nargin<1,options=struct();end
counts=[65,129,257];
if isfield(options,'counts'),counts=options.counts;end
assert(numel(counts)==3&&all(counts>=17)&&all(diff(counts)>0)&& ...
    counts(2)==2*(counts(1)-1)+1&&counts(3)==2*(counts(2)-1)+1, ...
    'ipm_c:TransportRefinementCounts','Use three nested factor-two tensor counts.');
R=4;alpha=3;records=cell(1,3);
for level=1:3
    n=counts(level);h=R*sinh(alpha/(n-1))/sinh(alpha);
    g=ipm_c.mesh.build(struct('nr',n,'ntheta',n,'rmax',R,'firstSpacing',h));
    [rho,rx,ry,rxx,rxy]=gaussian_primitive(g.X,g.Y);
    d=ipm_c.field.derivatives(rho,g);
    first=false(size(rho));first(2,:)=true;transition=g.R>=.25 & g.R<=4;
    origin=g.R==0;cases=cell(1,2);
    for k=1:2
        speeds=[1,1000];dilations=[0,.2];c=dilations(k);
        vx=speeds(k)+c*g.X;vy=c*g.Y;
        exact=-vx.*rx-vy.*ry;exactQ=-c*rx-vx.*rxx-vy.*rxy;
        B=ipm_c.field.transport(rho,vx,vy,g);
        q=ipm_c.field.derivatives(B,g);
        assert(all(isfinite(B),'all')&&all(isfinite(q.rx),'all'), ...
            'ipm_c:FiniteTransport');
        constant=ipm_c.field.transport(3.25*ones(size(rho)),vx,vy,g);
        assert(all(constant==0,'all'),'polar:ConstantTransport');
        assert(all(B(1,:)==B(1,1)),'polar:SingleOriginTransport');
        cases{k}=struct('speed',speeds(k),'dilation',c, ...
            'transport',errors(B,exact,first,transition,origin), ...
            'Q',errors(q.rx,exactQ,first,transition,origin));
    end
    records{level}=struct('nr',n,'ntheta',n,'firstSpacing',h, ...
        'source',errors(d.rx,rx,first,transition,origin),'cases',{cases});
end
originWall=origin_wall_continuity();
report=struct('kind','production_polar_direct_transport_manufacture_v4', ...
    'passed',true,'validationScope','finite actual path, exact constants, shared origin and origin-wall continuity', ...
    'records',{records},'retiredPolePartitionUsed',false,'nativePdeCalls',0,'poissonCalls',0, ...
    'physicalTimeAdvanced',0,'nativeTrajectoryQualificationClaimed',false, ...
    'originRayWallContinuity',originWall);
end

function report=origin_wall_continuity()
% A smooth signal immediately upstream of the origin must reach the first
% positive ring continuously as theta approaches the same physical wall.
% The C-infinity bump is supported only at X<0, so independent positive-ray
% extrapolation cannot see it, while the actual signed wall stencil can.
% This tests the public transport's physical continuity, without inspecting
% or reproducing its ghost values. Fixed radial spacing isolates the seam
% from ordinary radial-refinement truncation error.
records=cell(1,2);h=.03;nr=17;counts=[257,513];
for fitIndex=1:2
    degrees=[4,6];fitDegree=degrees(fitIndex);gaps=zeros(1,2);responses=zeros(1,2);
    for level=1:2
        g=ipm_c.mesh.build(struct('nr',nr,'ntheta',counts(level), ...
            'rmax',(nr-1)*h,'firstSpacing',h,'originFitDegree',fitDegree));
        coordinate=(g.X+2*h)/(1.75*h);bump=zeros(size(g.X));
        inside=abs(coordinate)<1;
        bump(inside)=exp(1-1./(1-coordinate(inside).^2));
        base=g.X;signal=base+.2*bump;
        vx=ones(size(base));vy=zeros(size(base));
        response=ipm_c.field.transport(signal,vx,vy,g)- ...
            ipm_c.field.transport(base,vx,vy,g);
        responses(level)=response(2,1);
        gaps(level)=abs(response(2,2)-response(2,1));
    end
    assert(min(abs(responses))>1e-3,'polar:OriginWallSignal', ...
        'The actual signed-wall transport must resolve the upstream signal.');
    assert(gaps(2)<.1*abs(responses(2))&&gaps(2)<.7*gaps(1), ...
        'polar:OriginRayWallSeam', ...
        ['The first interior ray must approach the signed-wall response ', ...
         'under angular refinement at fixed radial spacing.']);
    records{fitIndex}=struct('originFitDegree',fitDegree,'radialSpacing',h, ...
        'radialCount',nr,'angularCounts',counts,'wallResponses',responses, ...
        'firstInteriorRayWallGaps',gaps,'relativeFineGap',gaps(2)/abs(responses(2)), ...
        'gapReductionFactor',gaps(1)/gaps(2),'passed',true);
end
report=struct('definition','first_positive_ring_transport_continuity_at_signed_wall', ...
    'manufacturedDensity','X plus smooth compact bump supported entirely at negative X', ...
    'velocity',[1,0],'records',{records});
end

function e=errors(value,reference,first,transition,origin)
v=abs(value-reference);
e=struct('fullInfinity',max(v,[],'all'),'firstRingInfinity',max(v(first)), ...
    'transitionInfinity',max(v(transition)),'originInfinity',max(v(origin)));
end

function [rho,rx,ry,rxx,rxy]=gaussian_primitive(X,Y)
A=1;a=-.25;sx=.7;sy=.6;w=.25;z=(X-a)/sx;v=(z-.8)/1.25;
envelope=exp(-(Y/sy).^2);
rho=A*sx*sqrt(pi)/2*((1-w)*erfc(-z)+w*erfc(-v)).*envelope;
rx=A*((1-w)*exp(-z.^2)+(w/1.25)*exp(-v.^2)).*envelope;
ry=-(2*Y/sy^2).*rho;
rxx=-(2*A/sx)*((1-w)*z.*exp(-z.^2)+(w/1.25^2)*v.*exp(-v.^2)).*envelope;
rxy=-(2*Y/sy^2).*rx;
end
