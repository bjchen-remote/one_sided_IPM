classdef axes
    %AXES Coordinate maps, allocation and remesh geometry.
    methods (Static)
        function [out,info,varargout]=staircase(action,varargin)
            %STAIRAXES Smooth point allocation on the existing physical coordinate maps.
            switch action
                case 'core',[out,info]=core(varargin{:});
                case 'interfaceCore',[out,info]=interface_core(varargin{:});
                case 'resize',out=resize_axes(varargin{:});info=[];
                case 'radial',[out,info]=radial(varargin{:});
                case 'matchRadial',[out,info,varargout{1}]=match_radial(varargin{:});
                otherwise,error('ipm_cr:StairAxes','Unknown staircase axis operation.');
            end
        end

        function grid = polar(options)
            %POLAR_GRID Independent stretched polar tensor geometry on a half disk.
            % Arrays have size [nr,ntheta].  The repeated r=0 row denotes ONE point.
            % options: nr=129, ntheta=129, rmax=1e10, firstSpacing=.02,
            % domain='semicircle', radialMap='sinh', stencilWidth=7.
            % radialMap='uniform_core' adds a nearly uniform core.
            % The default r(s)=R*sinh(alpha*s)/sinh(alpha) has a smooth ODD continuation
            % across the signed physical wall.  Explicit radialMap='exp_gaps' retains
            % the old increasing-gap geometry for historical diagnostic comparisons.
            % This constructor returns axes, coordinates and the radial map metric;
            % it does not build origin jets or transport derivative operators.
            % No rectangle is replaced by a half disk silently.
            % Requires the original ipm_cr.mesh.interpolate.metricDerivative on the MATLAB search path.
            if nargin<1,options=struct();end
            validateattributes(options,{'struct'},{'scalar'},mfilename,'options');
            domain=char(option(options,'domain','semicircle'));
            if strcmp(domain,'rectangle')
                error('polar:RectangleDomainNotImplemented', ...
                    'A rectangular polar mapping is not implemented; choose semicircle explicitly.');
            end
            assert(strcmp(domain,'semicircle'),'polar:Domain', ...
                'The only implemented polar domain is semicircle.');
            nr=option(options,'nr',129);ntheta=option(options,'ntheta',129);
            rmax=option(options,'rmax',1e10);firstSpacing=option(options,'firstSpacing',.02);
            radialMap=char(option(options,'radialMap','sinh'));
            assert(ismember(radialMap,{'sinh','exp_gaps','uniform_core'}),'polar:RadialMap', ...
                'radialMap must be sinh, uniform_core, or the explicit historical exp_gaps option.');
            width=option(options,'stencilWidth',7);
            assert(~any(isfield(options,{'originRings','originFitDegree'})), ...
                'ipm_cr:RetiredOption','originRings/originFitDegree no longer affect this origin-free solver.');
            validateattributes(nr,{'numeric'},{'scalar','real','finite','integer','>=',7});
            validateattributes(ntheta,{'numeric'},{'scalar','real','finite','integer','>=',7});
            validateattributes(rmax,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(firstSpacing,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(width,{'numeric'},{'scalar','real','finite','integer','>=',3});
            assert(width<=min(nr,ntheta),'polar:StencilWidth', ...
                'The derivative stencil must fit both tensor axes.');
            intervals=nr-1;
            assert(log(firstSpacing)+log(intervals)<=log(rmax)+32*eps(max(1,abs(log(rmax)))), ...
                'polar:RadialSpacing','firstSpacing must not exceed the average radial spacing.');
            if strcmp(radialMap,'uniform_core')
                [r,uniformCoreInfo]=ipm_cr.mesh.axes.radialAxis(options);
                firstSpacing=uniformCoreInfo.firstSpacing;
                mapParameter=uniformCoreInfo.alpha;
            else
                logTarget=log(rmax)-log(firstSpacing);
                if abs(logTarget-log(intervals))<=32*eps(max(1,abs(logTarget)))
                    mapParameter=0;
                else
                    lower=0;upper=1;
                    while radial_ratio_log(intervals,upper,radialMap)<logTarget,upper=2*upper;end
                    for iteration=1:100
                        middle=lower+(upper-lower)/2;
                        if radial_ratio_log(intervals,middle,radialMap)<logTarget
                            lower=middle;
                        else
                            upper=middle;
                        end
                    end
                    mapParameter=lower+(upper-lower)/2;
                end
                r=zeros(nr,1);
                for k=1:intervals
                    if strcmp(radialMap,'exp_gaps')
                        r(k+1)=exp(log(firstSpacing)+geometric_logsum(k,mapParameter));
                    elseif mapParameter==0
                        r(k+1)=rmax*k/intervals;
                    else
                        r(k+1)=exp(log(rmax)+log_sinh(mapParameter*k/intervals)-log_sinh(mapParameter));
                    end
                end
            end
            r(2)=firstSpacing;r(end)=rmax;
            assert(all(isfinite(r))&&all(diff(r)>0),'polar:RadialRepresentability', ...
                'The requested radial axis must be finite and strictly increasing.');
            theta=linspace(0,pi,ntheta);
            cosTheta=cos(theta);sinTheta=sin(theta);
            cosTheta([1,end])=[1,-1];sinTheta([1,end])=0;
            if mod(ntheta,2)==1
                cosTheta((ntheta+1)/2)=0;sinTheta((ntheta+1)/2)=1;
            end
            X=r*cosTheta;Y=r*sinTheta;
            logical=linspace(0,1,nr)';radialMetric=ipm_cr.mesh.interpolate.metricDerivative(logical,1,7)*r;
            assert(all(radialMetric>0),'polar:GridMetric', ...
                'The mapped radial metric must remain positive.');
            grid=struct('kind','independent_polar_half_disk_tensor_v1','domain',domain, ...
                'nr',nr,'ntheta',ntheta,'rmax',rmax,'firstSpacing',firstSpacing, ...
                'radialMap',radialMap,'radialMapParameter',mapParameter, ...
                'signedWallMapSmoothAtOrigin',ismember(radialMap,{'sinh','uniform_core'}), ...
                'r',r,'theta',theta, ...
                'X',X,'Y',Y,'cosTheta',cosTheta,'sinTheta',sinTheta, ...
                'radialMetric',radialMetric,'stencilWidth',width, ...
                'originIsSinglePhysicalPoint',true,'originalRectanglePreserved',false);
            if strcmp(radialMap,'exp_gaps'),grid.radialGapRatio=exp(mapParameter);end
            if strcmp(radialMap,'uniform_core'),grid.uniformCoreInfo=uniformCoreInfo;end
        end

        function [r,drMetric,info]=redistributeRadial(original,options)
            %REDISTRIBUTERADIALAXIS Fixed smooth refinement of a uniform_core radial axis.
            % u is the OLD unit-index coordinate, 0<=u<=baseNr-1. New integer j samples
            % q(u)=u+c*integral_0^u chi(v)dv, chi=betainc(t,8,8) on the transition.
            % This is a one-time geometry construction, not ALE or time adaptation.
            % Returned drMetric differentiates r with respect to the NEW unit index.
            % Complete FD7 rows inside the exact-copy prefix retain old metric/N.
            if nargin<2,options=struct();end
            assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options), ...
                {'enabled','factor','preserveRadius','fullRefinementRadius'})), ...
                'ipm_cr:RadialRedistributionOptions','Unknown radial redistribution option.');
            enabled=redistributeRadial_option(options,'enabled',false);factor=redistributeRadial_option(options,'factor',2);
            preserve=redistributeRadial_option(options,'preserveRadius',4);full=redistributeRadial_option(options,'fullRefinementRadius',20);
            validateattributes(enabled,{'numeric','logical'},{'scalar','real','finite'});
            assert(any(enabled==[0,1]),'ipm_cr:RadialRedistributionOptions','enabled must be a logical flag.');
            assert(isstruct(original)&&isscalar(original)&&all(isfield(original, ...
                {'nr','r','radialMetric','radialMap'})), ...
                'ipm_cr:RadialRedistributionBase','Use a complete original ipm_cr.mesh.axes.polar geometry.');
            baseR=original.r(:);N=original.nr-1;
            assert(numel(baseR)==N+1&&N>=6&&baseR(1)==0&&all(isfinite(baseR))&& ...
                all(diff(baseR)>0)&&numel(original.radialMetric)==N+1, ...
                'ipm_cr:RadialRedistributionBase','The base radial axis must be finite and increasing from zero.');
            baseDr=original.radialMetric(:)/N;
            assert(all(isfinite(baseDr))&&all(baseDr>0),'ipm_cr:RadialRedistributionMetric');
            r=baseR;drMetric=baseDr;
            info=struct('enabled',logical(enabled),'baseNr',N+1,'actualNr',N+1, ...
                'monitorVersion','old_index_beta_8_8_integral_v1','geometryFixedInTime',true);
            if ~enabled,return;end
            validateattributes(factor,{'numeric'},{'scalar','real','finite','>=',1});
            validateattributes(preserve,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(full,{'numeric'},{'scalar','real','finite','>',preserve,'<=',baseR(end)});
            assert(strcmp(original.radialMap,'uniform_core')&&isfield(original,'uniformCoreInfo'), ...
                'ipm_cr:RadialRedistributionMap','Enabled redistribution requires the original uniform_core analytic map.');
            map=original.uniformCoreInfo;
            assert(all(isfield(map,{'a','alpha','power','tailCoefficient'}))&& ...
                all(isfinite([map.a,map.alpha,map.power,map.tailCoefficient]))&& ...
                map.a>0&&map.alpha>=0&&map.power==3&&map.tailCoefficient>0, ...
                'ipm_cr:RadialRedistributionMap','Invalid original uniform_core map parameters.');
            ua=inverse_radius(preserve);ub=inverse_radius(full);
            assert(0<ua&&ua<ub&&ub<=N,'ipm_cr:RadialRedistributionMap');
            integralN=N-(ua+ub)/2;
            K=round((factor-1)*integralN);Nnew=N+K;
            assert(isfinite(K)&&K>=0&&isfinite(Nnew)&&Nnew<flintmax, ...
                'ipm_cr:RadialRedistributionCount','The requested radial interval count is not representable.');
            c=K/integralN;protectedCount=floor(ua)+1;
            assert(protectedCount>=7,'ipm_cr:RadialRedistributionProtection', ...
                'The protected radial prefix must contain a complete seven-point stencil.');
            u=(0:Nnew)';preservedCount=protectedCount;
            if K==0
                % factor=1 or a rounded zero insertion is a literal geometry no-op.
                u=(0:N)';preservedCount=N+1;copiedMetricRows=(1:N+1)';
            else
                j=u;protected=j<=floor(ua);
                atUpper=ub+c*(ub-ua)/2;
                far=~protected&j>=atUpper;
                u(far)=(j(far)+c*(ua+ub)/2)/(1+c);
                transition=~protected&~far;
                lower=ua+zeros(nnz(transition),1);upper=ub+zeros(nnz(transition),1);
                target=j(transition);
                for iteration=1:80
                    middle=lower+(upper-lower)/2;
                    below=middle+c*monitor_integral(middle,ua,ub)<target;
                    lower(below)=middle(below);upper(~below)=middle(~below);
                end
                u(transition)=lower+(upper-lower)/2;
                u(1:protectedCount)=(0:protectedCount-1)';u(end)=N;
                r=analytic_radius(u);
                % Never reevaluate exact protected integer nodes through exp/log.
                r(1:protectedCount)=baseR(1:protectedCount);r(end)=baseR(end);
                assert(all(isfinite(r))&&all(diff(r)>0),'ipm_cr:RadialRedistributionAxis', ...
                    'The redistributed analytic radial axis must remain strictly increasing.');
                drMetric=ipm_cr.mesh.interpolate.metricDerivative((0:Nnew)',1,7)*r;
                rows=(1:Nnew+1)';stencilFirst=min(max(rows-3,1),Nnew+1-6);
                copiedMetricRows=find(stencilFirst+6<=protectedCount);
                drMetric(copiedMetricRows)=baseDr(copiedMetricRows);
            end
            assert(all(isfinite(drMetric))&&all(drMetric>0), ...
                'ipm_cr:RadialRedistributionMetric','The redistributed mapped radial metric must be finite and positive.');
            inverseResidual=max(abs(u+c*monitor_integral(u,ua,ub)-(0:Nnew)'));
            assert(inverseResidual<=256*eps(max(1,Nnew)), ...
                'ipm_cr:RadialRedistributionInverse','The smooth monitor inverse did not converge to machine brackets.');
            gaps=diff(r);
            info=struct('enabled',true,'kind','fixed_old_index_smooth_radial_refinement_v1', ...
                'monitorVersion','old_index_beta_8_8_integral_v1','geometryFixedInTime',true, ...
                'factor',factor,'preserveRadius',preserve,'fullRefinementRadius',full, ...
                'baseNr',N+1,'actualNr',Nnew+1,'baseIntervals',N,'actualIntervals',Nnew, ...
                'addedIntervals',K,'ua',ua,'ub',ub,'monitorAmplitude',c, ...
                'actualTailDensityFactor',1+c,'integratedMonitorAtOuter',integralN, ...
                'protectedPrefixCount',protectedCount,'preservedPrefixCount',preservedCount, ...
                'preservedLastRadius',r(preservedCount), ...
                'copiedMetricPrefixCount',numel(copiedMetricRows),'copiedMetricRows',copiedMetricRows, ...
                'metricDerivativeWidth',7,'oldIndexAtNewNodes',u, ...
                'baseR',baseR,'baseDrMetric',baseDr,'maximumInverseResidual',inverseResidual, ...
                'minimumRadialGap',min(gaps),'maximumRadialGap',max(gaps), ...
                'maximumAdjacentGapRatio',max([gaps(2:end)./gaps(1:end-1);gaps(1:end-1)./gaps(2:end)]), ...
                'analyticMap',struct('a',map.a,'alpha',map.alpha,'power',map.power, ...
                    'tailCoefficient',map.tailCoefficient), ...
                'metricDefinition','FD7 on new integer index; original radialMetric/baseIntervals on complete protected stencils');

                function value=analytic_radius(oldIndex)
                    s=oldIndex/N;
                    if map.alpha==0
                        tail=s.^map.power;
                    else
                        tail=zeros(size(s));positive=s>0;
                        tail(positive)=exp(map.power*(redistributeRadial_log_sinh(map.alpha*s(positive))-redistributeRadial_log_sinh(map.alpha)));
                    end
                    value=map.a*s+map.tailCoefficient*tail;
                    value(oldIndex==0)=0;value(oldIndex==N)=baseR(end);
                end

                function value=inverse_radius(radius)
                    if radius==baseR(end),value=N;return;end
                    lower=0;upper=N;
                    for iteration=1:80
                        middle=lower+(upper-lower)/2;sample=analytic_radius(middle);
                        if sample==radius,lower=middle;upper=middle;break;end
                        if sample<radius,lower=middle;else,upper=middle;end
                    end
                    value=lower+(upper-lower)/2;
                end
        end

        function [r,info]=radialAxis(options)
            %UNIFORMCOREAXIS Smooth odd linear core with an analytic stretched tail.
            % r(s)=a*s+(R-a)*(sinh(alpha*s)/sinh(alpha))^3, 0<=s<=1.
            % firstSpacing is a lower bound: increase it when the measured core would
            % occupy more than 60 percent of the logical axis. This retains enough tail
            % intervals to reach rmax smoothly instead of exhausting them as the core
            % grows. Set a so the first positive node equals the selected spacing.
            % Choose alpha so the tail at s=coreRadius/((nr-1)*spacing) is only the
            % declared fraction of coreRadius. Both positive terms increase and are odd.
            nr=radialAxis_option(options,'nr',129);R=radialAxis_option(options,'rmax',1e10);
            requestedSpacing=radialAxis_option(options,'firstSpacing',.02);
            core=radialAxis_option(options,'coreRadius',2);fraction=radialAxis_option(options,'coreTailFraction',.01);
            validateattributes(nr,{'numeric'},{'scalar','real','finite','integer','>=',7});
            validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(requestedSpacing,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(core,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(fraction,{'numeric'},{'scalar','real','finite','>',0,'<',1});
            N=nr-1;maximumCoreLogicalLocation=.6;
            h=max(requestedSpacing,core/(maximumCoreLogicalLocation*N));
            sFirst=1/N;sCore=core/(N*h);power=3;
            assert(N*h<R&&sCore<1&&core<R,'polar:UniformCoreRange', ...
                'The uniform core must fit inside the logical axis and physical radius.');
            targetTail=core*fraction;
            lower=0;upper=1;
            while tail_at_core(upper)>targetTail
                upper=2*upper;
                assert(isfinite(upper),'polar:UniformCoreParameter','The tail parameter must be finite.');
            end
            assert(tail_at_core(lower)>=targetTail,'polar:UniformCoreParameter', ...
                'The requested tail fraction lies outside the analytic map family.');
            for iteration=1:100
                middle=lower+(upper-lower)/2;
                if tail_at_core(middle)>targetTail,lower=middle;else,upper=middle;end
            end
            alpha=lower+(upper-lower)/2;
            a=linear_coefficient(alpha);tailCoefficient=R-a;
            assert(isfinite(a)&&a>0&&isfinite(tailCoefficient)&&tailCoefficient>0, ...
                'polar:UniformCoreParameter','The core and tail coefficients must both be positive.');
            s=(0:N)'/N;r=a*s+tailCoefficient*tail_ratio(s,alpha);
            r([1,2,end])=[0,h,R];
            assert(all(isfinite(r))&&all(diff(r)>0)&&abs(r(2)-h)<=64*eps(h), ...
                'polar:UniformCoreRepresentability','The analytic map must be increasing and retain its first spacing.');
            gaps=diff(r);inside=r(2:end)<=core;
            adjacentGapRatio=max([gaps(2:end)./gaps(1:end-1);gaps(1:end-1)./gaps(2:end)]);
            assert(any(inside),'polar:UniformCoreRange','The core must contain a positive radial node.');
            coreGaps=gaps(inside);observedTail=tail_at_core(alpha);
            info=struct('kind','positive_linear_plus_cubic_sinh_odd_axis_v2', ...
                'nr',nr,'rmax',R,'requestedFirstSpacing',requestedSpacing, ...
                'firstSpacing',h,'firstSpacingAdjusted',h>requestedSpacing, ...
                'a',a,'alpha',alpha,'power',power, ...
                'tailCoefficient',tailCoefficient,'coreRadius',core,'coreTailFraction',fraction, ...
                'coreLogicalLocation',sCore,'maximumCoreLogicalLocation',maximumCoreLogicalLocation, ...
                'positiveCoreNodeFraction',nnz(inside)/N, ...
                'tailIntervalCount',N-nnz(inside),'maximumAdjacentGapRatio',adjacentGapRatio, ...
                'tailAtCoreLogicalLocation',observedTail, ...
                'targetTailAtCoreLogicalLocation',targetTail, ...
                'coreTailRelativeResidual',abs(observedTail-targetTail)/targetTail, ...
                'positiveCoreNodeCount',nnz(inside),'coreGaps',coreGaps, ...
                'minimumCoreGap',min(coreGaps),'maximumCoreGap',max(coreGaps), ...
                'coreGapRatio',max(coreGaps)/min(coreGaps), ...
                'signedInverseMapOddAndSmooth',true,'originAndRadiusRetained',true, ...
                'mapDefinition','a*s+(rmax-a)*(sinh(alpha*s)/sinh(alpha))^3');

                function value=linear_coefficient(parameter)
                    eta=tail_ratio(sFirst,parameter);
                    value=(h-R*eta)/(sFirst-eta);
                end
                function value=tail_at_core(parameter)
                    value=(R-linear_coefficient(parameter))*tail_ratio(sCore,parameter);
                end
                function value=tail_ratio(location,parameter)
                    if parameter==0
                        value=location.^power;
                    else
                        value=zeros(size(location));positive=location>0;
                        value(positive)=exp(power*(radialAxis_log_sinh(parameter*location(positive))-radialAxis_log_sinh(parameter)));
                    end
                end
        end

        function axes=tensor(a,b,input)
            %TENSORAXIS Smooth sinh tensor axes with prescribed origin metric spacing.
            % x=a*sinh(alpha*xi)/sinh(alpha), xi in [-1,1];
            % y=b*sinh(beta*eta)/sinh(beta), eta in [0,1].
            % The requested spacing is J(0)*dlogical, not a secant approximation.
            % coreMap contains only constructor inputs and can be passed back to BUILD;
            % derived map parameters are stored separately. No PDE state is changed.
            if nargin<3,input=struct();end
            validateattributes(a,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(b,{'numeric'},{'scalar','real','finite','positive'});
            map=struct('enabled',true,'kind','sinh','intervals',[240,128], ...
                'minSpacing',[.01,.00375]);
            assert(isstruct(input)&&isscalar(input)&&all(ismember(fieldnames(input),fieldnames(map))), ...
                'ipm_cr:TensorMap','Unknown tensor map input.');
            names=fieldnames(input);
            for k=1:numel(names),map.(names{k})=input.(names{k});end
            validateattributes(map.enabled,{'logical'},{'scalar'});
            assert(map.enabled,'ipm_cr:TensorMap','Disabled mappings use the unchanged uniform BUILD path.');
            map.kind=char(validatestring(map.kind,{'sinh'}));
            validateattributes(map.intervals,{'numeric'},{'vector','numel',2,'real','finite','integer','>=',8});
            validateattributes(map.minSpacing,{'numeric'},{'vector','numel',2,'real','finite','positive'});
            map.intervals=reshape(map.intervals,1,2);map.minSpacing=reshape(map.minSpacing,1,2);
            assert(mod(map.intervals(1),2)==0,'ipm_cr:TensorMapOrigin', ...
                'An even number of x intervals is required for the exact origin.');
            [x,Jx,dx,dxi,xi,alpha]=axis_map(a,map.intervals(1),map.minSpacing(1),true);
            [y,Jy,dy,deta,eta,beta]=axis_map(b,map.intervals(2),map.minSpacing(2),false);
            parameters=struct('alpha',alpha,'beta',beta,'uniformAxes',[alpha==0,beta==0], ...
                'originMetricSpacing',[min(dx),min(dy)], ...
                'averagePhysicalSpacing',[2*a/map.intervals(1),b/map.intervals(2)], ...
                'logicalSpacing',[dxi,deta],'logicalBounds',[-1,1,0,1], ...
                'definition','analytic_J_times_logical_spacing');
            axes=struct('x',x,'y',y,'Jx',Jx,'Jy',Jy,'dxMetric',dx,'dyMetric',dy, ...
                'dxi',dxi,'deta',deta,'xi',xi,'eta',eta,'coreMap',map,'mapParameters',parameters);
        end

        function [out,info]=allocate(a,plan,halo)
            %STUDYALLOCATION Reallocate a coarse thin-staircase budget on native axes.
            % The first halo/ring coordinates and the box endpoints are fixed.
            validateattributes(plan,{'double'},{'vector','real','finite'});
            assert(ismember(numel(plan),[6,7,8]));
            nearOffset=0;
            if numel(plan)>=7
                nearOffset=plan(7);
                validateattributes(nearOffset,{'double'},{'integer','>=',-1,'<=',2});
            end
            plan=reshape(plan,1,[]);
            % Reuse one allocation in physical coordinates across the paired presets.
            % The 181-angle medium chart has twice the base coordinate intervals.
            refinement=1+double(numel(a.theta)==181);
            addedPairs=refinement*plan(1);removedRings=refinement*plan(2);center=plan(3);
            halfWidth=plan(4);warpStrength=plan(5);nearRings=refinement*plan(6);
            validateattributes([addedPairs,removedRings],{'double'},{'integer','>=',-8,'<=',30});
            validateattributes(nearRings,{'double'},{'integer','nonnegative','<=',30});
            assert(center>0&&halfWidth>0&&warpStrength>=0&&warpStrength<1);
            a.r=a.r(:).';a.drMetric=a.drMetric(:).';
            cap=min(a.x(end-halo),a.y(end-halo));
            capIndex=find(a.r<=cap,1,'last');
            assert(capIndex>halo,'ipm_cr:StudyGeometry','The inner ring must remain inside the donor box.');
            second=a.r(capIndex);
            baseCounts=[numel(a.x),numel(a.y),numel(a.r)-capIndex+halo+2,numel(a.theta)];
            tailAdjustment=0;
            if refinement==2
                % Count the coarse partner in the same physical frame. The three
                % ghost rows make the cropped sizes differ from the full-axis ratio.
                coarseR=a.r(end:-2:1);
                coarseCap=min(a.x(end-2*halo),a.y(end-2*halo));
                coarseNr=nnz(coarseR>coarseCap)+halo+2-plan(2);
                tailAdjustment=2*(coarseNr-1)+1-(baseCounts(3)-removedRings);
            end
            middle=(numel(a.x)+1)/2;x=a.x(middle:end);dx=a.dxMetric(middle:end);
            start=max(x(2),(center-halfWidth)*second);
            finish=min(x(end-halo-2),(center+halfWidth)*second);
            assert(finish>start,'ipm_cr:StudyGeometry','The target allocation interval is empty.');
            focus=min(finish,max(start,center*second));
            if warpStrength>0
                width=finish-start;t=(x-start)/width;active=t>0&t<1;
                z=t(active);
                bump=64*(z.*(1-z)).^3;
                bumpPrime=192*(z.*(1-z)).^2.*(1-2*z)/width;
                offset=x(active)-focus;
                x(active)=x(active)-warpStrength*offset.*bump;
                dx(active)=dx(active).*(1-warpStrength*(bump+offset.*bumpPrime));
            end
            if addedPairs>0
                z=interp1(x,0:numel(x)-1,[start,finish],'linear');
                [x,dx]=redistribute(x,dx,[z(1),z(2),addedPairs],addedPairs);
            elseif addedPairs<0
                % Keep the staircase and its active Cartesian columns inside the
                % unchanged prefix. Spend fewer columns on the unused tensor tail.
                tailStart=interp1(x,0:numel(x)-1,1.12*second,'linear');
                tailEnd=numel(x)-1;
                assert(tailEnd-tailStart>1.875*abs(addedPairs),'ipm_cr:StudyTailDensity', ...
                    'The Cartesian tail must have positive coordinate density.');
                [x,dx]=redistribute(x,dx,[tailStart,tailEnd,addedPairs],addedPairs);
            end
            y=a.y;dy=a.dyMetric;
            if numel(plan)==8
                % Refine the measured corner height with a C2 physical-coordinate map.
                strength=plan(8);
                validateattributes(strength,{'double'},{'scalar','nonnegative','<',1});
                lo=.06*second;hi=min(.44*second,y(end-halo-1));focusY=.24*second;
                t=(y-lo)/(hi-lo);inside=t>0&t<1;z=t(inside);
                bump=64*(z.*(1-z)).^3;
                bumpPrime=192*(z.*(1-z)).^2.*(1-2*z)/(hi-lo);
                offset=y(inside)-focusY;
                y(inside)=y(inside)-strength*offset.*bump;
                dy(inside)=dy(inside).*(1-strength*(bump+offset.*bumpPrime));
            end
            out=struct('x',[-fliplr(x(2:end)),x],'dxMetric',[fliplr(dx(2:end)),dx], ...
                'y',y,'dyMetric',dy,'r',a.r,'drMetric',a.drMetric,'theta',a.theta);
            if removedRings~=0||nearRings>0||tailAdjustment~=0
                n=numel(a.r)-1;
                windows=[capIndex+refinement*(2+nearOffset), ...
                    capIndex+refinement*(15+nearOffset),nearRings; ...
                    capIndex+refinement*20,n,tailAdjustment-removedRings-nearRings];
                [out.r,out.drMetric]=redistribute(a.r,a.drMetric,windows,tailAdjustment-removedRings);
            end
            counts=baseCounts+[2*addedPairs,0,tailAdjustment-removedRings,0];
            nominal=@(c)c(1)*c(2)+c(3)*c(4);
            assert(all(diff(out.x)>0)&&all(out.dxMetric>0)&& ...
                all(diff(out.r)>0)&&all(out.drMetric>0)&&all(diff(out.y)>0)&&all(out.dyMetric>0));
            newCap=min(out.x(end-halo),out.y(end-halo));
            assert(isequal(newCap,cap)&&isequal(out.r(1:capIndex+1),a.r(1:capIndex+1)));
            assert(isequal(out.x([1,end]),a.x([1,end]))&&out.r(end)==a.r(end));
            measure=struct('plan',plan,'targetBand',[start,finish],'secondActiveRadius',second, ...
                'referenceSizes',baseCounts,'sizes',counts,'referenceNominalCount',nominal(baseCounts), ...
                'nominalCount',nominal(counts),'allocationRefinement',refinement);
            info=struct('referenceX',a.x,'referenceDx',a.dxMetric,'measure',measure);
            if numel(plan)==8
                info.referenceY=a.y;
                info.referenceDy=a.dyMetric;
            end
        end

        function initial=reference(g)
            % Capture the fixed geometry and allocation used for later proposals.
            c=g.core;grid=g.gridOptions;
            stride=1;factor=1;
            if isfield(g,'axisSource')&&strcmp(g.axisSource.kind,'portable_preset')
                % Full radial map used to construct the reference mesh in the paper.
                stride=g.axisSource.stride;factor=2;
                grid=struct('nr',385,'ntheta',numel(g.polar.theta),'rmax',g.polar.r(end), ...
                    'firstSpacing',.03,'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01);
            elseif strcmp(grid.radialMap,'explicit_axes')
                grid.radialMap='uniform_core';
                grid.firstSpacing=min(min(g.polar.drMetric),grid.rmax/(4*(grid.nr-1)));
                grid.coreRadius=min(max(2,2*grid.firstSpacing),grid.rmax/4);
            elseif isfield(g,'radialRedistribution')
                factor=g.radialRedistribution.factor;
            end
            if isfield(c,'dxMetric')
                spacing=[c.dxMetric(c.originIndex),c.dyMetric(1)];
            else
                spacing=[c.dx,c.dy];
            end
            initial=struct('bounds',[c.x(end),c.y(end)], ...
                'coreMap',struct('intervals',[c.nx-1,c.ny-1],'minSpacing',spacing), ...
                'theta',g.polar.theta,'radialGrid',grid,'radialStride',stride,'tailFactor',factor);
            if isfield(g,'axisSource')&&isfield(g.axisSource,'referenceTheta')
                initial.theta=g.axisSource.referenceTheta;
            end
            if isfield(c,'stair')
                initial.stairAxes=struct('x',c.x,'y',c.y,'dxMetric',c.dxMetric,'dyMetric',c.dyMetric);
                if isfield(g.axisSource,'interfaceCore')
                    initial.stairAxes=g.axisSource.interfaceCore.referenceAxes;
                    initial.bounds=[initial.stairAxes.x(end),initial.stairAxes.y(end)];
                end
                if isfield(g.axisSource,'interfaceAllocation')
                    allocation=g.axisSource.interfaceAllocation;
                    initial.stairAxes.x=allocation.referenceX;
                    initial.stairAxes.dxMetric=allocation.referenceDx;
                    if isfield(allocation,'referenceY')
                        initial.stairAxes.y=allocation.referenceY;
                        initial.stairAxes.dyMetric=allocation.referenceDy;
                    end
                end
            end
        end

        function [axes,info]=propose(extents,memory)
            % Propose new patch axes from extents relative to the original mesh.
            initial=memory.initialMesh;
            ratio=extents./memory.initialExtents;
            bounds=initial.bounds.*ratio;
            if isfield(initial,'stairAxes')
                core=initial.stairAxes;
                core.x=ratio(1)*core.x;core.dxMetric=ratio(1)*core.dxMetric;
                core.y=ratio(2)*core.y;core.dyMetric=ratio(2)*core.dyMetric;
            else
                core=ipm_cr.mesh.axes.tensor(bounds(1),bounds(2),initial.coreMap);
            end

            grid=initial.radialGrid;
            innerScale=min([ratio,1]);
            grid.firstSpacing=grid.firstSpacing*innerScale;
            if isfield(grid,'coreRadius'),grid.coreRadius=grid.coreRadius*innerScale;end
            full=ipm_cr.mesh.axes.polar(grid);
            r=full.r;dr=full.radialMetric/(full.nr-1);
            if initial.tailFactor>1
                preserve=max(4*innerScale,2.5*hypot(bounds(1),bounds(2)));
                transition=min(grid.rmax,max(20*innerScale,4*preserve));
                [r,dr]=ipm_cr.mesh.axes.redistributeRadial(full,struct('enabled',true,'factor',initial.tailFactor, ...
                    'preserveRadius',preserve,'fullRefinementRadius',transition));
            end
            % Keep the common outer endpoint when halving the regenerated full axis.
            stride=initial.radialStride;keep=fliplr(numel(r):-stride:1);
            r=r(keep);dr=stride*dr(keep);
            requested=.6*min(bounds);
            first=find(r>=requested,1);
            if isfield(initial,'stairAxes')
                % Let the geometry constructor align the active staircase and annulus.
                first=find(r>0,1);requested=NaN;
            end
            axes=struct('x',core.x,'y',core.y,'r',r(first:end),'theta',initial.theta, ...
                'dxMetric',core.dxMetric,'dyMetric',core.dyMetric,'drMetric',dr(first:end));
            info=struct('initialBounds',initial.bounds,'initialExtents',memory.initialExtents, ...
                'measuredExtents',extents,'bounds',bounds,'requestedInnerRadius',requested, ...
                'innerRadius',axes.r(1),'precedingRadius',r(max(1,first-1)), ...
                'outerRadius',r(end),'radialStride',stride,'fullRadialNodes',numel(r));
        end

        function new=rebuild(old,axes)
            % Rebuild patches and their exchange plans from the proposed axes.
            geometry=struct('axes',axes,'ghostWidth',old.ghostWidth,'interpolationWidth',old.interpolationWidth);
            if isfield(old.core,'stair')
                geometry.staircase=true;
                geometry.staircaseFine=old.axisSource.staircaseFine;
                if isfield(old.axisSource,'k')
                    geometry.k=old.axisSource.k;
                elseif isfield(old.axisSource,'nodeCounts')
                    geometry.nodeCounts=old.axisSource.nodeCounts;
                end
                if isfield(old.axisSource,'radialRatioLimit')
                    geometry.radialRatioLimit=old.axisSource.radialRatioLimit;
                end
                if isfield(old.axisSource,'interfaceCore')
                    geometry.corePlateauIntervals=old.axisSource.interfaceCore.coreIntervals;
                end
                if isfield(old.axisSource,'interfaceAllocation')
                    geometry.interfaceAllocation=old.axisSource.interfaceAllocation.measure.plan;
                end
            end
            new=ipm_cr.mesh.build(geometry);
        end
    end
end

function a=resize_axes(a,counts)
% Resample the preset's Cartesian maps in their existing logical coordinates.
if isfield(counts,'nx')&&counts.nx~=numel(a.x)
    assert(mod(counts.nx,2)==1,'ipm_cr:MeshNodeCounts', ...
        'nx must be odd to retain symmetry and an exact Cartesian origin.');
    center=find(a.x==0,1);n=numel(a.x)-center;
    [x,dx]=hermite(a.x(center:end),a.dxMetric(center:end),linspace(0,n,(counts.nx+1)/2));
    dx=dx*(2*n/(counts.nx-1));
    a.x=[-fliplr(x(2:end)),x];a.dxMetric=[fliplr(dx(2:end)),dx];
end
if isfield(counts,'ny')&&counts.ny~=numel(a.y)
    n=numel(a.y)-1;
    [y,dy]=hermite(a.y,a.dyMetric,linspace(0,n,counts.ny));
    a.y=y(:);a.dyMetric=dy(:)*(n/(counts.ny-1));
end
% The angular and radial budgets are applied outside the Cartesian fit.
fields={'x','y','r','theta','dxMetric','dyMetric','drMetric'};
a=rmfield(a,setdiff(fieldnames(a),fields));
end

function [a,info]=interface_core(a,halo,coreIntervals)
% Keep the dense-region anchors and interval budget. Fit the log metric
% across a central plateau, a C2 ramp and a seam plateau.
reference=struct('x',a.x,'y',a.y,'dxMetric',a.dxMetric,'dyMetric',a.dyMetric);
validateattributes(coreIntervals,{'numeric'},{'vector','numel',2,'integer','positive'});
m=find(a.x==0,1);xp=a.x(m:end);
cap=min(xp(end-halo),a.y(end-halo));
limit=1.1+.1*double(numel(a.theta)==91);
[x,dx,xInfo]=fit_axis(xp,coreIntervals(1),cap,halo,limit);
[y,dy,yInfo]=fit_axis(a.y,coreIntervals(2),cap,halo,limit);
a.x=[-fliplr(x(2:end)),x];a.dxMetric=[fliplr(dx(2:end)),dx];
a.y=y(:);a.dyMetric=dy(:);
fields={'x','y','r','theta','dxMetric','dyMetric','drMetric'};
a=rmfield(a,setdiff(fieldnames(a),fields));
info=struct('coreIntervals',coreIntervals,'referenceRadius',cap,'adjacentStepRatioLimit',limit, ...
    'x',xInfo,'y',yInfo,'referenceAxes',reference);
end

function [x,metric,info]=fit_axis(v,inner,radius,halo,limit)
v=v(:).';n=numel(v)-1;
assert(inner<n,'ipm_cr:StairCoreSpacing','The dense-region quota must leave exterior intervals.');
% Long unused Cartesian extensions cost too many intervals under a
% spacing-ratio bound. Keep the physical donor box and trim excess padding.
outer=min(v(end),max(1.28*radius,1.08*v(inner+1)));
[x,h,centralEnd,rampEnd]=fit_head(v(inner+1),inner,outer,n,halo,limit);
metric=h(1:2:end).';
steps=diff(x);
adjacent=steps(2:end)./steps(1:end-1);
maxRatio=max([adjacent,1./adjacent]);
info=struct('centralBand',[0,x(centralEnd+1)], ...
    'seamBand',x([rampEnd+1,end]), ...
    'centralStep',h(1),'seamStep',h(end), ...
    'centralCells',centralEnd,'seamCells',[rampEnd+1,n], ...
    'maxAdjacentStepRatio',maxRatio);
end

function [x,h,centralEnd,rampEnd]=fit_head(coreEnd,inner,seam,seamIndex,halo,limit)
centralEnd=max(1,floor(.3*inner));
rampEnd=seamIndex-max(halo+1,round(.20*seamIndex));
assert(0<inner&&inner<seamIndex&&centralEnd<rampEnd&&seam>coreEnd, ...
    'ipm_cr:StairCoreSpacing','The Cartesian core and seam anchors must be ordered.');
t=(0:2*seamIndex)'/2;
z=max(0,min(1,(t-centralEnd)/(rampEnd-centralEnd)));
% Smooth the first and last fifth of a constant log-slope ramp.
shoulder=.2;
u=min(z/shoulder,1);
q=min((1-z)/shoulder,1);
endIntegral=@(s) shoulder*s.^3.*(1-s/2);
blend=(endIntegral(u)+max(z-shoulder,0)- ...
    (shoulder/2-endIntegral(q)))./(1-shoulder);
ratio=seam/coreEnd;
% The maximum derivative of blend is 1/((1-shoulder)*rampLength).
% Thus bounding the log metric derivative also bounds integrated steps.
changeLimit=min(32,(1-shoulder)*(rampEnd-centralEnd)*log(limit-.001));
minimum=integrated_ratio(-changeLimit,blend,inner);
maximum=integrated_ratio(changeLimit,blend,inner);
anchorReachable=minimum<ratio&&ratio<maximum;
if anchorReachable
    logChange=fzero(@(b) integrated_ratio(b,blend,inner)-ratio,[-32,32]);
else
    % Keep the tensor endpoints and point allocation at the slope limit.
    logChange=(2*double(ratio>=maximum)-1)*changeLimit;
end
h=exp(logChange*blend);
if anchorReachable
    h=h*(coreEnd/sum(integrate_cells(h(1:2*inner+1))));
else
    h=h*(seam/sum(integrate_cells(h)));
end
x=[0;cumsum(integrate_cells(h))].';
if anchorReachable,x(inner+1)=coreEnd;end
x(seamIndex+1)=seam;
end

function value=integrated_ratio(b,blend,inner)
h=integrate_cells(exp(b*blend));
value=sum(h)/sum(h(1:inner));
end

function [a,info]=core(a,scale,fine)
% The reference dense rectangle ends at 35% of each positive logical axis.
% A scalar scales the 2D share isotropically. A pair scales the X1 and X2
% interval shares independently. Each axis keeps its total interval budget.
cut=.35;
if isscalar(scale)
    validateattributes(scale,{'numeric'},{'real','finite','>=',.1,'<=',1});
    directionScale=sqrt(scale)*[1,1];
else
    validateattributes(scale,{'numeric'},{'vector','numel',2,'real','finite','>=',sqrt(.1),'<=',1.8});
    directionScale=reshape(scale,1,2);
end
extra=1.35*double(fine);baseTotal=1+extra*(1-cut);
m=(numel(a.x)+1)/2;
[x,dx,xcut]=redistribute(a.x(m:end),a.dxMetric(m:end),directionScale(1));
[y,dy,ycut]=redistribute(a.y,a.dyMetric,directionScale(2));
a.x=[-fliplr(x(2:end)),x];a.dxMetric=[fliplr(dx(2:end)),dx];
a.y=y(:);a.dyMetric=dy(:);
info=struct('targetRelativeCoreShare',prod(directionScale),'coreBounds',[xcut,ycut], ...
    'targetAxisScale',directionScale, ...
    'referenceLogicalCut',cut,'referenceAxisShare',cut/baseTotal, ...
    'targetAxisShare',directionScale*cut/baseTotal);

    function [v,d,physicalCut]=redistribute(v,d,axisScale)
        n=numel(v)-1;newN=round(n*baseTotal);
        if axisScale<=1
            outer=baseTotal*(1-axisScale)/(axisScale-baseTotal*cut^2);
            total=baseTotal+outer;
        else
            % Add a positive monitor over the core and taper it to zero.
            % This avoids driving the exterior node density toward zero.
            share=axisScale*cut/baseTotal;
            coreMass=1.5*cut;
            gain=(axisScale-1)*cut/(cut-share*coreMass);
            total=baseTotal+gain*coreMass;
        end
        target=linspace(0,total,newN+1);
        lo=zeros(size(target));hi=ones(size(target));
        for k=1:46
            s=(lo+hi)/2;
            mass=s+extra*max(0,s-cut).^3/(1-cut)^2;
            if axisScale<=1
                mass=mass+outer*s.^3;
            else
                mass=mass+gain*core_weight(s,cut);
            end
            left=mass<target;lo(left)=s(left);hi(~left)=s(~left);
        end
        s=(lo+hi)/2;s([1,end])=[0,1];
        [physicalCut,~]=hermite(v,d,n*cut);
        [v,d]=hermite(v,d,n*s);
        density=1+3*extra*max(0,s-cut).^2/(1-cut)^2;
        if axisScale<=1
            density=density+3*outer*s.^2;
        else
            [~,weight]=core_weight(s,cut);
            density=density+gain*weight;
        end
        d=d.*(n*total/newN)./density;
    end
end

function [mass,weight]=core_weight(s,cut)
t=max(0,min(1,(s-cut)/cut));
mass=s-cut*t.^4.*(2.5-3*t+t.^2)-max(0,s-2*cut);
weight=1-t.^3.*(10-15*t+6*t.^2);
end

function [r,dr]=radial(r,dr,interfaceRadius)
% Double resolution near the interface with a smooth transition to the
% original far axis. In reversed indices the enriched region is at the end.
n=numel(r)-1;
first=find(r>=1.18*interfaceRadius,1);
stop=max(1,n+1-first);start=max(0,stop-24);
span=stop-start;
added=ceil(n-(start+stop)/2);gain=added/(n-(start+stop)/2);
target=0:n+added;
lo=zeros(size(target));hi=n+lo;
for k=1:52
    z=(lo+hi)/2;
    t=max(0,min(1,(z-start)/span));
    mass=z+gain*(span*t.^4.*(2.5-3*t+t.^2)+max(0,z-stop));
    left=mass<target;lo(left)=z(left);hi(~left)=z(~left);
end
z=(lo+hi)/2;z(target<=start)=target(target<=start);z(end)=n;
t=max(0,min(1,(z-start)/span));
[rev,d]=hermite(-flip(r(:)),flip(dr(:)),z);
d=d./(1+gain*t.^3.*(10-15*t+6*t.^2));
r=-flip(rev(:));dr=flip(d(:));
end

function [r,dr,info]=match_radial(r,dr,core,cover,angle,halo,ratioLimit,requestedNr)
if nargin<8,requestedNr=[];end
% Along a ray, cos(theta)/dx + sin(theta)/dy counts Cartesian grid-line
% crossings per unit radius. Match the median reciprocal at the interface.
co=abs(cos(angle));si=sin(angle);
ix=discretize(cover*co,core.x);iy=discretize(cover*si,core.y);
dx=diff(core.x(:));dy=diff(core.y(:));
cartesianStep=1./(co./dx(ix)+si./dy(iy));
baseIndex=find(r<=cover,1,'last');
baseStep=r(baseIndex+1)-r(baseIndex);
h=min(median(cartesianStep),baseStep);
% Keep the first active ring inside the last full Cartesian cells.
% A radius equal to cover is already inside that staircase by construction.
inner=cover;
uniformCells=6;
rampStart=inner+uniformCells*h;
uniform=inner+(-halo:uniformCells)'*h;
if ~isempty(requestedNr)
    % Fit the requested number of intervals to the same six-cell seam and
    % outer endpoint. The log-metric ramp has zero slope at the seam.
    [r,dr,maxRatio]=radial_budget(uniform,h,r(end),requestedNr);
    fullRatio=maxRatio;
    endRadius=r(end);
else
    [r,dr,maxRatio,fullRatio,endRadius]=join_radial(r,dr,uniform,h,ratioLimit);
end
ratio=h./cartesianStep;
info=struct('uniformActiveCells',uniformCells,'innerActiveRadius',inner, ...
    'baseRadialStep',baseStep,'radialStep',h,'medianCoreRadialStep',median(cartesianStep), ...
    'radialToCoreStepRange',[min(ratio),max(ratio)], ...
    'maxAdjacentStepRatio',maxRatio,'maxWholeAxisStepRatio',fullRatio, ...
    'adjacentStepRatioLimit',ratioLimit, ...
    'transitionEndRadius',endRadius, ...
    'coreStepAtAngle',cartesianStep(:).','angles',angle(:).');
end

function [r,dr,localRatio,wholeRatio,endRadius]=join_radial(baseR,baseDr,prefix,h,target)
% Keep the first good geometric join at the target when available. Over
% target, use the positive transition with the smallest whole-axis ratio.
rampStart=prefix(end);
steps=diff(baseR);
adjacent=steps(2:end)./steps(1:end-1);
suffixRatio=flipud(cummax(flipud([max(adjacent,1./adjacent);1;1])));
first=find(baseR>rampStart,1);
r=[];dr=[];wholeRatio=Inf;
margin=min(.005,.05*(target-1));
if ~isempty(first)
    for join=first:numel(baseR)
        span=baseR(join)-rampStart;
        logChange=log(baseDr(join)/h);
        average=h;
        if logChange~=0,average=h*expm1(logChange)/logChange;end
        cells=round(span/average);
        if cells<6,continue;end
        t=(0:2*cells)'/(2*cells);
        metric=h*exp(logChange*t);
        taper=6*t.*(1-t);
        initial=integrate_cells(metric);
        correction=integrate_cells(metric.*taper);
        metric=metric.*(1+(span-sum(initial))/sum(correction)*taper);
        if any(~isfinite(metric)|metric<=0),continue;end
        ramp=rampStart+cumsum(integrate_cells(metric));
        ramp(end)=[];
        nearSteps=diff([prefix;ramp;baseR(join:min(join+1,numel(baseR)))]);
        if any(nearSteps<=0),continue;end
        adjacent=nearSteps(2:end)./nearSteps(1:end-1);
        ratio=max([adjacent;1./adjacent]);
        totalRatio=max(ratio,suffixRatio(join));
        atTarget=ratio<target-margin&&suffixRatio(join)<target;
        if atTarget||totalRatio<wholeRatio
            r=[prefix;ramp;baseR(join:end)];
            dr=[h*ones(size(prefix));metric(3:2:end-2);baseDr(join:end)];
            localRatio=ratio;wholeRatio=totalRatio;endRadius=baseR(join);
        end
        if atTarget,break;end
    end
end
assert(~isempty(r),'ipm_cr:StaircaseRadius', ...
    'The radial axis must admit a positive transition between the interface band and the outer grid.');
end

function [r,dr,maxRatio]=radial_budget(prefix,h,outer,nr)
cells=nr-numel(prefix);span=outer-prefix(end);
assert(cells>=1&&span>0,'ipm_cr:MeshRadialBudget', ...
    'nr must leave cells beyond the ghost rows and six uniform interface cells.');
z=(0:2*cells)'/2;
shoulder=min(12,.15*cells);t=min(1,z/shoulder);
ramp=shoulder*t.^3.*(1-.5*t)+max(0,z-shoulder);
spanInSteps=span/h;
% The seam endpoint contributes one sixth of its step to the first cell.
% Bracket the metric fit from the physical span and count alone.
assert(isfinite(spanInSteps)&&spanInSteps>1/6,'ipm_cr:MeshRadialBudget', ...
    'The space beyond the interface band must support a positive radial transition.');
lower=min(0,log((spanInSteps-1/6)/cells)/ramp(2));
upper=(log(max(1,spanInSteps))+log(12))/ramp(end);
distance=@(slope) sum(integrate_cells(exp(slope*ramp)));
slope=fzero(@(q) distance(q)-spanInSteps,[lower,upper]);
metric=h*exp(slope*ramp);
r=[prefix;prefix(end)+cumsum(integrate_cells(metric))];
r(end)=outer;
dr=[h*ones(size(prefix));metric(3:2:end)];
steps=diff(r);
assert(all(steps>0)&&all(isfinite(dr)&dr>0),'ipm_cr:MeshRadialBudget', ...
    'The radial budget must produce strictly increasing nodes and positive finite metrics.');
ratios=steps(2:end)./steps(1:end-1);
maxRatio=max([ratios;1./ratios]);
end

function step=integrate_cells(metric)
% Simpson quadrature on each unit index interval.
step=(metric(1:2:end-2)+4*metric(2:2:end-1)+metric(3:2:end))/6;
end

function [v,d]=hermite(v,d,z)
% Interpolate physical coordinates using the supplied metric as index slope.
v=v(:).';d=d(:).';n=numel(v)-1;
j=min(floor(z),n-1)+1;t=z-(j-1);
gap=v(j+1)-v(j);m0=d(j);m1=d(j+1);
c2=3*gap-2*m0-m1;c3=m0+m1-2*gap;
first=v(1);last=v(end);
v=v(j)+t.*(m0+t.*(c2+t.*c3));v(z==0)=first;v(z==n)=last;
d=m0+t.*(2*c2+3*t.*c3);
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end

function value=geometric_logsum(number,logRatio)
% log(1+q+...+q^(number-1)), evaluated without a large intermediate sum.
if logRatio==0
    value=log(number);
else
    value=log_expm1(number*logRatio)-log_expm1(logRatio);
end
end

function value=radial_ratio_log(intervals,parameter,kind)
if strcmp(kind,'exp_gaps')
    value=geometric_logsum(intervals,parameter);
elseif parameter==0
    value=log(intervals);
else
    value=log_sinh(parameter)-log_sinh(parameter/intervals);
end
end

function value=log_sinh(argument)
% log(sinh(argument)), without overflowing sinh for a very large box.
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=log_expm1(argument)
value=argument+log(-expm1(-argument));
end

function value=monitor_integral(u,ua,ub)
value=zeros(size(u));upper=u>=ub;value(upper)=u(upper)-(ua+ub)/2;
inside=u>ua&u<ub;
if any(inside)
    t=(u(inside)-ua)/(ub-ua);
    value(inside)=(ub-ua)*(t.*betainc(t,8,8)-.5*betainc(t,9,8));
end
end

function value=redistributeRadial_log_sinh(argument)
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=redistributeRadial_option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end

function value=radialAxis_log_sinh(argument)
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=radialAxis_option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end

function [physical,J,metric,h,logical,s]=axis_map(extent,n,requested,symmetric)
span=1+double(symmetric);average=span*extent/n;h=span/n;
if symmetric,logical=linspace(-1,1,n+1);else,logical=linspace(0,1,n+1)';end
if average<=requested
    s=0;
    % Match the original physical LINSPACE exactly when the map degenerates.
    if symmetric,physical=linspace(-extent,extent,n+1);
    else,physical=linspace(0,extent,n+1)';end
    J=extent*ones(size(logical));
else
    target=log(average)-log(requested);lo=0;hi=1;
    while log_sinhc(hi)<target,hi=2*hi;end
    for k=1:100
        mid=(lo+hi)/2;
        if mid==lo||mid==hi,break;end
        if log_sinhc(mid)<target,lo=mid;else,hi=mid;end
    end
    s=(lo+hi)/2;
    if s<40
        physical=extent*sinh(s*logical)/sinh(s);
        J=(extent*s/sinh(s))*cosh(s*logical);
    else
        % Avoid overflow for valid positive spacing ratios with large alpha.
        q=abs(logical);denominator=-expm1(-2*s);
        positive=exp(s*(q-1));negative=exp(-s*(q+1));
        physical=extent*sign(logical).*(positive-negative)/denominator;
        J=(extent*s/denominator)*(positive+negative);
    end
end
if symmetric
    physical([1,end])=[-extent,extent];physical(n/2+1)=0;
else
    physical(1)=0;physical(end)=extent;
end
metric=J*h;
assert(all(isfinite(physical))&&all(diff(physical)>0)&&all(isfinite(J))&&all(J>0)&& ...
    all(isfinite(metric))&&all(metric>0),'ipm_cr:TensorMapRange', ...
    'The mapped axis or analytic Jacobian is not representable.');
assert(min(metric)<=requested*(1+64*eps),'ipm_cr:TensorMapSpacing', ...
    'The analytic origin spacing exceeds the requested cap.');
end

function value=log_sinhc(s)
% log(sinh(s)/s), including its removable origin singularity.
if s==0,value=0;
elseif s<.1
    z=s*s;
    value=log1p(z/6+z^2/120+z^3/5040+z^4/362880+z^5/39916800+z^6/6227020800);
elseif s<40,value=log(sinh(s)/s);
else,value=s-log(2*s)+log1p(-exp(-2*s));end
end

function [v,d]=redistribute(v,d,windows,delta)
% Integrate a C2 quintic cumulative allocation and invert in index space.
n=numel(v)-1;target=0:n+delta;
lo=zeros(size(target));hi=n+lo;
for iteration=1:47
    mid=(lo+hi)/2;
    [mapped,~]=mass(mid,windows);
    below=mapped<target;
    lo(below)=mid(below);hi(~below)=mid(~below);
end
z=(lo+hi)/2;
% Restore every integer sample outside the windows exactly, including caps.
anchor=0:n;[position,~]=mass(anchor,windows);
fixed=position==round(position);
inside=false(size(anchor));
for j=1:size(windows,1)
    inside=inside|(anchor>windows(j,1)&anchor<windows(j,2));
end
keep=fixed&~inside;
z(position(keep)+1)=anchor(keep);
[~,density]=mass(z,windows);
j=min(floor(z),n-1)+1;t=z-j+1;
gap=v(j+1)-v(j);m0=d(j);m1=d(j+1);
c2=3*gap-2*m0-m1;c3=m0+m1-2*gap;
samples=v(j)+t.*(m0+t.*(c2+t.*c3));
samples(z==n)=v(end);
d=(m0+t.*(2*c2+3*t.*c3))./density;
v=samples;
end

function [value,derivative]=mass(z,windows)
value=z;derivative=ones(size(z));
for j=1:size(windows,1)
    row=windows(j,:);width=row(2)-row(1);
    t=min(1,max(0,(z-row(1))/width));
    value=value+row(3)*t.^3.*(10-15*t+6*t.^2);
    derivative=derivative+row(3)/width*30*t.^2.*(1-t).^2;
end
end
