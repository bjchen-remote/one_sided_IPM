classdef initial
    %INITIAL Immutable initial data, Gaussian evaluation and materialization.
    methods (Static)
        function rho=density(ops,physics)
            %INITIAL Registered Gaussian datum with an optional smooth radial cutoff.
            % The immutable wall authority admits this datum only. Keep its original
            % arithmetic: erfc avoids cancellation in 1+erf on the negative tail.
            assert(isfield(physics,'initialCondition')&& ...
                strcmp(physics.initialCondition,'localized_gaussian_primitive'), ...
                'ipm_cr:InitialCondition','The maintained solver uses the original Gaussian datum.');
            A=physics.gaussianAmplitude;a=physics.gaussianCenterX;
            sx=physics.gaussianWidthX;sy=physics.gaussianWidthY;
            w=physics.gaussianSkewWeight;z=(ops.X-a)/sx;
            primitive=(1-w)*erfc(-z)+w*erfc(-(z-.8)/1.25);
            rho=A*sx*sqrt(pi)/2*primitive.*exp(-(ops.Y/sy).^2);
            if isfield(physics,'gaussianCutoff')
                radii=physics.gaussianCutoff;
                validateattributes(radii,{'numeric'},{'real','finite','positive','increasing','vector','numel',2});
                rho=rho.*compact_cutoff(hypot(ops.X,ops.Y),radii(1),radii(2));
            elseif isfield(physics,'gaussianSourceCutoff')
                radii=physics.gaussianSourceCutoff;
                validateattributes(radii,{'numeric'},{'real','finite','positive','increasing','vector','numel',2});
                outside=abs(ops.X)>radii(1)|abs(ops.Y)>radii(1);
                if any(outside,'all')
                    model=struct('amplitude',A,'centerX',a,'widthX',sx,'widthY',sy,'skewWeight',w,'cutoff',radii);
                    value=ipm_cr.model.initial.compactGaussian(model,ops.X(outside),ops.Y(outside),false);
                    rho(outside)=value.rho;
                end
            end
        end

        function cfg=problem(initial)
            %PROBLEM Physical datum and fixed I/P/C reference declared in one manifest.
            % Both manifest and MAT paths are project-relative; neither may escape
            % through a symbolic link. The reference constructor validates the datum.
            if nargin<1,initial='data/initial.json';end
            folder=ipm_cr.io.identity('root');[manifestFile,path]=input_path(folder,initial);
            manifest=jsondecode(fileread(path));
            required={'format','file','sha256','physics','grid','transport','scaling'};
            assert(isstruct(manifest)&&isscalar(manifest)&&all(isfield(manifest,required))&& ...
                strcmp(manifest.format,'ipm_cr_initial_v1'), ...
                'ipm_cr:ProblemFormat','Use a complete ipm_cr_initial_v1 manifest.');
            [~,authorityFile]=input_path(folder,manifest.file);
            assert(((ischar(manifest.sha256)&&isrow(manifest.sha256))|| ...
                (isstring(manifest.sha256)&&isscalar(manifest.sha256)))&& ...
                ~isempty(regexp(char(manifest.sha256),'^[0-9a-fA-F]{64}$','once')), ...
                'ipm_cr:ProblemIdentity','The initial MAT SHA must contain 64 hexadecimal digits.');
            cfg=struct('physics',manifest.physics,'grid',manifest.grid,'transport',manifest.transport, ...
                'scaling',manifest.scaling,'authorityFile',authorityFile, ...
                'authoritySha',lower(char(manifest.sha256)), ...
                'manifestFile',manifestFile,'manifestSha',ipm_cr.io.identity('sha',path));
            assert(isfile(cfg.authorityFile)&&strcmp(ipm_cr.io.identity('sha',cfg.authorityFile),cfg.authoritySha), ...
                'ipm_cr:ProblemIdentity','The original initial datum is missing or changed.');
        end

        function v=compactGaussian(m,x,y,second)
            %COMPACTGAUSSIAN Smoothly continue the primitive to its horizontal limits.
            % The density can have a plateau. Its x derivative has compact support.
            [wx,wx1,wx2]=cutoff(x,m.cutoff,second);
            [wy,wy1]=cutoff(y,m.cutoff,false);
            sx=m.widthX;sy=m.widthY;w=m.skewWeight;A=m.amplitude;
            z=(x-m.centerX)/sx;z2=(z-.8)/1.25;
            amplitude=A*sx*sqrt(pi)/2;
            primitive=amplitude*((1-w)*erfc(-z)+w*erfc(-z2));
            tail=primitive;
            right=x>m.cutoff(1);
            tail(right)=-amplitude*((1-w)*erfc(z(right))+w*erfc(z2(right)));
            horizontal=wx.*primitive;
            horizontal(right)=2*amplitude+wx(right).*tail(right);
            E=exp(-(y/sy).^2);vertical=E.*wy;
            v=struct('rho',horizontal.*vertical, ...
                'correction',E.*((wy-1).*primitive+wy.*((wx-1).*tail)));
            e1=exp(-z.^2);e2=exp(-z2.^2);
            gx=A*((1-w)*e1+(w/1.25)*e2);
            horizontalX=wx.*gx+wx1.*tail;
            verticalY=E.*(wy1-2*(y/sy^2).*wy);
            v.x=horizontalX.*vertical;
            v.y=horizontal.*verticalY;
            if second
                gxx=(-2*A/sx)*((1-w)*(z.*e1)+(w/1.25^2)*(z2.*e2));
                v.xx=(wx.*gxx+2*wx1.*gx+wx2.*tail).*vertical;
                v.xy=horizontalX.*verticalY;
            end
        end

        function initial=materialize(c,rho0)
            %SET_INITIAL Materialize a density expression as an immutable input.
            % rho0=@(x,y)... must return matching finite real arrays (a scalar is allowed).
            root=ipm_cr.io.identity('root');
            [~,method]=ipm_cr.config(c,'run');g=ipm_cr.mesh.build(method.geometry);
            base=ipm_cr.model.initial.problem();old=load(base.authorityFile,'frame');
            assert(isa(rho0,'function_handle'),'ipm_cr:InitialFunction','Write rho0 as @(x,y) your_expression.');
            evaluate=rho0;description=func2str(rho0);
            [X,Y]=meshgrid(old.frame.x,old.frame.y);rho=values(X,Y);
            native=struct('geometryKey',g.geometryKey,'core',values(g.core.X,g.core.Y), ...
                'polar',values(g.polar.X,g.polar.Y),'origin',values(0,0));
            frame=struct('x',old.frame.x,'y',old.frame.y,'rho',rho, ...
                'step',0,'tau',0,'t',0,'Cx',1,'Cy',1,'Comega',1,'Xshift',0,'native',native); %#ok<NASGU>
            digest=java.security.MessageDigest.getInstance('SHA-256');
            digest.update(unicode2native([base.manifestSha,jsonencode(g.geometryKey)],'UTF-8'));
            for v={rho,native.core,native.polar,native.origin}
                digest.update(typecast(double(v{1}(:)),'uint8'));
            end
            hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
            name=['rho_',hash(1:24)];initial=['data/cases/',name,'/initial.json'];
            parent=fullfile(root,'data','cases');target=fullfile(parent,name);
            if isfolder(target)
                cfg=ipm_cr.model.initial.problem(initial);saved=load(cfg.authorityFile,'frame');
                assert(isequaln(saved.frame,frame),'ipm_cr:InitialIdentity','Existing initial data differ; preserve and inspect %s.',initial);
                return
            end
            if ~isfolder(parent),mkdir(parent);end
            stage=tempname(parent);mkdir(stage);cleanup=onCleanup(@()remove_stage(stage)); %#ok<NASGU>
            file=fullfile(stage,'t0.mat');save(file,'frame','-v7.3');
            manifest=jsondecode(fileread(fullfile(root,'data','initial.json')));
            manifest.physics=struct('initialCondition','sampled_density','symmetryMode','half_plane');
            manifest.file=['data/cases/',name,'/t0.mat'];manifest.sha256=ipm_cr.io.identity('sha',file);
            manifest.initialDescription=description;
            manifest.representation='Exact materialized native samples; transport mesh identity is pinned. No function handle required on restart.';
            for key={'referenceWallIntegral','referenceWallGradientCentroid','referenceWallPeak'}
                if isfield(manifest.scaling,key{1}),manifest.scaling=rmfield(manifest.scaling,key{1});end
            end
            cfg=base;cfg.physics=manifest.physics;cfg.scaling=manifest.scaling;
            cfg.authorityFile=file;cfg.authoritySha=manifest.sha256;
            [~,~,refs]=ipm_cr.model.reference(cfg);
            assert(all(isfinite([refs.I0,refs.P0,refs.C0]))&&refs.I0>0&&refs.P0>0, ...
                'ipm_cr:InitialGauge','This initial density does not admit the configured I/P/C normalization.');
            manifest.scaling.referenceWallIntegral=refs.I0;
            manifest.scaling.referenceWallGradientCentroid=refs.C0;
            ipm_cr.io.identity('json',fullfile(stage,'initial.json'),manifest);
            java.nio.file.Files.move(java.io.File(stage).toPath(),java.io.File(target).toPath(), ...
                javaArray('java.nio.file.CopyOption',0));
            fprintf('Initial density saved: %s\n',initial);

                function R=values(X,Y)
                    R=evaluate(X,Y);
                    if isnumeric(R)&&isscalar(R),R=repmat(R,size(X));end
                    assert(isnumeric(R)&&isreal(R)&&isequal(size(R),size(X))&&all(isfinite(R),'all'), ...
                        'ipm_cr:InitialFunction','rho0(x,y) must return matching finite real values.');
                    R=double(R);
                end

        end
    end
end

function weight=compact_cutoff(x,start,stop)
%COMPACTCUTOFF Smooth cutoff equal to one before START and zero after STOP.
u=(x-start)/(stop-start);
weight=double(u<=0);
inside=u>0&u<1;
z=u(inside);
rise=exp(-1./z);
fall=exp(-1./(1-z));
weight(inside)=fall./(rise+fall);
end

function [relative,path]=input_path(root,value)
assert((ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)), ...
    'ipm_cr:ProblemPath','Use a project-relative input path.');
relative=char(value);parts=strsplit(relative,'/','CollapseDelimiters',false);
assert(~isempty(relative)&&~startsWith(relative,'/')&&~contains(relative,'\')&& ...
    ~contains(relative,':')&&~any(ismember(parts,{'','..','.'})), ...
    'ipm_cr:ProblemPath','Use a canonical project-relative input path.');
base=char(java.io.File(root).getCanonicalPath());
path=char(java.io.File(fullfile(root,parts{:})).getCanonicalPath());
assert(startsWith(path,[base,filesep]),'ipm_cr:ProblemPath','An initial input resolves outside the project.');
end

function [w,d1,d2]=cutoff(x,radii,second)
start=radii(1);stop=radii(2);width=stop-start;
z=(abs(x)-start)/width;
w=double(z<=0);d1=zeros(size(z));d2=[];
inside=z>0&z<1;
u=z(inside);v=1-u;
rise=exp(-1./u);fall=exp(-1./v);total=rise+fall;
a=fall./total;b=rise./total;
L1=1./u.^2+1./v.^2;
w(inside)=a;
d1(inside)=-sign(x(inside)).*(a.*b).*L1/width;
if second
    d2=zeros(size(z));
    L2=-2./u.^3+2./v.^3;
    d2(inside)=(a.*b).*((1-2*a).*L1.^2-L2)/width^2;
end
end

function remove_stage(folder)
if isfolder(folder),rmdir(folder,'s');end
end
