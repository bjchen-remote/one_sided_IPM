function peakSoFar=studyTau40(folder,label,trial,gBefore,state,gAfter,meshInfo,peakSoFar)
% Passive active-node residual audit within a requested time window.
tau=state.scale.canonicalTime;
before=maxima(trial,gBefore);after=maxima(state,gAfter);
fid=fopen(fullfile(folder,[label,'_peaks.csv']),'a');
assert(fid>=0);closeFile=onCleanup(@()fclose(fid));
if ftell(fid)==0
    fprintf(fid,['step,tau,remeshed,preCoreQ,preCoreX,preCoreY,prePolarQ,prePolarX,prePolarY,' ...
        'postCoreQ,postCoreX,postCoreY,postPolarQ,postPolarX,postPolarY\n']);
end
fprintf(fid,'%d,%.17g,%d',state.step,tau,meshInfo.remeshed);
fprintf(fid,',%.17g',[before,after]);fprintf(fid,'\n');
clear closeFile
peak=max(abs(after([1,4])));
if meshInfo.remeshed
    pre=frame(trial,gBefore);post=frame(state,gAfter);
    write_frame(folder,sprintf('%s_remesh_%06d.mat',label,state.step),struct('pre',pre,'post',post,'info',meshInfo));
end
name=[label,'_peak.mat'];path=fullfile(folder,name);
if isempty(peakSoFar)
    peakSoFar=0;
    if isfile(path)
        previous=load(path,'record');
        peakSoFar=previous.record.activePeak;
    end
end
if peak>peakSoFar
    peakSoFar=peak;
    record=frame(state,gAfter);record.activePeak=peak;
    write_frame(folder,name,record);
end
end

function values=maxima(s,g)
values=zeros(1,6);
for j=1:2
    charts={'core','polar'};name=charts{j};
    p=g.(name);q=s.rhs.Q.(name);
    [~,i]=max(abs(q(p.activeIndices)));idx=p.activeIndices(i);
    values(3*j+(-2:0))=[q(idx),p.X(idx),p.Y(idx)];
end
end

function a=frame(s,g)
a=struct('tau',s.scale.canonicalTime,'step',s.step,'geometry',g.geometryKey,'Q',s.rhs.Q);
end

function write_frame(folder,name,record)
pending=fullfile(folder,[name,'.writing']);final=fullfile(folder,name);
save(pending,'record','-v7');[ok,why]=movefile(pending,final,'f');assert(ok,'%s',why);
end
