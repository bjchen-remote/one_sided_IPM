function atomicJson(path,value)
temp=[path,'.writing'];f=fopen(temp,'w');assert(f>=0,'ipm_c:JsonFile','Validation failed.');
try,fwrite(f,jsonencode(value,PrettyPrint=true),'char');fwrite(f,newline,'char');fclose(f);
catch problem,fclose(f);rethrow(problem);end
[ok,msg]=movefile(temp,path,'f');assert(ok,'ipm_c:JsonRename','%s',msg);
end
