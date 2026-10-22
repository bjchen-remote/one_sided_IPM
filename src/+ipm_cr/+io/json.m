function json(path,value)
%JSON Replace a JSON record atomically in its containing directory.
temporary=[char(path),'.writing'];fid=fopen(temporary,'w');
assert(fid>=0,'ipm_cr:JsonOpen','Cannot open JSON staging file.');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(value,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,path,'f');
assert(ok,'ipm_cr:JsonRename','%s',message);
end
