function json(path,value)
%JSON Replace a UTF-8 JSON record atomically in its containing directory.
temporary=[char(path),'.writing'];fid=fopen(temporary,'wb');
assert(fid>=0,'ipm_cr:JsonOpen','Cannot open JSON staging file.');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,unicode2native(jsonencode(value,'PrettyPrint',true),'UTF-8'),'uint8');clear cleanup
[ok,message]=movefile(temporary,path,'f');
assert(ok,'ipm_cr:JsonRename','%s',message);
end
