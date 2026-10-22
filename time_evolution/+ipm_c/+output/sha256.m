function value=sha256(path)
%SHA256 Streaming file identity, without reading a large checkpoint at once.
d=java.security.MessageDigest.getInstance('SHA-256');
f=fopen(path,'rb');assert(f>=0,'ipm_c:HashFile','Validation failed.');cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
while true
    bytes=fread(f,1048576,'*uint8');if isempty(bytes),break;end;d.update(bytes);
end
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2)',1,[]));
end
