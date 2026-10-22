function folder=root()
%ROOT Project root for the installed source layout, never the current folder.
folder=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
