function pth = osIZE(pth)
% path_var = CMBHOME.Database.osIZE(path_var);
%
% Takes a path from the database and formats it to work with your local
% setup. 

%pth = [dropboxPath  '/' pth];
pth = strrep(pth,'dropboxPath',dropboxPath);
pth = strrep(pth,'//','/');

if ispc == 1
    pth = strrep(pth,'/','\');
end 

end