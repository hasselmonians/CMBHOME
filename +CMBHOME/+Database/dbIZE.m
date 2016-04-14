function out = dbIZE(fp)
% Turns a filepath into a filepath relative to the dropbox folder.

if ~iscell(fp)
    fp = CMBHOME.Utils.path2array(fp);
end

out = [];
for l = 1:length(fp)
    out = [out '/' fp{l}];
end


end