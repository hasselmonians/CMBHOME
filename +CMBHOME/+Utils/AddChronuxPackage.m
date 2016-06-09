function AddChronuxPackage
% Add the chronux pacage
if ~isempty(findstr('chronux', path)), return; end % check that it isnt already there

a = what('CMBHOME');

if length(a) > 1
    addpath(path, genpath(fullfile(a(2).path(1:end-8), 'chronux'))); % adds chronux package
else 
    addpath(path, genpath(fullfile(a.path(1:end-8), 'chronux'))); % adds chronux package
end