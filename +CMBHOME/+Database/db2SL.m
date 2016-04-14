function SL = db2SL(db,state)
% Accepts a database statement as the input (especially from GUI) and creates 
% an SL from that. Initial fields for SL will be the the session path,
% [trode cell],fpath_CMB,and any fields check in the GUI.

% wchapman 2013.01.14

%% Initialized with Base information
import CMBHOME.Database.*

db.prepareStatement(state);
q = db.query();
fnames = id2fname(db,q);

[trode cut] = id2cut(db,q);

SL(length(fnames)).fname = []; %initialize

for i = 1:length(SL)
    SL(i).fname = osIZE(fnames{i});
    SL(i).cel = [trode(i) cut(i)];
end

%% Add fields that the user asked for
if isfield(q,'filepath_cmb')
    q = rmfield(q,'filepath_cmb');
end

fns = fieldnames(q);
for fn = 1:length(fns)
    for i = 1:length(SL)
        eval(['SL(i).' fns{fn} ' = q.' fns{fn} '(i);'])
    end
end

SL = SL(:);
end