function fnames = id2fname(db,q)
% fnames = CMBHOME.Database.id2fname(db,q)
%
% Takes in a struct with field q.session_id and returns the corresponding
% filepath_cmb

% wchapman 2013.01.14

sids = q.session_id;

for i = 1:length(sids)
    db.prepareStatement('SELECT filepath_cmb FROM session WHERE id = "{S}"',num2str(sids(i)));
    q = db.query();
    fnames{i} = q.filepath_cmb;
end


end