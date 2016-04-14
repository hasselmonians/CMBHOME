function add_epoch(db,session_id,labeln,epochn,notesn)
% Add an epoch. If nothing other than session_id passed, then add the
% default values for label & epoch.

import CMBHOME.Database.*


if ~exist('labeln','var')
    label = 'full session';
    epoch = '[0 inf]';
    notes = 'NULL';
else
    label = labeln;
    epoch = epochn;
    if exist('notes','var')
        notes = notesn;
    else
        notes = 'NULL';
    end
end

%%{
sv = sprintf('INSERT INTO epochs (label,session_id,epoch,notes) VALUES(''%s'',%s,%s,''%s'')', ...
    label, ...
    num2str(session_id), ...
    mat2str(epoch), ...
    notes');

sv = fns(sv);

db.prepareStatement(sv);
db.query();

end