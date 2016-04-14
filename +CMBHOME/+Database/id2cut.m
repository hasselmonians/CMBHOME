function [cut] = id2cut(db,q)
% root.cel = helpful.id2cut(db,cell_id,session_id);
%
% Takes in the database object, cell_id & session_id to return the "cut"
% identifications.

% wchapman 2013.01.14

import CMBHOME.Database.*

cids = q.cell_id;
sids = q.session_id;

for i = 1:length(cids)
% The cut_id:
    db.prepareStatement('SELECT cut_id FROM cell_session WHERE (cell_id = "{S}" AND session_id = "{S}")',num2str(cids(i)),num2str(sids(i)));
    a = db.query();
    cut_id = a.cut_id;

    % The trode:
    db.prepareStatement('SELECT trode from tetrode WHERE id IN (SELECT tetrode_id FROM cell WHERE (id = "{S}"))',num2str(cids(i)));
    a = db.query();
    trode = a.trode;

    cut{i} = [trode cut_id];
end