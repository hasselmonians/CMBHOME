function tetrode_id = add_tetrode(db,rat_id,trode,notes)
% Adds a tetrode for a given rat to the database, if it already exists in
% the database (identified by the same trode number in an existing rat),
% then returns the tetrode_id of the existing trode.  

import CMBHOME.Database.*

if ~isstr(trode)
    trode = num2str(trode);
end

if ~isstr(rat_id)
    rat_id = num2str(rat_id);
end

if ~exist('notes','var')
    notes = 'NULL';
end

check_state = sprintf('SELECT id FROM tetrode WHERE (rat_id = ''%s'' AND trode = ''%s'')',rat_id,trode);
db.prepareStatement(check_state);
a = db.query();

if isempty(a.id)

    tetrodestate = sprintf('INSERT INTO tetrode (trode,rat_id,notes) VALUES(''%s'',''%s'',''%s'')',trode,rat_id,notes);
    tetrodestate = fns(tetrodestate);
    
    db.prepareStatement(tetrodestate,10001);
    [a] = db.query();

    tetrode_id_state = sprintf('SELECT id FROM tetrode WHERE (rat_id = ''%s'' AND trode = ''%s'') ',rat_id,trode);
    db.prepareStatement(tetrode_id_state,10001);
    a = db.query();
    tetrode_id = num2str(a.id);

else
    fprintf('Tetrode Already exists, returning existing ID \n');
    tetrode_id = num2str(a.id);
end



end