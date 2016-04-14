function link_cells(db,id1,id2)
% 
% This is a script created to merge two cell_id's into a single entry. The
% later of the two id's will be deleted (all references to it will be
% changed to the earlier id. 

% Bill 2012.08.29

dbpe
import CMBHOME.Database.*

%% Check to make sure that they are on the same tetrode
state = sprintf('SELECT tetrode_id FROM cell WHERE id = ''%s''',num2str(id1));
db.prepareStatement(state);
[a] = db.query;
tid1 = a.tetrode_id(1);

state = sprintf('SELECT tetrode_id FROM cell WHERE id = ''%s''',num2str(id2));
db.prepareStatement(state);
[a] = db.query;
tid2 = a.tetrode_id(1);

if tid1 ~= tid2
    fprintf('ERROR! Cells are not on the same tetrode \n');
    return
else
    
    % get lower of the two cell ids
    cid = num2str(min([id1 id2]));
    cid_bad = num2str(max([id1 id2]));
    
    % Change cell_session references of the higher id
    try
        state2 = sprintf('UPDATE cell_session SET cell_id = ''%s'' WHERE cell_id = ''%s'' ',cid,cid_bad);
        db.prepareStatement(state2);
        [a2] = db.query();
    catch
        fprintf('These cells are in the same session, they can''t be the same cell! \n');
        return
    end
    
    % Change cell_epoch references of the higher id
    state1 = sprintf('UPDATE cell_epoch SET cell_id = ''%s'' WHERE cell_id = ''%s'' ',cid,cid_bad);
    db.prepareStatement(state1);
    [a1] = db.query();
    
    
    % Delete the cell entry of the higher id
    state3 = sprintf('DELETE from cell USING cell WHERE id = ''%s'' ',cid_bad);
    db.prepareStatement(state3);
    [a3] = db.query();
    
end
