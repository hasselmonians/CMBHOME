function cell_id = add_cell(db,tetrode_id,cut_id,subObj,sid)
% cell_id = add_cell(db,grp,rat,session,root.cel(i)(1),root.cel(i)(2))
%
% Adds a new cell to the CMBdatabase. If more than one input argument, uses
% the provided information to populate the cell information. 

% Bill 2012.08.30 

import CMBHOME.Database.*

%% Cell Statement
cellstate = sprintf('INSERT INTO cell (tetrode_id,type,notes,morphology) VALUES(''%s'',''%s'',''%s'',''%s'')',...
                    num2str(tetrode_id),...
                    subObj.type,...
                    subObj.notes,...
                    subObj.morphology);

                
cellstate = CMBHOME.Database.fns(cellstate);
 
db.prepareStatement(cellstate);
[~] = db.query();

db.prepareStatement('SELECT LAST_INSERT_ID()');
[a] = db.query();
cell_id = (a.LAST_INSERT_ID__);

%% Cell Session Statement
csstate = sprintf('INSERT INTO cell_session (cell_id,session_id,cut_id) VALUES(''%s'',''%s'',''%s'')',num2str(cell_id),num2str(sid),num2str(cut_id));
db.prepareStatement(csstate);
[~] = db.query();
    

