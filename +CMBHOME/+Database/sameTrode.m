function op = sameTrode(db,inv)
% Finds all cells that are on the same tetrode

% Inputs: inv: If CMBHOME.Session, then uses root.name_formatted to look up
%              If query result, uses inv.cell_id (resorts to inv.id if needed) 
% Outputs: op -- op.cid = cell.id of all cells on the tetrode

if strcmp(class(inv),'CMBHOME.Session')
    
    trode = inv.cel(1);
    fname_cmb = CMBHOME.Database.dbIZE(inv.name_formatted);
    
    state = ['SELECT rat_id FROM session WHERE filepath_cmb = ''' fname_cmb ''''];
    state = ['SELECT id FROM tetrode WHERE rat_id in (' state ')'];
    state = ['SELECT id FROM cell WHERE tetrode_id IN (' state ')'];
    
    db.prepareStatement(state);
    q = db.query();
    
elseif isstruct('inv')
    
    
    
else
   error('Unknown input type') 
end


end