function ratid = add_rat(db,subObj)
% Takes in a database object and possibly varargin to create a new rat in
% the CMB database. 
%
% If more than one input, then uses varargin to fill the fields, if only
% one field, then prompts a GUI. If the rat already exists, then inform the
% user and simply return the id for the existing rat. 

% Bill 2012.08.30

import CMBHOME.Database.*


%% Check to see if it's already in the database

check_state = sprintf('SELECT id FROM rat WHERE (name = ''%s'' AND grp_name = ''%s'')',subObj.ratname,subObj.grp_name);
db.prepareStatement(check_state);
a = db.query();

if isempty(a.id)

    ratstate = sprintf('INSERT INTO rat (name,grp_name) VALUES(''%s'',''%s'')',subObj.ratname,subObj.grp_name);

    db.prepareStatement(ratstate,10001);
    [~] = db.query();

    rat_id_state = sprintf('SELECT id FROM rat WHERE (name = ''%s'' AND grp_name = ''%s'')',subObj.ratname,subObj.grp_name);
    db.prepareStatement(rat_id_state,10001);
    a = db.query();
    ratid = num2str(a.id);

else
    fprintf('Your rat already exists, returning the ID for it \n');
    ratid = num2str(a.id);

end

end