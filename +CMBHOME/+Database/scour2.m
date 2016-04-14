function scour2(path,db,grp,aquisition_system)
% scour(path_to_search)
%
% Looks through a folder and all of it's subfolders for CMBObjects that
% have not yet been added to the database. When it finds one, adds the rat,
% session, tetrode,tetrode_session, cell, and cell_session if they don't
% already exist. 
addpath([dropboxPath '/Bill/CMBHOME_beta'])
import CMBHOME.Database.*

rat_folders = CollectFolders(path); 

for i = 1:length(rat_folders)
    i
    ratName = rat_folders{i};
    ratID = add_rat(db,ratName,grp);
    
    session_folder = CollectFolders([path '/' rat_folders{i}]);
    
    for p = 1:length(session_folder)
        p
        session_folder{p} = [path '/' rat_folders{i} '/' session_folder{p}];
        cmbsv = dir(fullfile(session_folder{p},'*CMB*.mat'));
    
        for k = 1:length(cmbsv)
            cmbs{k} = cmbsv(k).name;
        end

        if ~exist('cmbs','var')
            continue
        end

        % Add the session #### NEEDS TO CHECK TO SEE IF ALREADY IN####
        for k = 1:length(cmbs)
            cmbs{k} = [session_folder{p} '/' cmbs{k}];
            %cmbs{k} = strrep(cmbs{i},'unitrecordingdata/','');
            cmbs{k} = CMBHOME.Database.dbIZE(cmbs{k});
            sid = CMBHOME.Database.add_session(db,ratName,cmbs{k},aquisition_system);  

            % Load the CMBobject 
            fp = CMBHOME.Database.osIZE(cmbs{k});
            c = load(fp);
            root = c.root;
            clear c; 

            % Add tetrodes if needed & add cell while in there:
            if length(root.cells > 0)
                for m = 1:size(root.cells,1)
                    tid{m} = add_tetrode(db,ratID,root.cells(m,1),'NULL'); % Tetrode
                    cid{m} = add_cell(db,tid{m},sid,root.cells(m,2));
                end
            end

        end

        clear cmbs
        clear cmbsv
    end
    
    clear session_folder

end