function align_cells(db,dbp,rat_id)
% Creates a temporary .mat file of CMBcells on the same tetrode at
% different sessions, which may actually be the same cell. Does preliminary
% filtering on those cells to determine if possible matches. Creates a .mat
% file that the user should then later go through and select which pairs
% are the same cell, and pass that to "link_cells.m" to link the cell_ids
% together. 

% 1.) Find cells that meet: 
    %a.) Same tetrode 
    %b.) Different Session 
    %c.) Same depth. (Exclude zero or NULL depth tetrodes). 

    
% 2.) Grab additional information about each of those possible matches (cut
% information & session file), save in a struct

% 3.) Do initial filtering to decrease the amount of intervention required by the user.  

dbpe
import CMBHOME.Database.*


%% Part I: Find possible matches:
sf = sprintf('%s/Bill/Projects/database/possible_matches/possible_matches_%s',dbp,num2str(rat_id));
    
[a] = sprintf('SELECT id from tetrode WHERE rat_id = ''%s'' ',num2str(rat_id));
db.prepareStatement(a);
[a] = db.query();
tid = a.id;

a = sprintf('SELECT tetrode_id,session_id,depth from tetrode_session where tetrode_id in(');

for i = 1:length(tid)
    a = [a num2str(tid(i)) ','];
end
a = a(1:end-1);
a = [a ')'];

    
db.prepareStatement('SELECT tetrode_id,session_id,depth from tetrode_session');
[a] = db.query;
tv = [a.tetrode_id,a.session_id,a.depth];

tetrodes = unique(tv(:,1));
possible = cell([length(tetrodes)*3,1]);
er = 0;
for i = 1:length(tetrodes)
    ds = sprintf('SELECT DISTINCT(depth) FROM tetrode_session WHERE tetrode_id = ''%s'' ',num2str(tetrodes(i)));
    db.prepareStatement(ds);
    [a] = db.query();
    depth = a.depth;
    
    for k = 1:length(depth)
        a = sprintf('SELECT DISTINCT(session_id) from tetrode_session WHERE (tetrode_id = ''%s'' AND depth = ''%s'') ',num2str(tetrodes(i)), num2str(depth(k)));
        db.prepareStatement(a);
        [a] = db.query;
        session_list = a.session_id;
        er = er + 1;
        
        if length(session_list)>1 % if more than 1 session at a depth, then possible matches
           for l = 1:length(session_list)
               a = sprintf('SELECT cell_id FROM cell_session WHERE (session_id = ''%s'')',num2str(session_list(l)));
               db.prepareStatement(a);
               [a] = db.query;
               cells = a.cell_id;
               
               for m = 1:length(cells)
                  a =  sprintf('SELECT tetrode_id from cell WHERE (id = ''%s'')',num2str(cells(m)));
                  db.prepareStatement(a);
                  [a] = db.query();
                  tid2 = a.tetrode_id;
                  if (tid2 == tetrodes(i));possible{er}(l) = cells(m); end
               end
               
           end
        end

    end
    

end

% Get rid of emptys & save
for i = 1:length(possible)
    ev(i) = isempty(possible{i});
end
possible(ev) = [];

possible = struct('cid',possible);

save(sf,'possible');

clearvars -except db possible sf

%% Part II: Get additional information about each cell session of interest
rv = length(possible);
%{
possible.cel = cell(rv,1);
possible.fname = cell(rv,1);
%}

for i = 1:rv

    cv = length(possible(i).cid);

    for k = 1:cv
        
        a = sprintf('SELECT session_id,cut_id FROM cell_session WHERE (cell_id = ''%s'') ',num2str(possible(i).cid(k)));
        db.prepareStatement(a);
        a = db.query();
        sid = a.session_id(1);
        cut = a.cut_id(1);

        a = sprintf('SELECT tetrode_id FROM cell WHERE (id = ''%s'') ',num2str(possible(i).cid(k)));
        db.prepareStatement(a);
        a = db.query();
        tid = a.tetrode_id(1);


        a= sprintf('SELECT trode from tetrode WHERE (id = ''%s'') ',num2str(tid));
        db.prepareStatement(a);
        a = db.query();
        trode = a.trode(1);

        possible(i).cel{k} = [trode cut];

        a = sprintf('SELECT filepath_cmb FROM session WHERE(id = ''%s'') ',num2str(sid));
        db.prepareStatement(a);
        a = db.query();
        possible(i).fname(k) = a.filepath_cmb; 

    end


end

save(sf,'possible');

end





