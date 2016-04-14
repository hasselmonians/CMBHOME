function populate_db(db)
% Goes through and finds which epochs which don't yet have cell_epoch
% defined. Takes those sessions and then creates a condor submit file which
% will look for computers to populate those statistics.

dbpe
import CMBHOME.Database.*

%% Select epochs without cell_epoch
db.prepareStatement('SELECT * FROM epochs');
[a] = db.query();

labels = a.label;
sids = a.session_id;
epochs = a.epoch;
nr = [];
for i = 1:length(labels)
    state = sprintf('SELECT cell_id,session_id,epoch_label from cell_epoch WHERE(session_id =''%s'' AND epoch_label = ''%s'')', ...
        num2str(sids(i)),...
        labels{i});
    
    db.prepareStatement(state);
    [a] = db.query();
    
    if isempty(a.cell_id);nr = [nr,i];end % Need run is now an index of all of rows that need to be populated
    
end


%% Prepare the Submit files 
fid = fopen('to_workon.txt','w');

% The file to read from, each row is information needed for a seperate loop
for i = 1:length(nr)
    rv = nr(i);
    sid = num2str(sids(rv));
    label = labels{i};
    epoch = epochs{i};
    
    state = sprintf('SELECT filepath_cmb FROM session WHERE (id = ''%s'')',sid);
    db.prepareStatement(state);
    [a] = db.query();
    fname = a.filepath_cmb;
    
    fprintf(fid,'%s,%s,%s,%s\n',fname{1},label,epoch,sid);
end

fclose(fid);

%% Write the submit file
fid = fopen('mlbsub.submit','w');
fprintf(fid,'universe = vanilla\n');
fprintf(fid,'executable =  ./condor/matlab_jvm.$$(OpSys) \n');
fprintf(fid,'arguments = populate_cell_epoch($(PROCESS)) \n');
%fprintf(fid,'requirements = (TARGET.OpSys == "OSX" ) && (TARGET.Arch == "INTEL") && (TARGET.Memory >100) \n'); 
fprintf(fid,'queue %d\n',length(nr));
fclose(fid);

sprintf('text files ready. Please proceed to mlbsub phase \n')
