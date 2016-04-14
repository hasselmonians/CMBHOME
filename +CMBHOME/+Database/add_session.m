function session_id = add_session(db,subObj)
%
% Adds a new session (single recording period) to the central database
% Bill 2012.08.30

import CMBHOME.Database.*

%% rat_id
if subObj.newRat == 1
    rat_id = CMBHOME.Database.add_rat(db,subObj);
else
    state = sprintf('SELECT id FROM rat WHERE (name = ''%s'' AND grp_name = ''%s'')',subObj.ratname,subObj.grp_name);
    db.prepareStatement(state);
    [a] = db.query();

    if isempty(a.id)
        error('Could not locate the rat name that you specified, exiting');
    else
        rat_id = [num2str(a.id(1))];
    end
end

%% Load the CMBobject to get additional fields
c = load(subObj.filepath_cmb);
root = c.root;
clear c; 

% Get time and date from the CMBobject as well as time updated (set to now)
if ~isempty(root.date_created)
    datev = [datestr(root.date_created,'yymmdd')];
    timev = [datestr(root.date_created,'HHMM')];
else
    datev = 'NULL';
    timev = 'NULL';
end

dateUpdated = [datestr(now,31)];

% Duration
duration = [num2str(root.epoch(end) - root.epoch(1))];

% Filepath_blfp
if ~isempty(root.path_lfp)
    filepath_blfp = root.path_lfp{1};    
else
    filepath_blfp = 'NULL';    
end

% Filepath_raw
if ~isempty(root.path_raw_data)
    filepath_raw = root.path_raw_data{1}; 
else
    filepath_raw = 'NULL';
end

%% Copy over to RatDatabase & add all of the LFP information.
%{
fpN = [dropboxPath '/' subObj.grp_name '/' subObj.ratname '/' root.name{end}];
root.epoch = [-inf inf];
root.cel = [];
root.active_lfp = [];

mkdir([dropboxPath '/' subObj.grp_name '/' subObj.ratname '/']);

% Load all LFP data
root = root.LoadLFP(1:length(root.path_lfp));
root.active_lfp=[];
root.cel = [];

% Align spike,seesion,lfp all together
root = root.AlignSpike2LFP;
root = root.AlignSpike2Session;
root = root.AlignSpike2LFP;

save(fpN,'root'); % This is the version that will be in dropboxPath/RatDatabase
%}
fpN = subObj.filepath_cmb;
%% Enter the session:
if iscell(filepath_blfp);filepath_blfp = CMBHOME.Utils.array2path(filepath_blfp);end

enter_session = sprintf('INSERT INTO session (rat_id,date,time,dateUpdated,filepath_cmb,filepath_blfp,filepath_raw,experimenters,enclosure,task,manipulation,task_notes,ephys_notes,duration) VALUES(''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'')', ...
    rat_id, ...
    datev, ...
    timev, ...
    dateUpdated, ...
    fpN, ...
    filepath_blfp, ...
    subObj.filepath_cmb, ...
    subObj.experimenters, ...
    subObj.enclosure, ...
    subObj.task, ...
    subObj.manipulation, ...
    subObj.task_notes, ...
    subObj.ephys_notes, ...
    duration);

enter_session = CMBHOME.Database.fns(enter_session);

% Give the command to mySQL
db.prepareStatement(enter_session);
[a] = db.query();

% Now return the session_id just entered for the user's own knowledge.
db.prepareStatement('SELECT LAST_INSERT_ID()');
[a] = db.query();
session_id = a.LAST_INSERT_ID__;

%% Automatically add the full session epoch
sv = sprintf('INSERT INTO epochs (label,session_id,epoch,notes) VALUES(''full_session'',''%s'',''[0 inf]'',NULL)',num2str(session_id));
db.prepareStatement(sv);
db.query();

%% Add all of the tetrodes
dv = NaN(length(subObj.trode),2);

if subObj.newRat == 1
   for i = 1:length(subObj.trode)      
       dv(i,1) = i;
       dv(i,2) = str2num(CMBHOME.Database.add_tetrode(db,rat_id,i));
   end
end

%% Add all of the cells
[nc,~] = size(root.cells);
for i = 1:nc
   trode = root.cells(i,1);
   tetrode_id = dv((dv(:,1)==trode),2);
   cut_id = root.cells(i,2);
   
   cell_id(i) = CMBHOME.Database.add_cell(db,tetrode_id,cut_id,subObj.cells(i),session_id);
end

end
