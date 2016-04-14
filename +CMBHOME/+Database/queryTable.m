function SL = queryTable(db)

if nargin == 0
   try
       db = evalin('base','db');
       fprintf('Didnt pass in db variable, so taking from workspace. Make sure this is correct \n');
   catch
       fprintf('Error: Didnt pass in db variable, and no db in workspace. Exiting \n');
       return
   end
end


%% GUI
%import CMBHOME.*
%import CMBHOME.Database.*

t = which('CMBHOME.Database.queryTable');

if ispc
   slash = strfind(t,'\');
else
   slash = strfind(t,'/'); 
end

t = t(1:slash(end-2));
t = [t 'support' filesep];

od = pwd;
cd(t);

t = guidata(queryTable_gui);
cd(od);
%% Extract Information from GUI
for i =1;
    % Group
    grps = get(t.grp,'String');

    % Rat
    rat_name = get(t.rat_name,'String');
    rat_implantDateMin = (get(t.rat_implant1,'String'));
    rat_implantDateMax = (get(t.rat_implant2,'String'));
    rat_implantDepthMin = (get(t.rat_depth1,'String'));
    rat_implantDepthMax = (get(t.rat_depth2,'String'));

    % Tetrode
    tet_numMin = (get(t.tetrode_num1,'String'));
    tet_numMax = (get(t.tetrode_num2,'String'));
    tet_notes= get(t.tetrode_notes,'String');

    % Cell
    cell_numMin = (get(t.cell_num1,'String'));
    cell_numMax = (get(t.cell_num2,'String'));
    cell_morphology = get(t.cell_morphology,'String');
    cell_type = get(t.cell_type,'String');

    % Session
    session_filepath = get(t.session_filepath,'String');
    session_dateMin = (get(t.session_date1,'String'));
    session_dateMax = (get(t.session_date2,'String'));
    session_experimenter = get(t.session_experimenter,'String');
    session_manipulation = get(t.session_manipulation,'String');
    session_task = get(t.session_task,'String');
    session_taskNotes = get(t.session_taskNotes,'String');
    session_enclosure = get(t.session_enclosure,'String');
    session_notes= get(t.session_notes,'String');

    % Tetrode Session
    ts_depthMin = (get(t.ts_depth1,'String'));
    ts_depthMax = (get(t.ts_depth2,'String'));
    ts_anatomy = get(t.ts_anatomy,'String');
    ts_notes = get(t.ts_notes,'String');


    % Epoch
    epoch_Label = get(t.epoch_label,'String');
    epoch_notes = get(t.epoch_notes,'String');


    % Cell Epoch
    ce_limits = get(t.ce_table,'Data');
    ce_columnNames = get(t.ce_table,'RowName');

    
end


%% ------We Now begin parsing into queries for each table------%
%% Rat
ratState = 'SELECT id,name FROM rat WHERE (grp_name IN (';

for i = 1:length(grps)
    ratState = [ratState ,'''', grps{i},'''', ','];
end
ratState = ratState(1:end-1);
ratState = [ratState ')'];

if notEmpty(rat_name)
    ratState = [ratState 'AND name LIKE ''%' rat_name '%'')'] ;
end

if str2num(rat_implantDateMin) > 000101
    ratState= [ratState 'AND (implant_date >=' rat_implantDateMin ')'];
end

if str2num(rat_implantDateMax) < 991231
    ratState = [ratState 'AND (implant_date <=' rat_implantDateMax ')'];
end

if str2num(rat_implantDepthMin) > -1
    ratState = [ratState 'AND final_depth >=' rat_implantDepthMin ')'];
end

if str2num(rat_implantDepthMax) < 9999
   ratState = [ratState 'AND final_depth <=' rat_implantDepthMax ')'];
end

ratState = [ratState ')'];


%% Tetrode
tetState = 'SELECT id,trode FROM tetrode WHERE (';

if str2num(tet_numMin) > -1
    tetState = [tetState '(trode >=' tet_numMin ')'];
end

if str2num(tet_numMax) < 17
    tetState = [tetState 'AND (trode<=' tet_numMax ')'];
end

if ~strcmp(tet_notes,'%')
   tetState = [tetState 'AND (notes LIKE ''%' tet_notes '%'')'] ;
end

tetState = [tetState ')'];

if length(tetState) <= 37
   tetState = 'SELECT id,trode FROM tetrode';
end


%% Cell
cellState = 'SELECT id,morphology,tetrode_id FROM cell WHERE (';

if str2num(cell_numMin) > -1
   cellState = [cellState '(cut_id >=' cell_numMin ')'];
end

if str2num(cell_numMax) < 99999
   cellState = [cellState 'AND (cut_id <=' cell_numMax ')'];
end

if ~strcmp(cell_morphology,'%')
    cellState = [cellState 'AND (morphology LIKE ''%' cell_morphology '%)'];
end

cellState = [cellState ')'];

if length(cellState) <=50
    cellState = 'SELECT id,morphology,tetrode_id FROM cell';
end


%% Session
sessionState = 'SELECT filepath_cmb,id,rat_id FROM session WHERE (';
im = 0;

if str2num(session_dateMin) > 000101
   sessionState = [sessionState '(date >=' session_dateMin ')'];
end

if str2num(session_dateMax) < 991231
   if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
   sessionState = [sessionState '(date <=' session_dateMax ')'];
end

if notEmpty(session_filepath)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
   sessionState = [sessionState '(filepath_cmb LIKE ''%' session_filepath '%'')'];
end

if notEmpty(session_experimenter)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
   sessionState = [sessionState '(experimenter LIKE ''%' session_experimenter '%'')'];
end

if notEmpty(session_manipulation)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
   sessionState = [sessionState '(manipulation LIKE ''%' session_manipulation '%'')'];
end

if notEmpty(session_task)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
    sessionState = [sessionState '(task LIKE ''%' session_task '%'')'];
end

if notEmpty(session_taskNotes)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
    sessionState = [sessionState '(task_notes LIKE ''%' session_taskNotes '%'')'];
end

if notEmpty(session_enclosure)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
    sessionState = [sessionState '(enclosure LIKE ''%' session_enclosure '%'')'];
end

if notEmpty(session_notes)
    if im~=0;sessionState = [sessionState 'AND ']; else im =1;end
    sessionState = [sessionState '(notes LIKE ''%' session_notes '%'')'];
end

sessionState = [sessionState ')'];

if length(sessionState) <= 51
    sessionState = 'SELECT filepath_cmb,id,rat_id FROM session';
end

%% TetrodeSession
tsState = 'SELECT tetrode_id,session_id,depth,anatomy FROM tetrode_session WHERE (';

if str2num(ts_depthMin) > -1
   tsState = [tsState '(depth >=' ts_depthMin ')']; 
end

if str2num(ts_depthMax) < 9999
    tsState = [tsState '(AND depth <=' ts_depthMax ')'];
end

if notEmpty(ts_anatomy)
   tsState = [tsState 'AND (anatomy LIKE ''%''' ts_anatomy '%'')'];
end

if notEmpty(ts_notes)
   tsState = [tsState 'AND (notes LIKE ''%''' ts_notes '%'')'];
end

tsState = [tsState ')'];

if length(tsState) <= 72
    tsState = 'SELECT tetrode_id,session_id,depth,anatomy FROM tetrode_session';
end

%% Epoch 
epochState = 'SELECT session_id,label,epoch FROM epochs WHERE (';
epochState = [epochState 'label LIKE ''%' epoch_Label '%'')'];

if notEmpty(epoch_notes)
   epochState = [epochState 'AND notes LIKE ''%''' epoch_notes '%'')'];
end

%% Cell Epoch
state = 'SELECT cell_id,session_id,epoch_label,';
for i = 1:size(ce_limits,1)
    state = [state ce_columnNames{i} ',']; 
end

if get(t.ce_rm,'Value')
   state = [state 'ratemap,']; 
end

if get(t.ce_acorr,'Value')
   state = [state 'acorr,']; 
end
state = state(1:end-1);
state = [state ' FROM cell_epoch WHERE('];
state= [state '(epoch_label LIKE ''%' epoch_Label '%'') AND'];

stateb = state;

cv = 0; 
for i = 1:size(ce_limits,1)
    if ~isempty(ce_limits{i,1})
       cv = cv+1;
       state = [state ' (' ce_columnNames{i} '>=' (ce_limits{i,1}) ') AND'];
    end
    
    if ~isempty(ce_limits{i,2})
        cv = cv+1;
        state = [state ' (' ce_columnNames{i} '<=' (ce_limits{i,2}) ') AND'];
    end
end

state = [state(1:end-3) ')'];


% close it
close(t.figure1);

%% Query the database and filter results

db.prepareStatement(state);
cell_epoch = db.query(); 

db.prepareStatement(ratState);
rat = db.query();

db.prepareStatement(sessionState)
session = db.query();

db.prepareStatement(tetState)
tetrode = db.query();

db.prepareStatement(cellState)
cell = db.query();

db.prepareStatement(tsState)
tetrodeSession = db.query();

db.prepareStatement('SELECT cell_id,session_id,cut_id FROM cell_session');
cs = db.query();

% Filter session based on rat:
inds = ismember(session.rat_id,rat.id);
session = allFields(session,inds);

%filter tetrodeSession on session
inds = ismember(tetrodeSession.session_id,session.id);
tetrodeSession = allFields(tetrodeSession,inds);

% Filter cell on tetrodeSession
inds = ismember(cell.tetrode_id,tetrodeSession.tetrode_id);
cell = allFields(cell,inds);

% Filter cell_epoch on session
inds = ismember(cell_epoch.session_id,session.id);
cell_epoch = allFields(cell_epoch,inds);

% Filter cell_epoch on cells
inds = ismember(cell_epoch.cell_id,cell.id);
cell_epoch = allFields(cell_epoch,inds);

%% Create SL
fn = fieldnames(cell_epoch);
fn = fn(3:end);

for i = 1:length(cell_epoch.cell_id)
    SL(i).id.session = cell_epoch.session_id(i);
    SL(i).id.cell = cell_epoch.cell_id(i);
    SL(i).id.tetrode = cell.tetrode_id((cell.id == cell_epoch.cell_id(i)));
    SL(i).id.rat = session.rat_id(session.id == cell_epoch.session_id(i)) ;
    SL(i).id.epochLabel = cell_epoch.epoch_label{i};
    
    for k = 1:length(fn)
        eval(['SL(i).' fn{k} '= cell_epoch.' fn{k} '(i);']);
    end
    
    cm = find(cs.cell_id == SL(i).id.cell);
    sm = find(cs.session_id == SL(i).id.session);    
    c = intersect(cm,sm);
    
    SL(i).cel = [tetrode.trode(find(tetrode.id == SL(i).id.tetrode,1)) cs.cut_id(c)];
    SL(i).filepath_cmb = session.filepath_cmb{find(session.id == SL(i).id.session,1)};
end

SL = SL(:);

end


function fv = notEmpty(sv)

fv = strcmp(sv,'%');
fv = ~fv;

end

function sv = allFields(sv,inds)
    
fn = fieldnames(sv);

for i = 1:length(fn)
    eval(['sv.' fn{i} '= sv.' fn{i} '(inds);']);
end

end