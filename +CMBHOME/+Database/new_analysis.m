function new_analysis(db,analysis_file,cells,sessions,epochLabel,path)
% Makes use of both the database and Condor for an extensive analysis to be
% performed on a subset of cells. Your analysis will become a BLOB entry in the
% "Analysis" table of database.

% Input Arguments:
% db: The database object
% analysis_file: Link to the .m file that should be run on each cell.
% cells: A list of cellIDs which you want to run the analysis on. ('All'
%        for all cells).
% epochLabel: The label of the epoch (defaults to entire session). 
% path: The path relative to your Dropbox folder where we should store
%       temporary files used in the analysis.


state = ['CREATE TABLE IF NOT EXISTS `hasselmo`.`analysis1`',...
        '(`cell_id` MEDIUMINT(9) NOT NULL ,',...
        '`session_id` MEDIUMINT(9) NOT NULL ,',...
        '`epoch_label` VARCHAR(20) NOT NULL ,'...
        'PRIMARY KEY (`cell_id`, `session_id`, `epoch_label`) ,'...
        'CONSTRAINT `fk_analysis1_cell_epoch1`',...
        'FOREIGN KEY (`cell_id` , `session_id` , `epoch_label` )',...
        'REFERENCES `hasselmo`.`cell_epoch` (`cell_id` , `session_id` , `epoch_label` )',...
        'ON DELETE NO ACTION',...
        'ON UPDATE NO ACTION)',...
        'ENGINE = InnoDB'];



%{
dbpe
import CMBHOME.Database.*

lv = length(cells);

fid = fopen(fullfile(path,'To_workon.txt'),'w');



%% Fill out the list of what to work on
for i = 1:lv
    db.prepareStatement('SELECT filepath_cmb FROM session WHERE (id = ''%s'')',sessions(i));
    [a] = db.query();
    fpCMB = a.filepath_cmb;
    
    sv = sprintf('SELECT * FROM cell_session WHERE (cell_id =''%s'' AND session_id = ''%s'')',cells(i),sessions(i));
    db.prepareStatement(sv);
    [a] = db.query();
    cutID = a.cut_id;
    
    sv = sprintf('SELECT tetrode_id FROM cell WHERE (cell_id = ''%s'')',cells(i));
    db.prepareStatement(sv);
    [a] = db.query();
    tid = a.tetrode_id;
    
    sv = sprintf('SELECT trode FROM tetrode WHERE (id = ''%s'')',num2str(tid));
    db.prepareStatement(sv);
    [a] = db.query();
    trode = a.trode;
    
    fprintf(fid,'%s [%s %s]\n',fpCMB,num2str(trode),num2str(cutID));
end

fclose(fid);

%% Create the submit file:
fid = fopen(fullfile(path,'mlbsub.submit'),'w');    % Create the submit file, 'w' means create the file from scratch
fprintf(fid,'universe = vanilla\n'); % use the vanilla universe
fprintf(fid,'executable =  ./condor/matlab.$$(OpSys) \n'); % Use the matlab wrapper that we made
fprintf(fid,'arguments = myfunc2($(PROCESS)) \n'); % Thse are the commands passed to the matlab wrapper 
fprintf(fid,'requirements = (TARGET.OpSys == "WINDOWS" ) && (TARGET.Arch == "INTEL") && (TARGET.Memory >1000) \n'); % Use pretty much all of the computers
fprintf(fid,'transfer_input_files = %s \n',analysis_file); % In case the analysis file wasn't on the computer, transfer it.
fprintf(fid,'queue %d\n',lv); % Repeat the process for "lv" loops
fclose(fid);

%% Later
% Need to run condor_submit
% After the processes are done running, run a collector to put the results
% into a single file and submit that to the database. 
%}
