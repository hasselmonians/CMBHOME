function CondorSubmit(varargin)
% CMBHOME.Utils.CondorSubmit('f2run','myFunction','logs',1,'RAMreq',1000,'jvm',0,'user','wchapman','an','myFunction_test')
% Creates & submits a condor job. See parser for input descriptions

% Bill 2012.10.16

%% Parse Inputs
p = inputParser;

p.addParamValue('f2run',        'none',     @(x) ischar(x)); % The name of the function that you would like to run
p.addParamValue('logs',         1,          @(x) isnumeric(x)); % Whether or not to keep log files
p.addParamValue('RAMreq',       1000,       @(x) isnumeric(x)); % The amount of ram (MB)
p.addParamValue('holderVar',    [],         @(x) iscell(x)|isstruct(x)); % Variable array with input information for each iteration
p.addParamValue('jvm',          0,          @(x) isnumeric(x)); % Whether or not to use the java machine
p.addParamValue('useWindows',   0,          @(x) isnumeric(x)); % Whether of not to use Windows machines
p.addParamValue('user',         'lost',     @(x) ischar(x)); % The user where to keep the output files (relative to dropboxPath/condor/ )
p.addParamValue('an',           'lost',     @(x) ischar(x)); % Subfolder to keep the log files and output files in. (Stands for "analysis name")
p.addParamValue('superNice',    1,          @(x) isnumeric(x)); % By default, only submits to machines owner by you. Useful for prototyping.

p.parse(varargin{:});

f2run = p.Results.f2run;
an = p.Results.an;
holder_var = p.Results.holderVar;
jvm = p.Results.jvm;
logs = p.Results.logs;
RAMreq = p.Results.RAMreq;
superNice = p.Results.superNice;
user = p.Results.user;
useWindows = p.Results.useWindows;

nowStr = datestr(now,30); % Time stamp used to keep track of intermediate output files.

%% double check script is prepared to exit after call
% If an error is thrown here, please add "exit" at the line before "end"

fprintf('checking for exit command: ');
validScript = unix(['grep ''exit'' ' which(f2run)]);
if ~validScript
  %error('f2run must call ''exit'' at the end to avoid hanging the job.')
end


%% Create EXE
% If user wants to use windows computers, then the first thing to do is
% compile the exe to use on those computers

if useWindows == 1
    if ~ispc
        error('Must submit from a windows machine to execute on windows nodes');
    else 
        addpath(pwd);
        compileState = ['mcc -o ' f2run]; % mcc = "Compile" -o f2run defines which function we are compiling.
        compileState = [compileState ' WinMain:end -T link:exe -d' pwd ]; %Make the exe with name f2run.exe in current directory.
        compileState = [compileState ' -w enable:specified_file_mismatch -w enable:repeated_file  -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license']; %Bunch of options to help compile smoothly
        compileState = [compileState ' -v' pwd '\' f2run '.m']; % Point to the main file used
        compileState = [compileState ' -a' dropboxPath '\Bill\CMBHOME_package']; %Add the minimized CMBHOME package to the EXE.
        
        eval(compileState); %run it
        movefile([f2run '.exe'],[f2run '.WINDOWS']); % Final product: f2run.WINDOWS
    end 
end
 
%% Create the shell for mac & linux machines
fid = fopen([f2run '.OSX'],'w');
fprintf(fid,'#!/bin/bash \n # \n \n');
fprintf(fid,'source /condor/pathv \n');

if jvm==0
    fprintf(fid,'exec matlab -nojvm -nodisplay -nosplash -r "$*" \n'); % requires that /Applications/MATLAB_VERSION/bin be in your path
else
     fprintf(fid,'exec matlab -nodisplay -r "$*" \n');
end
fclose(fid);

%% create holder_var file
holder_var_name = [pwd '/' 'holder_var_',nowStr,'.mat'];
save(holder_var_name,'holder_var');

%% Create Submit file
fid = fopen([f2run,'.submit'],'w');
fprintf(fid,'universe = vanilla \n');
fprintf(fid,'executable = %s.$$(OpSys) \n',f2run);
argument = ['arguments = ' f2run '($(PROCESS),''' holder_var_name ''')'];
fprintf(fid,'%s \n',argument);
fprintf(fid,'transfer_input_files = %s \n',which(f2run));

req = 'requirements = ((OpSys == "OSX")'; % Allow macs
req = [req '|| (Target.OpSys == "LINUX")'];

if useWindows ==1
    req = [req '||(OpSys == "WINDOWS"))'];
else
    req = [req ')'];
end

if superNice == 1
    %req = [req ')'];
    req = [req '&& (TARGET.Memory >=' num2str(RAMreq) ,') '];
    req = [req '&& ((TARGET.OWNER == "' user '")) \n'];
else
    req = [req ' \n'];
end

fprintf(fid,req);

if logs == 1
% Sets the paths to keep the logs in, based on inputs

    mkdir([dropboxPath '/condor/' user '/' an])

    sv = ['output = ' dropboxPath '/condor/' user '/' an '/$(CLUSTER).$(PROCESS).out \n'];
    fprintf(fid,sv);
    
    sv = ['log = ' dropboxPath '/condor/' user '/' an '/$(CLUSTER).$(PROCESS).log \n'];
    fprintf(fid,sv);
    
    sv = ['error = ' dropboxPath '/condor/' user '/' an '/$(CLUSTER).$(PROCESS).err \n'];
    fprintf(fid,sv);
end

lengthv = length(holder_var);
fprintf(fid,'queue %d',lengthv);

fclose(fid);

%% Submit to manager:
if ispc
    dos(['condor_submit' f2run '.submit']); 
else
    if isempty(getenv('CONDOR_CONFIG'))
        setenv('CONDOR_CONFIG','/condor/etc/condor_config');
        setenv('PATH',['/condor/bin:/condor/sbin:',getenv('PATH')]);
    end
    unix(['chmod 777 ' f2run '.OSX']);
    unix(['cp ',f2run,'.OSX ',f2run,'.LINUX']);
    unix(['condor_submit ' f2run '.submit']);
end
