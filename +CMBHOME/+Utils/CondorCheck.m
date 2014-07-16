function CondorCheck(varargin)
% CMBHOME.Utils.CondorCheck()

% Bill 2012.12.11

%% Parse Inputs
p = inputParser;

p.addParamValue('user',      'lost',     @(x) ischar(x));
p.addParamValue('an',        'lost',     @(x) ischar(x));
p.addParamValue('endline',   'Done...',  @(x) ischar(x));

origPath = pwd; % So that we can come back at the end of the function.

cd([dropboxPath '/condor/' user '/' an])

%% find all relevant files
outFiles = dir('*.out');
logFiles = dir('*.log');
errFiles = dir('*.err');


% parse files
for i = 1:length(outFiles)

  % look for finalOutStr in *.out file
  fid = fopen(fullfile(condorPath,outFiles(i).name)); % out file  
  comparisons = false;
  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    comparisons = comparisons || strcmp(tline,finalOutStr);
  end
  fclose(fid);

  if comparisons ~=0
    delete(outFiles(i).name);
    delete(logFiles(i).name);
    delete(errFiles(i).name);
  end
  
end

cd(origPath)
end
