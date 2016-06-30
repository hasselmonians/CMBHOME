function sessions = ImportGUI()

    %
    % Will displace the location of the saved CMBHOME Objects, which can also be loaded easily as:
    % load(sessions(1).mat).

    %addpath C:\Users\Holger\Dropbox' (hasselmonians)'\hdannenb\InvivoImport\
    %{
    sessions = invivoImport('basePath', 'C:\Users\Holger\Dropbox (hasselmonians)\hdannenb\UnitRecordingData\Testmouse_1',...
    'cutFiles', {'160304_1_2.cut'}, ...
    'setFiles', {'160304_1.set','160304_2.set','160307_1.set','160308_1.set'});
    %}         

    basePath = uigetdir(pwd,'Select Base Dir (Contains .set)');
    od=pwd;
    cd(basePath);
    cutFiles = uigetfile({'*.cut;*.clu.*'}, 'Select intended cut files (1 per tetrode)','MultiSelect','on');
    if ~iscell(cutFiles), cutFiles={cutFiles};end
    setFiles = uigetfile('*.set', 'Select intended set files (1 per session)','MultiSelect','on');
    if ~iscell(setFiles), setFiles={setFiles};end
    cd(od);

    disp('--------------------------------------------------------------------')
    disp('|                         WORKING                                  |')
    disp('--------------------------------------------------------------------')

    sessions = invivoImport('basePath',basePath,...
                            'cutFiles',cutFiles,...
                            'setFiles',setFiles);

    for i = 1:length(sessions)
        disp(sessions(i).mat)
    end


end