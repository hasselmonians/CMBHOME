function root = AlignSpike2LFP(root)
% root = root.AppendThetatoSpike;
%
% For every root.b_lfp, aligns root.b_lfp.ts to spike times for all cells
% in root.spike.ts.
%
% andrew 28 sept 2010

import CMBHOME.Utils.*

root.cell_thresh = 0;

if isempty(root.b_lfp)
    
    disp('No LFP data was loaded, adding all root.path_lfp');
    
    root = root.LoadLFP(1:length(root.path_lfp));
    
    if isempty(root.b_lfp), error('No LFP files were loaded'); end
    
end

for i = 1:length(root.b_lfp)
        
    if isempty(root.b_lfp(i).signal), continue, end
    
    for j = 1:size(root.cells,1)
        
        spk_ts = root.spike(root.cells(j,1), root.cells(j,2)).ts;
        
        root.spike(root.cells(j,1),root.cells(j,2)).i_lfp{i} = SpkInds(root.b_lfp(i).ts, spk_ts);
        
    end
    
end