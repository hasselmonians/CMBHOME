function eeg_files = dacq_eeg(basePath, prefix)
	%Provides list of eeg files with the defined prefix int eh defined datapath
	%folder.
    
	file_list = listfiles(basePath);
    file_list = alphabetizeEEGfnames(file_list);

	required = [prefix,'.eeg'];
	count=1;
    
	for i = 1:length(file_list)
		iseeg = strfind(file_list{i},required);
		
		if ~isempty(iseeg)
			eeg_files{count} = [basePath filesep file_list{i}];
			count = count+1;
		end
    end
        
end