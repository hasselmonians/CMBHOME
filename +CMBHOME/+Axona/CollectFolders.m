function folders = CollectFolders(base_path)
% returns full path to directories within base_path

    in_base_path = dir(base_path);
    
    folders = [];
    counter = 1;
    for folderi = 1:length(in_base_path)
        if in_base_path(folderi).isdir && ~strcmp(in_base_path(folderi).name,'.') && ~strcmp(in_base_path(folderi).name,'..')
            folders{counter} = fullfile(base_path, in_base_path(folderi).name);
            counter = counter+1;
        end
    end            
            
end