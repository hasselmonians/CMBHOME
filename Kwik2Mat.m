function Kwik2Mat(fname, vid_ts)

    hinfo = hdf5info(fname);
    for i = 1:length(hinfo.GroupHierarchy.Groups)
        gn{i} = hinfo.GroupHierarchy.Groups(i).Name;
    end

    ind = find(strcmp(gn,'/channel_groups'));
    nGroups =  length(hinfo.GroupHierarchy.Groups(ind).Groups);

    sr = h5readatt(fname, ['/recordings/0'], 'sample_rate');
    
    for i = 1:nGroups
        try
            ts = (1/sr)*double(h5read(fname, ['/channel_groups/' num2str(i-1) '/spikes/time_samples'])) + ...
                (1/sr)*double(h5read(fname, ['/channel_groups/' num2str(i-1) '/spikes/time_fractional']))/255;

            clust = double(h5read(fname, ['/channel_groups/' num2str(i-1) '/spikes/clusters/main']));
            clust_auto = double(h5read(fname, ['/channel_groups/' num2str(i-1) '/spikes/clusters/original']));

            cels = unique(clust) +1;
            for k = 1:length(cels)
               inds = clust+1 == k  ;
               Spike(i,cels(k)) = CMBHOME.Spike('ts', ts(inds), 'vid_ts', vid_ts);
            end
            
            cels = unique(clust_auto) + 1;
            for k = 1:length(cels)
                inds = clust_auto+1 == k;
                Spike_Auto(i,cels(k)) = CMBHOME.Spike('ts', ts(inds), 'vid_ts', vid_ts);
            end
            
        catch ME
            if strcmp(ME.identifier, 'MATLAB:imagesci:h5read:libraryError')
                warning(['Channel Group ' num2str(i-1) ' is not clustered']);
            end
        end
    end

    keyboard
    
end