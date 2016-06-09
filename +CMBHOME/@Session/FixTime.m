function root = FixTime(root)
% Aligns the session to start at time 0.
%
% root = root.FixTime;
    
    for i = 1:numel(root.spike)
        root.spike(i).ts = root.spike(i).ts - root.b_ts(1);
    end
    
    for i = 1:numel(root.b_lfp)
        root.b_lfp(i).ts = root.b_lfp(i).ts - root.b_ts(1);
    end
    
    root.b_ts = root.b_ts - root.b_ts(1);
    
    root = root.AlignSpike2Session;
    root = root.AlignSpike2LFP;
    
    root.epoch = [-inf inf];

end