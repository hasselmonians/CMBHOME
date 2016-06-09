function n = selectLocalTetrode(self,cel)
% Returns the index of the lfp that corresponds to the local LFP of your cell.

% wchapman 20130403

for i = 1:length(self.b_lfp)    
    t{i} = self.b_lfp(i).channel_name{end};
end

if strfind(t{1},'ncs') %Neuralynx
    for i = 1:length(self.b_lfp)
        u = t{i}(4:end);
        u2 = strfind(u,'.ncs')-1;
        tn(i) = str2num(u(1:u2));
    end
elseif strfind(t{1},'eeg') %Axona
    for i = 1:length(self.b_lfp)
        
    end
end

n = find(tn == cel(1));
