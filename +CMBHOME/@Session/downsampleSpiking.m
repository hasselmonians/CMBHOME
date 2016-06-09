function self = downsampleSpiking(self,cel,f)
% Randomly removes spikes to achieve a lower mean firing rate
%
% root = root.downsampleSpiking(cel, newMean)

if ~exist('cel','var')
    if isempty(cel)
        cel = self.cel;
    end
else
    self.cel = cel;
end

  % duration of current epoch selection
  t = sum(diff(self.epoch,[],2));  
  
  % number of spikes needed given desired spike rate
  nSpk = round(t * f); 
  
  % only downsample if the new number of required spikes is lower than
  % already present
  spk_i = CMBHOME.Utils.ContinuizeEpochs(self.cel_i);
  spk_ts = CMBHOME.Utils.ContinuizeEpochs(self.cel_ts);

  if nSpk > length(spk_i)
    return;
  else
    % randomly select among existing spikes
    randOrder = randperm(length(spk_i));
    self.spike(cel(1),cel(2)).i = spk_i(randOrder(1:nSpk));
    self.spike(cel(1),cel(2)).ts = spk_ts(randOrder(1:nSpk));
    self.cel = [];
    self.cel = cel;
    self = self.AlignSpike2LFP;
    self = self.AlignSpike2Session;
  end
end

