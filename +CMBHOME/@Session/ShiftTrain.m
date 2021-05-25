function self = ShiftTrain(self, delta)
    % Circularly shifts ALL spike trains by a given amount
    % 
    % input:
    %   self: The original object, with self.cel set the cell desired to be
    %         shuffled.
    %
    %
    % parameters:
    %   delta: Amount to shift the spike trains by. Can be positive or
    %   negative
    %   
    % output:
    %   self: The modified object with shuffled spike times.
    
    dur = self.b_ts(end) - self.b_ts(1);
    
    for i = 1:size(self.cells,1)
        t = mod(self.spike(self.cells(i,1)).ts + delta, self.b_ts(self.cells(i,2)), dur);
        spk(self.cells(i,1), self.cells(i,2)) = CMBHOME.Spike('ts',t, 'vid_ts', self.b_ts);
    end
       
    self.spike = spk;
    self = self.AlignSpike2Session;
    
    old = self.cel;
    self.cel = [1 1];
    self.cel = old;
    
end