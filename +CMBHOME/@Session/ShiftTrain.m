function self = ShiftTrain(self, delta)
    % Circularly shifts the spike train by a given amount
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
    
    keyboard
    dur = self.b_ts(end) - self.b_ts(1);
    
    for i = 1:numel(self.spike)
        if ~isempty(self.spike(i).ts)
            t = mod(self.spike(i).ts + delta, self.b_ts(end));
            
        end
    end
    
    for i = 1:n
        delta = range(1) + (range(2)-range(1))*rand(length(self.cel_x{1}),1);
        ts = mod(self.spike(self.cel(1), self.cel(2)).ts + delta, self.b_ts(end));
        Spk(1,i) = CMBHOME.Spike('ts',ts, 'vid_ts', self.b_ts);
    end
    
    self.spike = Spk;
    self = self.AlignSpike2Session;
    
    self.cel = [1 1];
end