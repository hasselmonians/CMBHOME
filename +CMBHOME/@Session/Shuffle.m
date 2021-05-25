function self = Shuffle(self, varargin)
    % Shuffles spike times ALL cells by random amount.
    % 
    % input:
    %   self: The original object, with self.cel set the cell desired to be
    %         shuffled.
    %
    %
    % parameters:
    %   n: Number of shuffles to perform (1)
    %      Note: If >1 then returns self.spike as larger
    %
    %   range: How large of a shuffle to perform(seconds). 
    %          ([0.003, Duration-0.003])
    %
    %   
    % output:
    %   self: The modified object with shuffled spike times.
    %
    %  
    % See methods in Hinman et al, Neuron, 2016.
    % NOTE: Time must start at 0
    % root2 = root.Shuffle();

    if self.b_ts(1) ~= 0
        error('Run root.FixTime first');
    end
    dur = self.b_ts(end) - self.b_ts(1);
    p = inputParser;
    p.addParamValue('range', [0.003 self.b_ts(end)-0.003], @(x) numel(x)==2)

    p.parse(varargin{:})
    range = p.Results.range;
    
    for i = 1:size(self.cells,1)
        ii = self.cells(i,1);
        jj = self.cells(i,2);
        delta = range(1) + (range(2)-range(1))*rand(length(self.spike(ii,jj).ts),1);
        t = mod(self.spike(ii,jj).ts + delta, dur);
        spk(ii,jj) = CMBHOME.Spike('ts',t, 'vid_ts', self.b_ts);
    end
    
    self.spike = spk;
    self = self.AlignSpike2Session;
    
    old = self.cel;
    self.cel = [1 1];
    self.cel = old;
end