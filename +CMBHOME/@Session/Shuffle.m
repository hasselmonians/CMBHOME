function self = Shuffle(self, varargin)
    % Shuffles spike times of current cell by random amount.
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

    p = inputParser;
    p.addParamValue('n', 1,  @(x) numel(x)==1);
    p.addParamValue('range', [0.003 self.b_ts(end)-0.003], @(x) numel(x)==2)

    p.parse(varargin{:})
    n = p.Results.n;
    range = p.Results.range;
    
    for i = 1:n
        delta = range(1) + (range(2)-range(1))*rand(length(self.cel_x{1}),1);
        ts = mod(self.spike(self.cel(1), self.cel(2)).ts + delta, self.b_ts(end));
        Spk(1,i) = CMBHOME.Spike('ts',ts, 'vid_ts', self.b_ts);
    end
    
    self.spike = Spk;
    self = self.AlignSpike2Session;
    
end