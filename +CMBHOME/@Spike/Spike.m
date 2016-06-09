% This is a standard data format spike times class for the CMB

classdef Spike
    
    properties
        i           % index in root.b_ts aligned from more accurate spike times 'ts' below. not directly used by end user
        i_lfp       % cell array of indices in root.b_lfp(i).ts aligned from higher sampled spike times below. not directly used by end user
        ts          % spike timestamps
        tet_name    % tetrode name or index
        prop        % structure with fields like watsonu2, gridness etc. firing properties
    end 
    
    methods (Static)
                  
        [cor, lag] = CrossCorr(ts1, varargin);
                
        function [cor, lag] = Spike_XCorr(x, y, max_lag_t, t_bin)
            %   returns the unbiased (see doc) cross correlation for x and y (autocorr if
            %   x==y). bins in width t_bin (seconds) if passed. subtracts the mean, and
            %   returns the tau_lag, too.
            
            %   assumes 1000hz sampling is sufficient for spiking data
            
            %   andrew march 31 2010
           
            fs = 1000; % hz
                  
            if ~exist('max_lag_t', 'var')
                max_lag_t = 5;
            end
                        
            if exist('t_bin', 'var')
                fs = t_bin^-1; 
            end
            
            max_lag = round(max_lag_t * fs);  % samples / max_lag_t seconds
            
            t_start = min([min(x), min(y)]);
            t_stop = max([max(x), max(y)]);
            
            x_pp = histc(x, t_start:fs^-1:t_stop);
            y_pp = histc(y, t_start:fs^-1:t_stop);
            
%             x_pp = x_pp - mean(x_pp);
%             y_pp = y_pp - mean(y_pp);
            
            cor = xcorr(x_pp, y_pp, max_lag, 'unbiased'); %, 'unbiased');           
            
            lag = linspace(-max_lag_t, max_lag_t, length(cor));
            
            cor = cor(:);
            
            lag = lag(:);
            
        end  
        
%         function [cor, lag] = CrossCorr(x, y, max_lag_t, t_bin)
%             %   returns the unbiased (see doc) cross correlation for x and y (autocorr if
%             %   x==y). bins in width t_bin (seconds) if passed. subtracts the mean, and
%             %   returns the tau_lag, too.
%             
%             %   assumes 1000hz sampling is sufficient for spiking data
%             
%             %   andrew march 31 2010
%            
%             fs = 1000; % hz
%                   
%             if ~exist('max_lag_t', 'var')
%                 max_lag_t = 5;
%             end
%                         
%             if exist('t_bin', 'var')
%                 fs = t_bin^-1; 
%             end
%             
%             max_lag = round(max_lag_t * fs);  % samples / max_lag_t seconds
%             
%             x_pp = histc(x, min(x):fs^-1:max(x));
%             y_pp = histc(y, min(y):fs^-1:max(y));
%             
%             cor = xcorr(x_pp, y_pp, max_lag, 'unbiased');           
%             
%             lag = linspace(-max_lag_t, max_lag_t, length(cor));
%             
%         end   
        
        function [st_bins, t_bin] = TS2Bins(ts, t_range, dt)
        % [st_bins, t_bin] = CMBHOME.Spike.TS2Bins(ts, t_range, dt) 
        %
        % Converts vector of spike times to vector of binned observations
        % at binsize dt and from times t_range(1) to t_range(2)
        %
        % If t_range and dt are not passed as arguments, dt defaults to 2ms
        % and t_range defaults to the first and last spikes
        %
        % andrew 25 may 2010
        
        if ~exist('t_range', 'var'), t_range = [min(ts) max(ts)]; end % default t_range to first and last spikes
        
        if ~exist('dt', 'var'), dt = .002; end % default dt to 2ms
        
        bin_edges = t_range(1):dt:t_range(2);
        
        st_bins = histc(ts, bin_edges);
        
        st_bins(end-1) = st_bins(end-1) + st_bins(end); % last bin in histc is inclusive values to
        
        st_bins(end) = [];
        
        t_bin = filter(ones(2, 1)/2, 1, bin_edges);
        
        t_bin(1) = [];
        
        end
        
        function watsonsu2 = WatsonsU2(x, y)
        % watsonsu2 = root.HDWatsonsU2(cel)
        %
        % Returns U2 score for vectors of directional data x and y
        %
        % watsonsu2: Watson's U^2 test for uniformity in circular data. See
        % redish,
        % hippocampus 2005 to learn how it applies to head direction cells
        %
        % see
        % http://www.ncbi.nlm.nih.gov/pmc/articles/mid/NIHMS71281/table/T2/
        % for algorithm
        %
        % also see merriam and davis 2008 pg 10
        %
        % a common way to apply this to head direction/spiking data is to
        % pass the array of head direction samples and spike head direction
        % angles. the null hypothesis is thus that spikes head directions come 
        % from the same distribution of all head direction angles. the null 
        % hypothesis is rejected if U^2 is 'large'. >10 seems to be the 
        % neuroscience HD cell literature standard. - andrew

        % ehren newman march 2010

        if iscell(x); x = CMBHOME.Utils.ContinuizeEpochs(x);end
        
        % remove any nan's from x and y
              x(isnan(x)) = [];
              y(isnan(y)) = [];

              if isempty(x) || isempty(y), watsonsu2 = NaN; return; end
              
              z = cat(1, x, y);
              [z,I] = sort(z,'ascend');

              N = length(z);
              n1= length(x);
              n2= length(y);

              z2 = cat(1, ones(n1,1)/n1, -1*ones(n2,1)/n2);
              z2 = z2(I); % sort in same order as z

              dk = cumsum(z2);

              correctIt = 0;
              if correctIt
                % correct for only counting each repeat as coming from one dist
                if n1<n2
                  overlap = sum(ismember(x,y));
                else
                  overlap = sum(ismember(y,x));
                end      
                correction = (1/n2) * ones(overlap,1);

                dbar = (sum(dk) - sum(correction))/N;    
                watsonsu2 = ((n1*n2) / N^2) * (sum( (dk-dbar) .^ 2 ) - sum(correction.^2));

              else

                dbar = sum(dk./N);
                watsonsu2 = ((n1*n2) / N^2) * sum( (dk-dbar) .^ 2 );

              end    
        end            
        
    end
    
    methods
        
        function obj = Spike(varargin)
        % spike_object = Spike(params);
        %
        % Builds CMBHOME Spike Object
        %
        % PARAMS: Spike('param_name', param_value)
        %
        %   i           indeces aligned to root.b_ts (video timestamps). if
        %               passed, vid_ts is ignored
        %   ts          spike timestamps (often much more accurate sampling
        %               than video samples)
        %   vid_ts      timestamps in behavioral, Session object data. used
        %               to align 'ts' above to video samples. uses Utils
        %               function SpkInds
        %   tet_name    name of tetrode or tetrode file... whatevers
        %               important. must be char
        %   prop        struct used for custom fields. whatever you
        %               want. must be a struct
            
            if nargin>0
                
                p = inputParser;

                p.addParamValue('i',   [], @(x) isvector(x));
                p.addParamValue('ts',  [], @(x) isscalar(x)||sum(size(x)~=1)==1) 
                p.addParamValue('vid_ts',  [], @(x) sum(size(x)~=1)==1 || isempty(x));
                p.addParamValue('tet_name',  [], @(x) ischar(x) || isempty(x));
                p.addParamValue('prop', struct([]),@(x) isstruct(x) || isempty(x));

                p.parse(varargin{:});
                
                i = p.Results.i;
                ts = p.Results.ts;
                vid_ts = p.Results.vid_ts;
                tet_name = p.Results.tet_name;
                prop = p.Results.prop;
                
                if isempty(i) && ~isempty(vid_ts)
                    
                    import CMBHOME.Utils.*
                    
                    i = SpkInds(vid_ts, ts); 
                    
                end
                
                if isempty(i), error('Make sure to pass either aligned indices or video timestamps when importing Spike object.'); end

                obj.i = i(:);
                obj.ts = ts(:);
                obj.tet_name = tet_name;
                obj.prop = prop;
            end

        end
        
        function self = set.i(self, i)
            
            i = sort(i);
            
            if isnumeric(i) && all(i>0)
                self.i = i(:);
            else
                error('Improper data format for spike indeces')
            end
        end
        
        function self = set.ts(self, ts)
            
            ts = sort(ts);
            
            if any(diff(ts)==0), warning('CMBH:notify', 'Be aware: duplicate spike timestamps exist. No further action taken.'); end
            
            if isnumeric(ts) %&& all(diff(ts)>0) % must be numeric, but no longer increasing
                self.ts = ts(:);
            else
                error('Improper data format for spike timestamps')
            end
        end
             
    end           
    
end