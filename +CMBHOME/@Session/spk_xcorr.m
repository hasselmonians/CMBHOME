 function [spk_xcorr, lags, epoch, std_error] = spk_xcorr(self, cells, max_lag, t_bin, average_epochs, norm)
    % Returns the unbiased crosscorrelation, or autocorrelation
    %
    % ARGUMENTS
    %   cells               a Nx2 (N=1 or 2) array of cells for which to
    %                       calculate xcorr. if N=1, than autocorrelation is returned
    %   max_lag             maximum lag for which xcorr is calculated. epochs that are not this long are removed
    %   t_bin               binsize in seconds
    %   average_epochs      1 or 0 (default 0). If 1, cor is an array
    %                       calculated by the average of correlations for each epoch
    %                       weighted by epoch duration. If 0, cor is a
    %                       cell array of correlation vectors for each
    %                       epoch in the return variable epoch, which
    %                       omits root.epoch for those shorter than
    %                       max_lag, and marges any overlapping epochs
    %                       in root.epoch so that we don't double count
    %
    % RETURNS
    %   cor                 either a vector (if one epoch, or average_epochs==1, or a column array where each column is 
    %                       the xcorr signal for epochs with spikes and
    %                       as long as max_lag
    %   lag                 same as above, but lags in seconds
    %   epoch               returns epochs used in analysis (in case
    %                       any of root.epoch were < max_lag, or had no
    %                       spikes
    %   std_error           if average_epochs=1, then the average error
    %                       for each lag timepoint is returned for the averaged cor vector
    %
    % [cor, lag, epoch, std_error] = root.spk_xcorr(cel, max_lag,t_bin, average_epochs);


    import CMBHOME.* % has the spk_xcorr function
    import CMBHOME.Utils.*

    std_error = [];

    if nargin<2
        error('spk_xcorr: Not enough arguments');
    end

    if nargin<4
        t_bin=.001;
    end

    if nargin<3
        max_lag = 5; %seconds
    end

    if size(cells,2)~=2
        error('spk_xcorr: cells must be like [tetrode ind, cell ind]');
    end

    if size(cells,1)==1
        cells(2,:) = cells(1,:);
    end

    if ~exist('norm', 'var'), norm = 'count'; end

    if ~exist('average_epochs', 'var'), average_epochs = 0; end

    if ~exist('speed_dur', 'var'), speed_dur = 1; end %#ok<NASGU>

    if average_epochs, self = MergeEpochs(self); end % if we are averaging epochs, make sure that we aren't double counting

    epoch = self.epoch;

    tooshort = ( epoch(:,2)-epoch(:,1) < max_lag ); % remove any epochs which are shorter than max_lag
    self.epoch(tooshort, :) = []; epoch = self.epoch;

    if average_epochs, self = MergeEpochs(self); end % if we are averaging epochs, make sure that we aren't double counting

    if isempty(self.epoch)

        spk_xcorr = []; lags = []; epoch = []; std_error = [];

        warning('CMBH:error', 'No epochs meet this running speed requirement');

        return

    end

    %x = self.cel_ts(cells(1,:));
    %y = self.cel_ts(cells(2,:));
    self.cel=cells(1,:); x = self.cel_ts;
    self.cel=cells(2,:); y = self.cel_ts;
    
    if iscell(x)

        empty_inds = cellfun(@isempty, x); % find empty epochs and remove them

        empty_inds = empty_inds+cellfun(@isempty, y);

        empty_inds = logical(empty_inds);

        x(empty_inds) = []; % remove empty epochs

        y(empty_inds) = []; % remove empty epochs

        if isempty(x) % if no epochs with spikes

            lags = zeros(nbins, 1);

            spk_xcorr = zeros(nbins, 1);

            return

        end                    

        epoch = epoch(~empty_inds, :); % remove epochs without spikes

        nbins = length(-max_lag+t_bin/2:t_bin:max_lag-t_bin/2);

        spk_xcorr = zeros(nbins, length(x));

        lags = zeros(nbins, length(x));

        for i = 1 : length(x)
            [spk_xcorr(:, i), lags(:, i)] = Spike.CrossCorr(x{i}, y{i}, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', norm);
        end

        if average_epochs % cross correlation across all epochs

            spk_xcorr = sum(spk_xcorr, 2);

            lags = lags(:,1);


        end

    else          

        [spk_xcorr, lags] = Spike.CrossCorr(x, y, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', norm);

    end

end