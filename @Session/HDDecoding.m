function decoded_hd = HDDecoding(self, cels, tr_epochs, de_epochs, varargin)
% root.HDDecoding(cels, tr_epochs, de_epochs, varargin)
%
% Performs Head Direction Reconstruction.
% The model is trained during tr_epochs and reconstructed during de_epochs,
% unless 'TestModel' is set to 1, in which case training epochs are reconstructed as
% well.
%
% Arguments:
%
% Returns:
% decoded_hd -> an Nx4 cell array of {Prln, ts, bins, Ft, hd; ...} for as many
% epochs are decoded. Ft is matrix of cell firing rates for each decoding
% epoch (t x cells)
%
% Bayes method implemented (it's the only non-vector algorithm which is insensitive
% to bimodal tuning curves). Continuity constraint is such that the standard
% deviation if confined to 
%
% see Johnson, Seeland and Redish 2005 (Hippocampus)
% coded as per Zhang et al. 1998
%
% andrew 26 may 2010

import CMBHOME.Utils.*

p = inputParser;

p.addRequired('self');
p.addRequired('cels');
p.addRequired('tr_epochs');
p.addRequired('de_epochs');
p.addParamValue('ScreenPrint', 1, @(x) numel(x)==1);
p.addParamValue('SavePDF', 0, @(x) numel(x)==1);
p.addParamValue('SavePrefix', 'CMBHDecoding', @ischar);
p.addParamValue('Tau',  1, @(x) numel(x)==1);
p.addParamValue('BinWidth', 10, @(x) numel(x)==1);
p.addParamValue('TestModel', 0, @(x) numel(x)==1);
p.addParamValue('Coherency', 0, @(x) numel(x)==1);

p.parse(self, cels, tr_epochs, de_epochs, varargin{:});

ScreenPrint = p.Results.ScreenPrint;
SavePDF = p.Results.SavePDF;
SavePrefix = p.Results.SavePrefix;
tau = p.Results.Tau;
bin_angle = p.Results.BinWidth;
test_model = p.Results.TestModel;
Coherency = p.Results.Coherency;

load '/Users/abogaard/repos/abogaard/+CMBHOME/+Utils/cmap.mat'; % load colormap

epsilon = .001; % arbitrarily small bump to tuning curves so probability doesnt blow up when a cel fires and the tuning map = 0

bins = -180+bin_angle/2 : bin_angle : 180-bin_angle/2;

dt = .1; % seconds, window increment

C = []; % initialize return variables
dist_I = [];

if ~exist('tr_epochs', 'var')
    tr_epochs = cat(2, (self.epoch(1):30:self.epoch(2)-15)', (self.epoch(1)+15:30:self.epoch(2))'); % training epochs
    de_epochs = cat(2, (self.epoch(1)+15:30:self.epoch(2)-15)', (self.epoch(1)+30:30:self.epoch(2))');% decoding epochs
end

%%%%%%%% Training phase

self.epoch = tr_epochs;

if ~exist('cels', 'var'), cels = self.cells; end

[mFr, theta] = self.DirectionalTuningFcn(cels, 'binsize', bin_angle, 'Continuize', 1); % bins x cells, tuning functions (mFr = model F(r))

mFr = mFr + epsilon;

Pr = self.DirectionalOccupancy(bin_angle, 1) ./ sum(tr_epochs(:,2)-tr_epochs(:,1)); % vector bins x 1, probability of animal heading for all training epochs

if Coherency, [dist_I, I_bins] = BootstrapIncoherency(self, mFr, cels, tr_epochs, bins, tau, dt); end % distribution of I x angle bins padded by NaNs

%%%%%%%% Decoding phase

decoded_hd = cell(size(de_epochs,1),7);

if test_model, de_epochs = cat(1, tr_epochs, de_epochs); end % perform decoding on training epochs, as well

if ScreenPrint, InitFigures(size(de_epochs, 1) + size(mFr, 2)); end % initialize figures

for i = 1:size(de_epochs,1) % run through all decoding epochs
    
    epoch(:,1) = de_epochs(i,1):dt:de_epochs(i,2)-tau;
    epoch(:,2) = de_epochs(i,1)+tau:dt:de_epochs(i,2);
           
    self.epoch = epoch;
    
    N = cellfun(@numel, self.spk_ts(cels)); % matrix epochs x cells of spike counts
    
    expon = repmat( N, [1, 1, length(bins)]); % number of spikes each cell fired, epochs x cells x spike count repeated as per # bins

    Fr = repmat( reshape(mFr', [1, size(mFr,2), size(mFr, 1)]), [size(expon, 1), 1, 1]); % epochs x cells x bins
    
    Prln = prod(Fr .^ expon, 2) .* exp(-tau * sum(Fr, 2)); % time (epochs) x 1 x prob (bins), conditional probability of head direction given spiking activity

    Prln = repmat(Pr, [1 size(Prln,1)]) .* shiftdim(Prln, 2); ...
    Prln = Prln ./ repmat(sum(Prln, 1), size(Prln,1), 1); % reshape to prob (bins) x time for plotting & multiply by probability of head direction

    hd = ProbableHD(Prln, bins);
    
    if Coherency
        I = Incoherency(self, hd, bins, mFr, self.epoch, cels);
        C = CoherencyFcn(dist_I, I, hd, bins, I_bins);
    end
    
    decoded_hd(i,:) = {Prln, mean(self.epoch,2), bins, N/tau, hd, C, dist_I}; % normalize and save decoding
    
    if ScreenPrint, UpdateFigures(self, h_fig, h_axes, decoded_hd{i}); end
    
    clear epoch
    
end

if ScreenPrint, for j = 1:size(mFr,2), set(h_fig, 'CurrentAxes', h_axes(i+j)); plot(bins, mFr(:,j)); end; %% also plot tuning curves

ManageSP(h_fig, h_axes, 'MaxCols', 3, 'MaxRows', 5, 'Save', 1, 'CMap', cmap, 'Save', SavePDF, 'SavePrefix', SavePrefix); % in CMBHOME.Utils.*

end

end

function UpdateFigures(self, h_fig, h_axes, Prln)

set(h_fig, 'CurrentAxes', h_axes); % make current our axes
    
imagesc(self.epoch(:,1), bins, Prln), hold on % plot probability distribution

self.epoch = [self.epoch(1,1), self.epoch(end,2)]; % reset epoch so that we can see all headdir samples

line(self.ts, self.headdir, 'Color', 'r', 'LineWidth', 2); % plot recorded head direction

end

function [h_fig, h_axes] = InitFigures(nAxes)

h_axes = zeros(nAxes, 1);

h_fig = figure('Visible', 'off'); % initialize plotting functions

for i = 1:nAxes, h_axes(i) = axes('Visible', 'off'); end % make axes object array

end

function vec_hd = ProbableHD(Prln, bins)

    [~, ind] = max(Prln, [], 1);

    vec_hd = bins(ind);

end

function I = Incoherency(root, hd, bins, mFr, de_epochs, cells)
% A subfunction within HDDecoding which estimates the reconstruction
% accuracy, or error angle, as well as the Representational Quality, or
% coherence (as coined by Johnson and Redish).
%
% The Reconstruction Error is the difference between decoded
% angle and actual angle. Both median error and absolute median error are
% returned. Median error is good for bias, while absolute median error is
% good to determine the accuracy.
%
% Arguments:
%
% mFr -> bins x cells matrix of orientation tuning function
%
% de_epochs -> Nx2 array of timestamps for decoded epochs
%
% cells -> Mx2 array of cells
%
% Returns:
%
% I ->  vector of N elements (size(de_epochs,1)), indicating metric for
%       difference between expected and measured activity packets
%
% The Reconstruction Coherency is the probability that the actual and
% expected activity of the system are the same. Algorithm below:
%
% Calculate activity packets for actual and expected as the sum of the
% tuning curves of the cells, weighted by their firing rates. The actual
% activity packet is based on the actual firing rates, while the expected
% activity packet is based onthe expected firing rates, given the
% constructed orientation of the population at time t.
%
% Solve for incoherency between activity packets at each moment by taking:
%
%         angvar(A(phi,t) - A_est(phi,t))
% I(t) = -------------------------------- ,
%          integral( A_est(phi, t) dphi
%
% where angvar is the variance of the difference between activity packets
% of phi at each time point t.
%
% andrew 28/6/2010
    
% for every window, and thus estimation of hd, find number of spikes,
% divide by seconds for F_f, weight the tuning curves, and make the
% activity packet A

% BUILD A

root.epoch = de_epochs;

spikes = cellfun(@length, root.spk_ts(cells)); % matrix of spike counts epochs x cells

f = spikes ./ repmat(de_epochs(:,2) - de_epochs(:,1), 1, size(spikes,2)); % firing frequency epochs x cells

A = repmat(mFr, [1, 1, size(de_epochs,1)]) .* ... % weight the tuning curves by f 
    repmat(reshape(f', 1, size(f, 2), size(f, 1)), size(mFr,1), 1); % A is bins x cells x epoch

A = reshape( sum(A,2) ./ repmat(sum(mFr, 2), [1 1 size(A, 3)]), size(A, 1), size(A, 3) ); % A is bins x epoch

% BUILD Aexp

tmp = repmat(hd(:)', length(bins), 1) - repmat(bins(:), 1, length(hd));
    
[~, tmp] = min(abs(tmp), [], 1);

hd = bins(tmp);

f = mFr(tmp, :); % mFr is bins x cells, f is epochs x cells

Aexp = repmat(mFr, [1, 1, size(de_epochs,1)]) .* ... % weight the tuning curves by f 
    repmat(reshape(f', 1, size(f, 2), size(f, 1)), size(mFr,1), 1); % A is bins x cells x epoch

Aexp = reshape( sum(Aexp,2) ./ repmat(sum(mFr, 2), [1 1 size(Aexp, 3)]), size(Aexp, 1), size(Aexp, 3) ); % Aexp is bins x epoch

% SOLVE FOR I

I = var(A - Aexp) ./ sum(Aexp,1); % I is 1 x epochs Incoherency

%I = sqrt(sum((A - Aexp).^2, 1))./sum(Aexp,1); % RMSE

end

function [dist_I, I_bins] = BootstrapIncoherency(self, mFr, cels, tr_epochs, bins, tau, dt)
% for all training epochs, we calculate  the average head direction in
% windows of lngth tau, just as in the decoding phase, calculate A and
% Aexp, and come up with a distribution of I for every average head
% direction in bins as per bin_angle
%
% Returns dist_I, which is the CDF of I for training epochs

HDvI = cell(size(tr_epochs,1), 1);

for i = 1:size(tr_epochs,1) % run through all decoding epochs
    
    epoch(:,1) = tr_epochs(i,1):dt:tr_epochs(i,2)-tau;
    epoch(:,2) = tr_epochs(i,1)+tau:dt:tr_epochs(i,2);
           
    self.epoch = epoch;
    
    HD = self.headdir; % cell array of 1 x epochs of headdirection vectors
    
    HD = cellfun(@deg2rad, HD, 'Unif', false); % convert to radians
    
    HD = cellfun(@(c) atan2(mean(sin(c)), mean(cos(c))), HD); % average HD, epochs x cells matrix
    
    HD = HD*180/pi; % vector of average HDs
    
    I = Incoherency(self, HD, bins, mFr, self.epoch, cels);
    
    tmp = repmat(HD(:)', length(bins), 1) - repmat(bins(:), 1, length(HD));
    
    [~, tmp] = min(abs(tmp), [], 1);
    
    HDn = bins(tmp);
    
    HDvI{i} = cat(2, HDn(:), I(:));
    
    clear epoch

end

HDvI = vertcat(HDvI{:});

I_bins = 0:.001:3;

dist_I = zeros(length(I_bins), length(bins));

for i = 1:length(bins)
   
    I = HDvI( HDvI(:,1)==bins(i), 2 ); % cell array of I vectors indexed by angle bin  
    
    ['max and min I' num2str([min(I), max(I)])]; % for checking purposes

    I = smooth( hist(I, I_bins) , 10);  % smooth I once

%     figure
%     plot(I_bins, I)
    
    dist_I(:,i) = cumsum(I) / sum(I); % smooth I once CDF
       
%     figure
%     plot(I_bins, dist_I(:,i))
%     
%     keyboard
end

end

function C = CoherencyFcn(dist_I, I, hd, bins, I_bins)
% Interpolates the probability of I for each estimated HD sample according
% to the CDFs in dist_I, indexed by bins

tmp = repmat(hd(:)', length(bins), 1) - repmat(bins(:), 1, length(hd));
    
[~, ind] = min(abs(tmp), [], 1);

C = zeros(length(hd), 1);

for i = 1:length(ind)
   
    C(i) = 1 - interp1(I_bins, dist_I(:,ind(i)), I(i));
    
end

end