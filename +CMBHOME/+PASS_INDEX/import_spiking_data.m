function [pos_ts,pos,lfp_ts,lfp_sig,spk_ts] = import_spiking_data(filepath, tetrode, unit)
% [pos_ts,pos,lfp_ts,lfp_sig,spk_ts] = import_spiking_data(filepath, tetrode, unit)
%
% Extracts behavioral, lfp and spiking data from filepath.
%
% Inputs:
%   filepath: Full path to the directory for the target session
%   tetrode: Which tetrode to search (used to find corrext .ncs & .ntt'
%   unit: Which cell(cluster) from the ntt to use
% 
% Outputs:
%   pos_ts: Nx1 vector of time stamps from behavioral data (seconds)
%   pos: Nx2 vector of [x y] from video (pixels)
%   lfp_ts: Nx1 vector of lfp time stamps (seconds)
%   lfp_sig: Nx1 vector of lfp signal (unfiltered)
%   spk_ts: Nx1 vector of time stamps where desired unit fired (seconds)
%
%  2013-09-01 - Written by G. William Chapman IV
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com) 

%% set up & Check for existence 
% Check that videoFile and tetrode file exist
videoFile = fullfile(filepath,'VT1.nvt');
lfpFile = fullfile(filepath,['CSC' num2str(tetrode) '.ncs']);
spkFile = fullfile(filepath,['TT' num2str(tetrode) '.ntt']);
spkTxtFile = fullfile(filepath,['TT' num2str(tetrode) '.txt']);

if ~(exist(videoFile,'file') && exist(lfpFile,'file') && exist(spkFile,'file'))
    error('At least one of your files does not exist') % Throw error if files not present
end

%% Import the behavioral data (from Video)
FieldSelectionArray = [1 1 1 1 0 0];
ExtractHeader = 0;
ExtractionMode = 1;
[ts, x, y,~] = Nlx2MatVT(videoFile, FieldSelectionArray,ExtractHeader,ExtractionMode);
x = x(:); y = y(:); ts = ts(:);
[~, x, y] = FixPos(ts,x,y,10);

pos_ts = ts / 1000000; %Convert to seconds from microseconds
pos = [x y];

%% Import the lfp data
[tmp_ts, lfp_sig] = Nlx2MatCSC(lfpFile,[1 0 0 0 1], 0, 1);
lfp_ts = NaN(size(lfp_sig,1) * size(lfp_sig,2),1);
tmp_ts(end+1) = tmp_ts(end) + mean(diff(tmp_ts));
tmp_ts = tmp_ts(:);

for i = 1:length(tmp_ts)-1
    dt = ((tmp_ts(i+1)-tmp_ts(i)) / size(lfp_sig,1)) * (0:1:size(lfp_sig,1));
    dt = dt(1:end-1);
    lfp_ts( (i-1) * size(lfp_sig,1) + 1 : i * size(lfp_sig,1)) = tmp_ts(i) + dt;
end

lfp_ts = lfp_ts / 1000000;

lfp_sig = lfp_sig(:); 

%% Import spike data
% Uses either .ntt files or exported txt files to idenify spike times and
% which spikes belong to cell of interest.
%
% Format of text file: [cell_ind spike_time] 

if exist(spkTxtFile,'file') % Look for exported text files first
    spkTxt = load(spkTxtFile);
    cell_ind = spkTxt(:,1);
    spk_ts = spkTxt(:,2);
else                          % Resort to .ntt files
    FieldSelectionArray = [1 0 1 0 0];
    ExtractHeader = 0;
    ExtractionMode = 1;
    [spk_ts, cell_ind] = Nlx2MatSpike(spkFile, FieldSelectionArray,ExtractHeader,ExtractionMode);
    spk_ts = spk_ts / 1000000;
end

spk_ts = spk_ts(cell_ind == unit); %Filter down to the desired unit

end


function [ts,x,y] = FixPos(ts,x,y, max_allowed_flips)
% [ts,x,y] = FixPos(ts,x,y, max_allowed_flips)
%
% Takes all 0,0s and large jumps in position (greater than
% jitter_threshold pixels/sample) that persist for less than
% max_allowed_flips (default = 5) and linearly interpolates the missing
% data. Smooths conservatively afterward, as well (convolution with a gaussian, standard
% deviation = 2 samples).
%
% andrew december 2009
% update andrew june 15 2011

warning('off', 'MATLAB:interp1:NaNinY');

jitter_threshold = 20;

bads = (x==0 | y==0);

x(bads) = NaN;
y(bads) = NaN;

if ~exist('max_allowed_flips', 'var')
    max_allowed_flips = 5; % samples
end

flips = findOnsetsAndOffsets(isnan(x));

flips(:,2) = flips(:,2)+1;

flips = cat(1, 1, flips(:), find([0; sqrt(diff(x).^2 + diff(y).^2)]>jitter_threshold), length(x));

flips = sort(unique(flips)); % indeces of NaNs or jumps

flips = [flips(1:end-1), flips(2:end)];  % epochs formation
    
flips(:,2) = flips(:,2)-1; % adjust for diff shift

flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = [];

flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps

flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices

tso = ts;

x([flips{:}]) = []; % remove samples in ts and x
y([flips{:}]) = []; % remove samples in ts and x
ts([flips{:}]) = [];

x = interp1(ts, x, tso);
y = interp1(ts, y, tso);

x = ndnanfilter(x, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
y = ndnanfilter(y, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);


end

function [OnOffs] = findOnsetsAndOffsets(boolVec)
% function [OnOffs] = findOnsetsAndOffsets(boolVec)
%
% Returns list of aligned start and stops of chunks of 1's in a vector of
% 1's and 0's.
%
% INPUTS
%  boolVec - vector of 1's and 0's 
%
% OUTPUTS
%  startEnds - Nx2 list of indices of the first and last 1's for the N
%              contiguous blocks of 1's.

% created 6/16/11 arb

boolVec = boolVec(:)';

starts = find(diff(boolVec)==1);
ends = find(diff(boolVec)==-1);

% if starts out going correct speed, add 1 to starts
if boolVec(1)
  starts = [0 starts];
end

% if finishes going correct speed, add final value to ends
if boolVec(end)
  ends = [ends length(boolVec)];
end

OnOffs = [starts(:)+1 ends(:)];
end
