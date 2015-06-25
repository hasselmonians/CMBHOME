function [epochs, sig, p, F, df1, df2 ] = sig_theta( signal,varargin )
%[epochs, sig, p, F, df1, df2] = CMBHOME.LFP.sig_theta(root.lfp.signal,
%'ts', root.lfp.ts, 'FS', root.lfp.fs);
%
% INPUT
%   signal - The signal
% 
% PARAMETERS
%   ts ([]) - The time stamps for the signal (sec). If [], generates from Fs.
%   Fs ([]) - The sampling frequency (Hz). If [], generates from ts or uses
%       250
%   low_f (2) - The cutoff for "line noise" highpass filter.
%   examination_window (0.5) - The width of the examination window (secs).
%   sigLevel (0.05) - Rejection region for test
%   correction ('benjamini-hochberg') - How to correct for multiple
%       comparisons. Can be 'none', 'bonferonni', or 'benjamini-hochberg'
% 
% RETURNS
%   sig - The test at all time stamps in the signal, padded with NaNs where
%       the window runs into the edge.
%   p - The P value of the F-tests
%   F - The F values
%   df1 - The first degrees of freedom
%   df2 - The second degrees of freedom
%
% 2014-05-12
% Jason Climer - jason.r.climer@gmail.com

%% Parse input
signal = signal(:);

ip = inputParser;
ip.addParamValue('ts',[]);
ip.addParamValue('Fs',[]);
ip.addParamValue('low_f',2);
ip.addParamValue('examination_window',0.5);
ip.addParamValue('sigLevel',0.05);
ip.addParamValue('correction','benjamini-hochberg');
ip.parse(varargin{:});
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

if isempty('Fs')
   if ~isempty(ts), Fs = 1/mean(diff(ts)); else Fs = 250; end
end

if isempty(ts), ts = ((1:numel(signal))-1)/Fs;end


%%
% Filtr out line noise (adds a lot of RSS)
[bhigh,ahigh] = butter(3,low_f/(Fs/2),'high');
signal_line = filtfilt(bhigh,ahigh,signal);

% Do FFTs in overrepresented spectrogram
L = round(examination_window*Fs/2)*2;
NFFT = 2^nextpow2(L);
y = signal_line(repmat((1:L)'-1,[1 numel(signal)-L])+repmat(1:numel(signal)-L,[L 1])).*repmat(hann(L),[1 numel(signal)-L]);
Y = fft(y,NFFT);
Y = Y(1:NFFT/2+1,:);

% Get PSD
faxis = Fs/2*linspace(0,1,NFFT/2+1);
P = (1/(Fs*NFFT))*abs(Y).^2;
P(2:end-1,:) = 2*P(2:end-1,:);

% Resudual sum of squares and F tests
S1 = sum(P);
S2 = sum(P(faxis<4|faxis>13,:));

df1 = sum(faxis>=4&faxis<=13)*2;
df2 = L-df1-1;

F = ((S1-S2)/df1)./...
    (S2/df2);

p = 1-fcdf(F,df1,df2);

% Correct for multiple comparisons
switch correction
    case 'none'
        %
        sig = p<sigLevel;
    case 'bonferonni'
        %
        sig = p<sigLevel/numel(signal);
    case 'benjamini-hochberg'
        %
        sig = sort(p);sig = sig(:);
        i = 1:numel(sig);i = i(:);
        i = find(sig<i/numel(sig)*sigLevel,1,'last');
        sig = p<sig(i);
end

sig = [NaN(L/2,1);sig(:);NaN(L/2,1)];

% 
%% epoching
epochs = CMBHOME.Utils.OverThresholdDetect(sig,0.5,1,1)/Fs;