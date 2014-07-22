function filtered_signal = BandpassLFP(self, ind, varargin)
% (1) instantaneous_phase = root.InstantaneousPhase(ind);
% (2) instantaneous_phase = root.InstantaneousPhase(ind, params);
%
% ARGUMENTS:
% 'ind' required, index in root.b_lfp. Cannot exceed
% length(root.b_lfp)
% 'Fpass' 1x2 vector defining passband [f_low, f_high]
%
% (1) Returns vector with filtered LFP signal that exists in root.b_lfp(ind).signal,
% where ind is an integer and is less than or equal to size of structure
% array root.b_lfp.
%
% (2) takes parameters:
% 'field' - include a string referencing a
% field in root.b_lfp that will be analyzed instead of
% root.b_lfp(ind).signal
%
% 'Fpass' - include a 1x2v vector [f_low f_high] to bandpass signal. Signal
% is root.b_lfp(ind).signal unless 'field' property above is passed.
%
% 'Passband' - either 'theta' or 'delta' or 'ripple' to bandpass filter the signal.
% Define your own passband with Fpass
%
% Remember that root.fs is the LFP sampling frequency, and this function
% depends on it.
%
% andrew 17 may 2010

import CMBHOME.*
import CMBHOME.LFP

p = inputParser;

p.addRequired('self');
p.addRequired('ind');
p.addParamValue('field', 'signal', @(x) ischar(x)); 
p.addParamValue('Fpass',  [], @(x) numel(x)==2);
p.addParamValue('Passband', [], @ischar);

p.parse(self, ind, varargin{:});

field = p.Results.field;
Fpass = p.Results.Fpass;
Passband = p.Results.Passband;

if ind>length(self.b_lfp)
    error('ind exceeds length of root.b_lfp');
end

if isfield(self.b_lfp(ind), field)
    signal = self.b_lfp(ind).(field);
else
    error('Field does not exist in root.b_lfp');
end

if Passband
    switch lower(Passband)
        case 'theta'
            Fpass = [6 10];
        case 'delta'
            Fpass = [2 4];
        case 'ripple'
            Fpass = [80 250];
        otherwise
            error('Passband can be either theta or delta or ripple. Use Fpass to define your own passband');
    end
end

filtered_signal = []; % initialize returns

if isempty('Fpass', 'var')
    help CMBHOME.Session.BandpassLFP
    error('You must either define Fpass or Passband.');
end

if ind>length(self.b_lfp) 
    help CMBHOME.Session.BandpassLFP
    disp('ind exceeds number of lfp signals');
    return
    
elseif ~isfield(self.b_lfp(ind), field) % if signal doesnt exist
    help CMBHOME.Session.BandpassLFP
    disp([field 'is not a field in b_lfp.']);
    return
    
elseif isempty(self.b_lfp(ind).(field)) % if signal is empty
    help CMBHOME.Session.BandpassLFP
    disp('signal does not exist for this index');
    return
    
end   

filtered_signal = LFP.BandpassFilter(self.b_lfp(ind).signal, self.fs, Fpass);