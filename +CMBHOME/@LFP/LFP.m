% This is a standard data format class for in vivo lfp data at the CMB

classdef LFP
    
    properties
       
        signal          % Raw LFP signal 
        ts              % Timestamps of LFP signals=
        fs              % Sampling rate
        channel_name    % Name of the channel
        user_def        % struct for user defined fields
        
    end
    
    properties (Dependent=true, SetAccess=private)
       
        theta            % LFP signal, filtered to theta range
        theta_phase      % Hilbert phase of theta filtered signal
        theta_amplitude  % Hilbert amplitude of theta filtered signal
        myvar            % Other time varying signal
    end
    
    properties (Hidden)
        
        b_myvar
        b_theta                 % prop to store theta filtered signal
        b_theta_phase           % prop " phase
        b_theta_amplitude       % prop " amplitude
        
    end

    methods
       
function self = LFP(signal, ts, fs, channel_name, b_theta, b_theta_phase, b_theta_amplitude, user_def, b_myvar)            
            
            if ~exist('b_theta', 'var'), b_theta = []; end
            if ~exist('b_theta_phase', 'var'), b_theta_phase = []; end
            if ~exist('b_theta_amplitude', 'var'), b_theta_amplitude = []; end
            if ~exist('user_def', 'var'), user_def = []; end
            if ~exist('b_myvar','var'), b_myvar = []; end
            
            if nargin==0, self.signal = []; self.ts = []; self.fs = []; self.channel_name = []; return; end
            
            if length(signal)~=length(ts), error('signal must be length of ts'); end
            if length(b_myvar)~=length(ts) && ~isempty(b_myvar), error('myvar must be length of ts'); end
            if length(b_theta)~=length(ts) && ~isempty(b_theta), error('theta must be length of ts'); end
            if length(b_theta_phase)~=length(ts) && ~isempty(b_theta_phase), error('theta phase must be length of ts'); end
            if length(b_theta_amplitude)~=length(ts) && ~isempty(b_theta_amplitude), error('theta amplitude must be length of ts'); end
            
            if ~(ischar(channel_name) || iscellstr(channel_name)) && ~isempty(channel_name), error('channel name must be a string'); end
            
            self.signal = signal(:);
            self.ts = ts(:);
            self.b_theta = b_theta;
            self.b_theta_phase = b_theta_phase;
            self.b_theta_amplitude = b_theta_amplitude;
            self.fs = fs;
            self.channel_name = channel_name;
            self.b_myvar = b_myvar;
            self.user_def = user_def;
            
end
        
        function myvar = get.myvar(self) 
            %ehren 2012.06.27
            myvar = self.b_myvar;
        end

        
        function theta = get.theta(self) % retrieves theta filtered LFP signal, from b_theta if it exists. Calculates if not
        % andrew 26 sept 2010
        
            if isempty(self.b_theta)
                
                if ~isempty(self.signal) && ~isempty(self.fs)
            
                    theta = self.ThetaFilter(self.signal, self.fs);
                
                else
                    theta = [];
                end
                
            else
                
                theta = self.b_theta;
                
            end
            
        end
        
        function theta_phase = get.theta_phase(self) % retrieves theta phase LFP signal, from b_theta_phase if it exists. Calculates if not
        % andrew 26 sept 2010
        
            if isempty(self.b_theta_phase)
            
                theta = self.theta;
                
                theta_phase = self.InstPhase(theta);
                
            else
                
                theta_phase = self.b_theta_phase;
                
            end
            
        end  
        
        function theta_amplitude = get.theta_amplitude(self) 
        % retrieves theta amplitude LFP signal, from b_theta_amplitude if it exists. Calculates if not
        % eric zilli 3 nov 2010
        
            if isempty(self.b_theta_amplitude)
            
                theta = self.theta;
                
                theta_amplitude = self.InstAmplitude(theta);
                
            else
                
                theta_amplitude = self.b_theta_amplitude;
                
            end
            
        end
        
        function self = AppendTheta(self)
            % Runs filter and populates LFP.theta
            if numel(self)>1
                for i=1:numel(self)
                    self(i) = self(i).AppendTheta;
                end
            else
                
                if ~isempty(self.b_theta)
                    disp('It appears the theta filtered signal has already been calculated. Clear lfp.b_theta to override.');
                    return;
                end
                
                self.b_theta = self.ThetaFilter(self.signal, self.fs);                
                self = self.AppendThetaPhase;
                self = self.AppendThetaAmplitude;
            end
        end
        
        function self = AppendThetaPhase(self)
            % Runs Hilbert transform and populates LFP.ThetaPhase
            if numel(self)>1
                for i=1:numel(self)
                    self(i) = self(i).AppendThetaPhase;
                end
            else
                
            if ~isempty(self.b_theta_phase)
                disp('It appears theta phase has already been calculated. Clear lfp.b_theta_phase to override.'); 
                return; 
            end
            
            self.b_theta_phase = self.theta_phase;
            end
        end
        
        function self = AppendThetaAmplitude(self)
            % Runs Hilbert transform and populates LFP.Amplitude
            if numel(self)>1
                for i=1:numel(self)
                    self(i) = self(i).AppendThetaAmplitude;
                end
            else
                
                if ~isempty(self.b_theta_amplitude)
                    disp('It appears theta amplitude has already been calculated. Clear lfp.b_theta_amplitude to override.');
                    return;
                end
                
                self.b_theta_amplitude = self.theta_amplitude;
            end
        end
        
        function is = thetaAppended(self)
            % Returns if theta was appended
            is = ~isempty(self.b_theta);
        end
        
        function is = theta_amplitudeAppended(self)
            % Returns if theta amplitude was appended
            is = ~isempty(self.b_theta_amplitude);
        end
        
        function is = theta_phaseAppended(self)
            % Returns if theta phase was appended
            is = ~isempty(self.b_theta_phase);
        end
        
    end
    
    methods (Static)
        
        signal_filtered = BandpassFilter(signal, Fs, Fpass)
        
        signal_filtered = DeltaFilter(signal, Fs)
        
        signal_filtered = ThetaFilter(signal, Fs)
        
        signal_phase = InstPhase(signal)
        
        signal_phase = InstAmplitude(signal)
        
    end  
    
end