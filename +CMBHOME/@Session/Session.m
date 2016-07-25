% This is a standard data format class for in vivo spiking, video tracking,
% and lfp data at the CMB at Boston University. See
%
% ACCESSING DOCUMENTATION
%
%   >> doc CMBHOME.Session; % for documentation on root objects
%   >> doc CMBHOME.LFP; % for documentation on LFP objects
%   >> doc CMBHOME.Spike; % for documentation on Spike objects (units)

classdef Session
    
    properties (Hidden)
        
        ops = {'ge', 'le'}  % operators for epochs behavior
        
    end
    
    properties (Access=private)
       
        p_ind               % same as root.ind, but set by set.epoch
        p_lfp_ind
        p_cel_ind           % indices in root.b_ts and root.b_lfp(active_lfp).ts for each cell spike. set with set.cel
        p_cel_spkind        % indices in root.spike(cel(1), cel(2)).ts for each spike in the epoch. set with set.cel
        p_cel_lfp_ind
        
    end
    
    properties
        
        date_created        % date and time created as a serial date number (see 'doc now')
        notes               % cell array of notes and dates {datestamp, str_note}
        spike               % object array of class 'Spike', contains unit data  
        name                % 1xM cellstr of directory and filename of current object.
        cell_thresh = [0,0] % threshold of minimum spike rate [minimum spike frequency, check drift]
        fs_video            % sampling rate from video (samples/sec)
        raw_pos=1           % boolean, is this uncorrected position?
        raw_headdir=1       % boolean, is this uncorrected head direction?
        event = cell(1, 2)  % event flags, Nx2 cell array, {'label', time}            
        path_lfp            % cell array of fnames to lfp recording data
        path_raw_data       % path to original recording data
        spatial_scale=.5    % cm/pixel
         
        rotate=0            % deg rotation
        shift=[0 0]         % cm translation
        
        epoch_group         % vector, Nx1, N is number of epochs with integer indicating like-epochs
        user_def            % property that may be whatever user defines, ex. struct of relevant
        version = 1.5         % version control
        b_lfp               % object array of structs with fields ts, signal, path_lfp_index ex
        myvar2
    end
    
    properties (Access=private, Hidden)
        
        p_epoch
        p_cel
        p_active_lfp
        p_b_x               
        p_b_y                 
        p_b_vel               
        p_b_headdir           
        p_b_ts               
        p_b_myvar
        
    end
    
    properties (Dependent)
        
        epoch               % array, Nx2, start and stop times that exists (can be inf or -inf)
        cel                 % Nx2 array of [tetrode, cell index]
        active_lfp          % index to root.path_lfp indicating which LFP is loaded (or to be loaded)
        
    end
    
    properties (Dependent, Hidden)
    
        b_x                 % x position
        b_y                 % y position
        b_vel               % user defined velocity, if undefined, the derivative of b_x and b_y is used
        b_headdir           % head direction
        b_ts                % timestamps
        b_myvar             % other time varying variable
    
    end
    
    properties (Dependent=true, SetAccess=private)
        
        vel                 % velocity, within start and stop epoch (pixels/sec)
        svel                % velocity, (cm/sec)
        x                   % x position (pixels)
        y                   % y position (pixels)
        sx                  % x position (cm)        
        sy                  % y position (cm)
        sheaddir            % headdir (degrees,-180:180)
        ind                 % index in root.b_ts
        ts                  % timestamps (seconds)
        myvar               % Arbitrary time varying variable
        b_epoch_group       % vector, Nx1, N is number of epochs with integer indicating like-epochs
        headdir             % head direction, (degrees, -180:180)
        cells               % Nx2 array of cells which meet specs_cell above,  [tetrode, cell;...]
        active_event        % an Nx2 element cell array with event labels (if exist) that match current epoch times.
                            % N is the number of rows in root.epochs
        lfp                 % object array of LFP objects of fields ts and signal
        spk                 % struct with fields representing data at spike times of root.cel
        name_formatted      % this is the root.name field formatted for the current platform
        cel_x               % x position at spike times of root.cel
        cel_sx              % scaled x position at spike times of root.cel
        cel_y               % y position at spike times of root.cel
        cel_sy              % y position at spike times of root.cel
        cel_vel             % velocity at spike times of root.cel
        cel_svel            % scaled velocity position at spike times of root.cel
        cel_ts              % timestamps at spike times of root.cel
        cel_headdir         % headdir position at spike times of root.cel
        cel_sheaddir        % headdir position at spike times of root.cel
        cel_i               % indices to root.b_VAR at spike times of root.cel
        cel_theta           % Theta filtered LFP power at spike times
        cel_thetaphase      % Thetaphase position at spike times of root.cel
        cel_lfpmyvar        % root.b_lfp(active_lfp).b_myvar at spiketimes of root.cel
        cel_myvar           % root.myvar at spike times of root.cel
        cel_myvar2          % root.myvar2 at spike times of root.cel
                
    end
    
    methods (Static)
        
        [t,x,y,vx,vy,ax,ay] = KalmanVel(posx,posy,post,order,Q,R);
        
        function root = MergeSessions(cellarray_filename) %#ok<INUSD>
        % (1) root = CMBHOME.Session.MergeSessions
        % (2) root = CMBHOME.Session.MergeSessions(cellarray_filenames)
        %
        % Merges multiple Session objects by stitching together time
        % indeces. Saves original Session times as events
        %
        % (1) Prompts user to select files to merge, stiches them together
        % in the order selected
        % (2) Accepts a cell array of file names to be merged. Stitches
        % them together in the order of cellarray_filenames
        %
        % andrew 24 may 2010
        
        warning('CMBH:notify', '%s', 'Warning: If you want to merge the LFP as well you must load it, resave the files, then run this function');
        
        import CMBHOME.*
        f_ind=1;
        if ~exist('cellarray_filenames', 'var')
            while 1 % loop until user cancels, and we get no files

                [load_files, base_path] = uigetfile('*.mat','Please select Session object mat files to merge. Exit to finish.',pwd, 'MultiSelect', 'on');

                if iscell(load_files)
                        for i = 1:length(load_files)
                            fopen(f_ind) = {fullfile(base_path, load_files{i})}; %#ok<AGROW>
                            f_ind = f_ind+1;
                        end
                elseif ischar(load_files)
                        fopen(f_ind) = {fullfile(base_path, load_files)}; %#ok<AGROW>
                        f_ind = f_ind+1;
                else
                    break; % exit loop collecting files
                end

            end
        else
            fopen = cellarray_filenames;
        end
        
        clear x, clear y, clear ts, clear hd, clear vel
        clear spike, clear even, clear event, clear lfp
        for i = 1:length(fopen)
            tmp = load(fopen{i});
            x{i}=tmp.root.b_x;
            y{i}=tmp.root.b_y;
            ts{i}=tmp.root.b_ts;
            hd{i}=tmp.root.b_headdir;
            vel{i}=tmp.root.vel;
            spike{i}=tmp.root.spike;
            event{i}=tmp.root.event;
            lfp{i}=tmp.root.b_lfp;
            fs = tmp.root.fs_video;
            ss = tmp.root.spatial_scale;           
        end
        
        x=CMBHOME.Utils.ContinuizeEpochs(x(:)); y=CMBHOME.Utils.ContinuizeEpochs(y(:));
        hd=CMBHOME.Utils.ContinuizeEpochs(hd(:)); vel=CMBHOME.Utils.ContinuizeEpochs(vel(:));
        event=CMBHOME.Utils.ContinuizeEpochs(event(:));
        
        %ts
        firsts = cellfun(@(x) x(1), ts);
        ts=cellfun(@(x) x-x(1), ts,'UniformOutput',0);
        ends = cellfun(@(x) x(end), ts);
        offset = [0 ends(1:end-1)+ts{1}(1)];
        offset = cumsum(offset);
        ts=arrayfun(@(x,y) x{1}+y, ts,offset,'UniformOutput',0);
        ts=CMBHOME.Utils.ContinuizeEpochs(ts(:));
        
        %spike
        spk_ts = cell(size(spike{1}));
        for i = 1:numel(spike)
            for k = 1:numel(spike{i})
                spk_ts{k} = [spk_ts{k}; spike{i}(k).ts + offset(i)-firsts(i)];
            end
        end
        
        for i = 1:numel(spk_ts)
            if ~isempty(spk_ts{i})
                spike2(i) = CMBHOME.Spike('ts',spk_ts{i},'vid_ts',ts);
            else
                spike2(i) = CMBHOME.Spike();
            end
        end
        sz=max(CMBHOME.Utils.ContinuizeEpochs(cellfun(@(x) size(x), spike(:),'UniformOutput',0)));
        spike=reshape(spike2,sz);
        
        %LFP
        lfp_fs=lfp{1}.fs; 
        lfp_sig=cell(size(lfp{1})); lfp_ts=cell(size(lfp{1}));
        for i = 1:numel(lfp)
            for k = 1:numel(lfp{i})
                lfp_ts{k} = [lfp_ts{k}; lfp{i}(k).ts+offset(i)-firsts(i)];
                lfp_sig{k} = [lfp_sig{k}; lfp{i}(k).signal];
            end
        end
        
        for i = 1:numel(lfp_ts)
            lfp2(i) = CMBHOME.LFP(lfp_sig{i}, lfp_ts{i},lfp_fs,'merged');
        end
        
        % Make it:
        root = Session('b_x',x,'b_y',y,'b_headdir',hd,'b_ts',ts,...
                       'fs_video',fs, 'spatial_scale',ss,...
                        'b_lfp',lfp2, 'spike', spike);
        
        end
    end
    
    methods
                
        function self = Session(varargin)
            
            import CMBHOME.Spike
            
            p = inputParser;

            p.addParamValue('name',          'default session name', @(x) ischar(x)||iscellstr(x));
            p.addParamValue('b_x',          [], @(x) any(size(  x)<=1)); 
            p.addParamValue('b_y',          [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_headdir',        [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_ts',         [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_lfp',         CMBHOME.LFP, @(x) isstruct(x)||isa(x,'CMBHOME.LFP'));
            p.addParamValue('fs_video',     [], @(x) length(x)==1);
            p.addParamValue('raw_pos',      [], @(x) length(x)==1);
            p.addParamValue('raw_headdir',      [], @(x) length(x)==1);
            p.addParamValue('event',       cell(1,2), @(x) iscell(x));
            p.addParamValue('rem_ts',       [], @(x) length(x)==1);
            p.addParamValue('ripple_ts',    [], @(x) length(x)==1);
            p.addParamValue('path_lfp',     []);
            p.addParamValue('path_raw_data',[]);
            p.addParamValue('date_created', now, @(x) isnumeric(x));
            p.addParamValue('spatial_scale',.5, @(x) length(x)==1);
            p.addParamValue('epoch',        [-inf inf], @(x) numel(x)==2);
            p.addParamValue('spike',        Spike, @(x) isa(x, 'CMBHOME.Spike'));
                 
            p.parse(varargin{:});
            
            self.name = p.Results.name;
            self.notes = { date, ['Session object created, and saved to ', p.Results.name] };
            self.b_x = p.Results.b_x;
            self.b_y = p.Results.b_y;
            self.b_headdir = p.Results.b_headdir;
            self.b_ts = p.Results.b_ts;
            self.b_lfp = p.Results.b_lfp;
            self.fs_video = p.Results.fs_video;
            %self.fs = p.Results.fs;
            self.raw_pos = p.Results.raw_pos;
            self.raw_headdir = p.Results.raw_headdir;
            self.event = p.Results.event;
            self.path_lfp = p.Results.path_lfp;
            self.path_raw_data = p.Results.path_raw_data;
            self.date_created = p.Results.date_created;
            self.spatial_scale = p.Results.spatial_scale;
            
            if ~isempty(self.b_ts)
                self.epoch = p.Results.epoch;
            end
            
            if isempty(self.fs_video)&&~isempty(self.b_ts)
                self.fs_video = 1/mean(diff(self.b_ts));
            end
            
            if isempty(self.b_ts)&&~isempty(self.fs_video)
                self.b_ts = (0:numel(self.b_x)-1)/self.fs_video;
            end
                
            
            self.spike = p.Results.spike;
        end   
        
        
        %% DEPENDENT STATE VARIABLES
        
        function self = set.epoch(self, epoch)
        % makes sure that all epochs are valid periods of time, and that
        % the array is formatted correctly. 
               
            if isempty(self.b_ts), error('root.b_ts must be set to use epochs'); end
        
            %tmin = self.b_ts(1);
            %tmax = self.b_ts(end);
            
            if numel(epoch)==2, epoch = epoch(:)'; end

            if size(epoch,2)~=2
                error('Epoch must be a matrix with columns [start times, stop times]');
            end

            if any(epoch(:,2)-epoch(:,1)<=0) % check that all starts and stops are start<stop
                error('epoch must be in format [tstart, tstop]');
            end
            
            epoch(epoch==inf) = self.b_ts(end);
            epoch(epoch==-inf) = self.b_ts(1);

            if isempty(self.b_vel)
                %self = AppendKalmanVel(self);
                self.b_vel = sqrt(([0;diff(self.b_x)]).^2+([0;diff(self.b_y)]).^2) * self.fs_video;
                warning('No vel set, setting simple derivative. Recommended: root=root.AppendKalmanVel   (slow)');
            end
            
            if all(size(epoch)==size(self.p_epoch))
                
                if all(epoch==self.p_epoch) % you are reassigning the same thing!! dont run all the scripts below
                    
                    return
                    
                end
                
            end

            self.p_epoch = epoch;
            
            self.p_ind = IsolateEpoch(self.b_ts, (1:length(self.b_ts))', self.epoch);
            
            self = SetLFPInds(self); % update the aligned lfp inds
            
            self = SetCelInds(self); % update the indices for the cells which are set
            
            if ~isempty(self.myvar2)
                self.myvar2 = self.myvar2.setEpoch(self.ind,self.b_ts);
            end
            
        end  
        
        function epoch = get.epoch(self)
            
            epoch = self.p_epoch;
            
        end
        
        function self = set.active_lfp(self, active_lfp)
            
            if active_lfp==self.p_active_lfp; return; end % its already set, don't rerun the indexing
            
            if numel(active_lfp) == 1 && active_lfp>0 && active_lfp<=length(self.b_lfp)
                self.p_active_lfp = active_lfp;
            elseif isempty(active_lfp)
                self.p_active_lfp = [];
            else
                disp('root.active_lfp not set because it was not a valid index for root.b_lfp')
                self.p_active_lfp = [];
            end
            
            self = SetCelInds(self);
            
        end
        
        function active_lfp = get.active_lfp(self)
            
            active_lfp = self.p_active_lfp;
            
        end
        
        function self = set.cel(self, cel)
        % sets the root.cel property (Nx2, [tetrode, cell_ind])

            cells = ValidCells(self);
            %cells = cel;
            if isempty(cel), self.p_cel = []; self = SetCelInds(self); return, end
        
            if size(cel, 2)~=2, warning('CMBH:error', 'root.cel must be Nx2, [tetrode, cell_ind]'); return; end
        
            if ~isnumeric(cel), warning('CMBH:error', 'root.cel must be Nx2, [tetrode, cell_ind]'); return; end
            
            if any(~ismember(cel, cells, 'rows'))
                warning('CMBH:error', 'error, some of root.cel are not cells, or do not pass root.cell_thresh. removed the following cells:');
                disp(num2str(cel(find(~ismember(cel, cells, 'rows')),:)));
                cel(find(~ismember(cel, cells, 'rows')),:) = [];
            end
            
            if all(size(cel)==size(self.p_cel))
                
                if all(cel==self.p_cel) % you're trying to reset the same cells, so dont rerun the indexing
                    
                    return
                    
                end
                
            end
            
            self.p_cel = cel;
            
            self = SetCelInds(self); % update the indices for the cells which are set
            
        end
        
        function cel = get.cel(self)
            
            cel = self.p_cel;
            
        end
        
        
        %% DEPENDENT BASE VARIABLES
        
        function self = set.b_x(self, b_x)
        % sets p_b_x. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_x) && ~isvector(b_x)) && ~isempty(b_x)
                error('X data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_x)~=l) && ~isempty(b_x)
                
                error('b_x must be the same length as all other base variables');
                
            end
            
            self.p_b_x = b_x(:);
            
        end
        
        function self = set.b_y(self, b_y)
        % sets p_b_y. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_y) && ~isvector(b_y)) && ~isempty(b_y)
                error('y data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_y)~=l) && ~isempty(b_y)
                
                error('b_y must be the same length as all other base variables');
                
            end
            
            self.p_b_y = b_y(:);
            
        end
        
        function self = set.b_vel(self, b_vel)
        % sets p_b_vel. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_vel) && ~isvector(b_vel)) && ~isempty(b_vel)
                error('vel data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_vel)~=l) && ~isempty(b_vel)
                
                error('b_vel must be the same length as all other base variables');
                
            end
            
            self.p_b_vel = b_vel(:);
            
        end
        
        function self = set.b_headdir(self, b_headdir)
        % sets p_b_headdir. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_headdir) && ~isvector(b_headdir)) && ~isempty(b_headdir)
                error('headdir data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_headdir)~=l) && ~isempty(b_headdir)
                
                error('b_headdir must be the same length as all other base variables');
                
            end
            
            self.p_b_headdir = b_headdir(:);
            
        end
        
        function self = set.b_ts(self, b_ts)
        % sets p_b_ts. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_ts) && ~isvector(b_ts)) && ~isempty(b_ts)
                error('ts data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_ts)~=l) && ~isempty(b_ts)
                
                error('b_ts must be the same length as all other base variables');
                
            end
            
            self.p_b_ts = b_ts(:);
            
        end
        
        function self = set.b_myvar(self, b_myvar)
        % sets p_b_myvar. must be either empty or the same length as other base variables
        
                if (~isnumeric(b_myvar) && ~isvector(b_myvar)) && ~isempty(b_myvar)
                    error('y data must be numerical vector');
                end

                l = CheckBaseVarLength(self);

                if (l~=0 && length(b_myvar)~=l) && ~isempty(b_myvar)

                    error('b_myvar must be the same length as all other base variables');

                end

                self.p_b_myvar = b_myvar(:);
            
        end
        
        
        function self = set.myvar2(self,varargin2)
            % Because we can't have multiple constructors in matlab use the
            % varargin2 variable. Possible ways of setting root.myvar2 are:
            %
            % 1) root.myvar2 = []; clear root.myvar2
            %
            % 2) temp = CMBHOME.Utils.MyvarHolder('x',x,t); 
            %    root.myvar2 = temp; % Use a temporary variable
            %
            % 3) root.myvar2 = {'x',x,t}; % Uses a cell array to pass the necessary information
            %
            % Possibly rewrite in the future for clarity.
            
            if ~exist('varargin2','var'),varargin2 = [];end
            
            if iscell(varargin2)
               name = varargin2{1}; %#ok<PROP>
               p_var = varargin2{2};
               p_ts = varargin2{3};
               tsVar = self.b_ts; %#ok<MCSUP>
               self.myvar2 = CMBHOME.Myvar(name,p_var,p_ts,tsVar); %#ok<PROP> 
               
            elseif (isa(varargin2,'CMBHOME.Myvar'))
                self.myvar2 = varargin2;
                
            elseif (isa(varargin2,'CMBHOME.Utils.MyvarHolder'))
                self.myvar2 = CMBHOME.Myvar(varargin2);   
                
            elseif (isempty(varargin2))
                self.myvar2 = CMBHOME.Myvar(); %clear the root.myvar2
            else
                 error('mv2 must be a CMBHOME.myvar or the variables required to construct one')
            end
        end
        
        
        function b_x = get.b_x(self)
            
            b_x = self.p_b_x;
            
        end
        
        function b_y = get.b_y(self)
            
            b_y = self.p_b_y;
            
        end
        
        function b_ts = get.b_ts(self)
            
            b_ts = self.p_b_ts;
            
        end
        
        function b_vel = get.b_vel(self)
            
            b_vel = self.p_b_vel;
            
        end
        
        function b_headdir = get.b_headdir(self)
            
            b_headdir = self.p_b_headdir;
            
        end
        
        function b_myvar = get.b_myvar(self)
            
            b_myvar = self.p_b_myvar;
            
        end
        
        
        %% DEPENDENT NON-BASE VARIABLES
        
        function cells = get.cells(self)
            cells = ValidCells(self);
        end
        
        function name_formatted = get.name_formatted(self)
        % returns the filename of the current object formatted for the current platform
        
            if isempty(self.name)

                name_formatted = 'CMBHno_name.mat';

            elseif iscellstr(self.name)

                name_formatted = fullfile(self.name{:});

            else

                name_formatted = self.name;

            end
            
        end
        
        function active_event = get.active_event(self)
        % e = root.active_event;
        %
        % Returns strings which correspond to the start and stop labels for
        % each epoch time.
        %
        % RETURNS
        %
        % e -> Nx2 cell array of strings (same size as root.epoch, or if
        % root.epoch_group is set, N is the number of groups, and it 
        % defaults to the group's first label)
        %
        % andrew 18 june 2010
        
            if size(self.event,2)<2
                
                active_event = {'No Events', 'No Events'};
            
            else
            
                active_event = cell(size(self.epoch,1), 2);
                counter = 1;
                for i = self.b_epoch_group'
                    
                    inds = find(self.b_epoch_group==i);
                    
                    if sum(vertcat(self.event{:,2})==self.epoch(inds(1), 1))==1
                        str1 = self.event{vertcat(self.event{:,2})==self.epoch(inds(1), 1),1};
                    else
                        str1 = ['Unknown t = ' num2str(self.epoch(inds(1), 1))];
                    end

                    if sum(vertcat(self.event{:,2})==self.epoch(inds(1), 2))==1
                        str2 = self.event{vertcat(self.event{:,2})==self.epoch(inds(1), 2),1};
                    else
                        str2 = ['Unknown t = ' num2str(self.epoch(inds(1), 2))];
                    end

                    active_event(counter, 1) = {str1};
                    active_event(counter, 2) = {str2};
                    
                    counter = counter+1;
                end
            
            end
            
        end
        
        function vel = get.vel(self)
        % returns pixels/sec
            if isempty(self.b_vel)
                %self = AppendKalmanVel(self);
                self.b_vel = sqrt(([0;diff(self.b_x)]).^2+([0;diff(self.b_y)]).^2) * self.fs_video;
                warning('No vel set, setting simple derivative. Recommended: root=root.AppendKalmanVel   (slow)');
            end
            
            vel = [];

            if iscell(self.p_ind)
                vel = cellfun(@(c) self.b_vel(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_vel)
                vel=self.b_vel(self.p_ind);
            end
        end
                
        function ind = get.ind(self)
        % These are indices in the base vectors (b_*) used to return data dynamic fields (ex. root.x = root.b_x(root.ind))
        %
        % This is useful if you have a user defined vector that isnt
        % included in the set of properties already defined. For example, a
        % linearized position in the TMaze could be saved in
        % root.user_def.linear_t. Then for any epoch the linear position on
        % the track is root.user_def.linear_t(root.ind), so long as
        % length(linear_t)==length(root.b_ts)        
        
            ind = self.p_ind;
                        
        end
        
        function x = get.x(self)
        % x position in pixels
            
            x = [];
                        
            if iscell(self.p_ind)
                x = cellfun(@(c) self.b_x(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_x)
                x=self.b_x(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No x data for this session.');
            end
            
        end
        
        function sx = get.sx(self)
           pos = self.spos;
           if iscell(pos)
               sx = cellfun(@(x)x(:,1),pos,'UniformOutput',false);
               sx = sx(:);
           else
               sx = pos(:,1);
           end
        end
        
        function svel = get.svel(self)
           svel = self.vel;
           if iscell(svel)
              svel = cellfun(@(x)x*self.spatial_scale,svel,'UniformOutput',false); %eln 131002
              svel = svel(:);
           else
              svel = svel*self.spatial_scale; 
           end
        end
        
        function sy = get.sy(self)
           pos = self.spos;
           if iscell(pos)
               sy = cellfun(@(x)x(:,2),pos,'UniformOutput',false);
               sy = sy(:);
           else
               sy = pos(:,2);
           end
        end
        
        function y = get.y(self)
        % y position in pixels

            y = [];
            
            if iscell(self.p_ind)
                y = cellfun(@(c) self.b_y(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_y)
                y=self.b_y(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No y data for this session');
            end
            
        end
        
        function headdir = get.headdir(self)
        % headdirection in degrees (although the user could input radians, 
        % all Toolbox functionality assumes degrees)  

            headdir = [];
                     
                if iscell(self.p_ind)
                    headdir = cellfun(@(c) self.b_headdir(c), self.p_ind, 'UniformOutput', false);
                elseif ~isempty(self.b_headdir)
                    headdir=self.b_headdir(self.p_ind);
                else
                    warning('CMBH:notify', '%s', 'Head Direction vector does not exist for this session');
                end
        end
        
        function sheaddir = get.sheaddir(self)
            sheaddir = self.headdir;
            
            if iscell(sheaddir)
               sheaddir = cellfun(@(x)mod(x+self.rotate+180,360)-180,sheaddir,'UniformOutput',false); 
            else
                sheaddir = mod(sheaddir+self.rotate+180,360)-180;
            end
        end
        
        function ts = get.ts(self)
        % timestamps in seconds    
            
            if iscell(self.p_ind)
                ts = cellfun(@(c) self.b_ts(c), self.p_ind, 'UniformOutput', false);
            else
                ts=self.b_ts(self.p_ind);
            end
            
        end
        
        function myvar = get.myvar(self)
        % epoch data from b_myvar

            if isempty(self.b_myvar)
                myvar = [];
                return
            elseif length(self.b_myvar)~=length(self.b_ts) && ~isempty(self.b_myvar)
                disp('root.b_myvar must be the same length as root.b_ts');
                myvar = [];
                return

            end

            myvar = [];

            if iscell(self.p_ind)
                myvar = cellfun(@(c) self.b_myvar(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_myvar)
                myvar=self.b_myvar(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No myvar data for this session');
            end

        end
        
        function lfp = get.lfp(self)
            % check to see if b_lfp is populated. if not, check that there
            % is a file by the index at root.active_lfp
            
            % if there is, load lfp it
            
            % then call self.lfp.signal and self.lfp.time and get chucks as
            % per self.epoch
            
            if exist('self.b_lfp') % check if it's loaded yet
                if length(self.b_lfp)<self.active_lfp
                    
                    error('root.active_lfp exceeds length of root.b_lfp');
                    
                end
            elseif ~isempty(self.active_lfp)
                
                signal = self.b_lfp(self.active_lfp).signal;
                ts = self.b_lfp(self.active_lfp).ts;
                              
                b_theta = self.b_lfp(self.active_lfp).b_theta;

                b_theta_phase = self.b_lfp(self.active_lfp).b_theta_phase;

                b_theta_amplitude = self.b_lfp(self.active_lfp).b_theta_amplitude;
                
                b_myvar = self.b_lfp(self.active_lfp).b_myvar;
                
                if isempty(b_theta), b_theta = nan(length(signal), 1); end
                
                if isempty(b_theta_phase), b_theta_phase = nan(length(signal), 1); end
                
                if isempty(b_theta_amplitude), b_theta_amplitude = nan(length(signal), 1); end
            
            else
                
                lfp = [];
                return
            
            end
            %% grab data from epochs
            
            if size(self.epoch,1)>1
                
                tmp_signal = cell(size(self.epoch,1), 1);
                tmp_ts = cell(size(self.epoch,1), 1);
                tmp_b_theta = cell(size(self.epoch,1), 1);
                tmp_b_theta_phase = cell(size(self.epoch,1), 1);
                tmp_b_theta_amplitude = cell(size(self.epoch,1), 1);
            
                for i = 1:size(self.epoch,1)
                    
                    inds = ts>=self.epoch(i,1) & ts<=self.epoch(i,2);
                                        
                    tmp_signal{i} = signal(inds);
                    tmp_ts{i} = ts(inds);
                    tmp_b_theta{i} = b_theta(inds);
                    tmp_b_theta_phase{i} = b_theta_phase(inds);
                    tmp_b_theta_amplitude{i} = b_theta_amplitude(inds);
                    
                    if ~isempty(b_myvar)
                        tmp_b_myvar{i} = b_myvar(inds,:); %#ok<AGROW>
                    end
                    
                end
                
                signal = tmp_signal;
                ts = tmp_ts;
                b_theta = tmp_b_theta;
                b_theta_phase = tmp_b_theta_phase;
                b_theta_amplitude = tmp_b_theta_amplitude;
                
                if exist('tmp_b_myvar','var')
                    b_myvar = tmp_b_myvar(:);
                end
            else
                
                inds = ts>=self.epoch(1) & ts<=self.epoch(2);
                
                signal = signal(inds);
                ts = ts(inds);
                b_theta = b_theta(inds);
                b_theta_phase = b_theta_phase(inds);
                b_theta_amplitude = b_theta_amplitude(inds);
                if ~isempty(b_myvar)
                    b_myvar = b_myvar(inds,:);
                else
                    b_myvar = [];
                end
                
            end
            
            user_def = self.b_lfp(self.active_lfp).user_def; %#ok<PROP>
            
            fs = self.b_lfp(self.active_lfp).fs;
            
            channel_name = self.b_lfp(self.active_lfp).channel_name;
            
            lfp = CMBHOME.LFP(signal, ts, fs, channel_name, b_theta, b_theta_phase, b_theta_amplitude, user_def, b_myvar); %#ok<PROP>
            
        end
        
        function b_epoch_group = get.b_epoch_group(self)
        % Checks to see if root.epoch_group is set by user. If it is, and 
        % is valid (length=number of epochs), then root.b_epoch_groups ==
        % root.epoch_groups
        
            if length(self.epoch_group) == size(self.epoch,1) && all(self.epoch_group>0)

                b_epoch_group = self.epoch_group;

            else

                b_epoch_group = 1:size(self.epoch,1);

                b_epoch_group = b_epoch_group';

            end            
            
        end
        
        function cel_x = get.cel_x(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_x = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_x = cellfun(@(c) self.b_x(c), self.p_cel_ind, 'unif', 0);
            
        end
                        
        function cel_y = get.cel_y(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_y = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_y = cellfun(@(c) self.b_y(c), self.p_cel_ind, 'unif', 0);
            
        end

        function cel_sx = get.cel_sx(self)
            
            if isempty(self.p_cel_ind)
                disp('Set root.cel.');
                cel_sx = cell(size(self.epoch,1),1);
                return
            end
            
            cel_x = self.cel_x;
            cel_y = self.cel_y;
            cel_sx = self.spos(cel_x,cel_y);
            
            if ~isempty(self.cel)
               cel_x = self.cel_x;
               cel_y = self.cel_y;

               if iscell(cel_sx)
                  cel_sx = cellfun(@(x)x(:,1),cel_sx,'UniformOutput',false); 
               else
                  cel_sx = cel_sx(:,1);
               end 
            end
        end
        
        function cel_sy = get.cel_sy(self)
            if isempty(self.p_cel_ind)
                disp('Set root.cel.')
                cel_sy = cell(size(self.epoch,1),1);
                return
            end
            
            cel_x = self.cel_x;
            cel_y = self.cel_y;
            
            cel_sy = self.spos(cel_x,cel_y);
            
            if iscell(cel_sy)
              cel_sy = cellfun(@(x)x(:,2),cel_sy,'UniformOutput',false); 
            else
              cel_sy = cel_sy(:,2);
            end
            
        end
        
        function cel_vel = get.cel_vel(self)
            
            if isempty(self.p_cel_ind)         
                disp('Set root.cel.'); 
                cel_vel = cell(size(self.epoch,1), 1);
                return
            end
            
            if ~isempty(self.b_vel),
                cel_vel = cellfun(@(c) self.b_vel(c), self.p_cel_ind, 'unif', 0);
            else
                vel = cat(1, NaN, sqrt(diff(self.b_x).^2+diff(self.b_y).^2))*self.fs_video;
                cel_vel = cellfun(@(c) vel(c), self.p_cel_ind, 'unif', 0);  
            end

        end
        
        function cel_svel = get.cel_svel(self)
            if isempty(self.p_cel_ind)
                disp('Set root.cel.');
                cel_svel = cell(size(self.epoch,1),1);
                return
            end
            
           cel_svel = self.cel_vel;
           
           if iscell(cel_svel)
              cel_svel = cellfun(@(x) x*self.spatial_scale,cel_svel,'UniformOutput',false); 
           else
              cel_svel = cel_svel*self.spatial_scale;
           end
        end
        
        function cel_ts = get.cel_ts(self)
            %if isempty(self.p_cel_ind)
            %    disp('Set root.cel.'); 
            %    cel_ts = cell(size(self.epoch,1), 1);
            %    return
            %end
            
            if isempty(self.cel)
                disp('set root.cel.')
                cel_ts = cell(size(self.epoch,1),1);
                return
            end
            
            cel_ts = cell(size(self.epoch,1),size(self.cel,1));
            
            for i = 1:size(self.cel,1)
                cel = self.cel(i,:);
                spks_ts = self.spike(cel(1), cel(2)).ts;
                
                for k = 1:size(self.epoch,1)
                    curInds = spks_ts >= self.epoch(k,1) & spks_ts <= self.epoch(k,2);
                    cel_ts{k,i} = spks_ts(curInds);
                end
                %cel = self.cel(i,:);
                %cel_ts(:,i) = cellfun(@(c) self.spike(cel(1), cel(2)).ts(c), self.p_cel_spkind(:,i), 'unif', 0); %#ok<AGROW>
            end  
        end

        function cel_headdir = get.cel_headdir(self)
            
            if isempty(self.p_cel_ind)
                disp('Set root.cel.'); 
                cel_headdir = cell(size(self.epoch,1), 1);
                return
            end
            cel_headdir = cellfun(@(c) self.b_headdir(c), self.p_cel_ind, 'unif', 0);
        end

        function cel_sheaddir = get.cel_sheaddir(self)
           
            if isempty(self.p_cel_ind)
                disp('Set root.cel.');
                cel_sheaddir = cell(size(self.epoch,1),1);
                return
            end
            
            cel_sheaddir=self.cel_headdir;
            
            if iscell(cel_sheaddir)
               cel_sheaddir = cellfun(@(x)mod(x+self.rotate+180,360)-180,cel_sheaddir,'UniformOutput',false); 
            else
               cel_sheaddir = mod(cel_sheaddir+self.rotate+180,360)-180;
            end
            
        end
        
        function cel_i = get.cel_i(self)
            if isempty(self.p_cel_ind)
                disp('Set root.cel.'); 
                cel_i = self.p_cel_ind;
                return
            end
            
            cel_i = self.p_cel_ind;
        end

        function cel_myvar = get.cel_myvar(self)
                            
            if isempty(self.p_cel_ind)    
                disp('Set root.cel.'); 
                cel_myvar = cell(size(self.epoch,1), 1);  
                return  
            end
            
            % eln 131002
            if isempty(self.b_myvar)
              cel_myvar = cell(size(self.epoch,1), 1); 
              return             
            end
            
            cel_myvar = cellfun(@(c) self.b_myvar(c), self.p_cel_ind, 'unif', 0);
        end
        
        function cel_myvar2 = get.cel_myvar2(self)
            % cel_myvar_xvar(:) = root.cel_myvar2.get('xvar');
            %
            % Then cel_myvar_xvar will be cell array of "xvar" at spike
            % times.
            %
            % root.cel_myvar2 takes the form of a Myvar, but only holds
            % spike data. 
            
            % wchapman 22/07/2013
            if isempty(self.myvar2)
                cel_myvar2 = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.p_cel_ind)
                cel_myvar2 = cell(size(self.epoch,1), 1);
                return
            end
            
            cel_myvar2 = self.myvar2.setEpoch(self.p_cel_ind,self.b_ts);
            
        end
        
            
        function cel_theta = get.cel_theta(self)
           
            if isempty(self.cel) % if no cel is set, return empty
                disp('Set root.cel');
                cel_theta = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.active_lfp) % if no active lfp set 
                disp('Set root.active_lfp');
                cel_theta = cell(size(self.epoch,1), size(self.cel, 1));
                return                
            end
            
            cel_theta = cellfun(@(c) self.b_lfp(self.active_lfp).theta(c), self.p_cel_lfp_ind, 'unif', 0);
                
        end
        
        function cel_thetaphase = get.cel_thetaphase(self)
           
            if isempty(self.cel) % if no cel is set, return empty
                disp('Set root.cel');
                cel_thetaphase = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.active_lfp) % if no active lfp set 
                disp('Set root.active_lfp');
                cel_thetaphase = cell(size(self.epoch,1), size(self.cel, 1));
                return                
            end
            
            cel_thetaphase = cellfun(@(c) self.b_lfp(self.active_lfp).theta_phase(c), self.p_cel_lfp_ind, 'unif', 0);
                
        end
        
        function cel_lfpmyvar = get.cel_lfpmyvar(self)
           
            if isempty(self.cel) % if no cel is set, return empty
                disp('Set root.cel');
                cel_lfpmyvar = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.active_lfp) % if no active lfp set 
                disp('Set root.active_lfp');
                cel_theta = cell(size(self.epoch,1), size(self.cel, 1));
                return                
            end
            
            cel_lfpmyvar = cellfun(@(c) self.b_lfp(self.active_lfp).myvar(c), self.p_cel_lfp_ind, 'unif', 0);
                
        end
                
        
        %% INDEPENDENT PROPERTIES
        
        function self = set.cell_thresh(self, cell_thresh)
            
            if isnumeric(cell_thresh)
                
                if numel(cell_thresh)==1
                    cell_thresh = [cell_thresh(1), 0];
                elseif numel(cell_thresh)==2
                    if cell_thresh(2)~=1 && cell_thresh(2)~=0
                        warning('CMBH:error', 'cell_thresh(2) must be boolean: average Fmin (0) or sustained Fmin (1)? Defaulting to 0.')
                        cell_thresh(2) = 0;
                    end
                else
                    warning('CMBH:error', 'Improper number of elements on cell_thresh');
                end
                
                switch cell_thresh(2)
                case 1
                    str_cell = ' Hz continuously throughout the entire session';
                case 0
                    str_cell = ' Hz on average during the session';
                end
                
                self.cell_thresh = cell_thresh;

                % warning('CMBH:notify', '%s', ['Cell Threshold: F_min is ' num2str(self.cell_thresh(1)) str_cell '.']);
                
            end
        end
        
        function self = set.epoch_group(self, epoch_group)
           
            if all(epoch_group>0) && all(rem(epoch_group, 1)==0) && isvector(epoch_group)
               
                self.epoch_group = epoch_group(:);
                
            end
            
        end
                        
        function self = set.event(self, event)
            if ~iscell(event)
                error('session.event must be cell array')
            else
                self.event=event;
            end
            
            if size(event, 2)==2 && size(event,1)>1 % sort events if numeric
            
                isnumbers = cellfun(@isnumeric, event(:,2));
                isemptys = cellfun(@isempty, event(:,2));
                
                if any(~isnumbers) || any(isemptys)
                    
                    warning('CMBH:notify', '%s', 'There are values in the timestamps column (2) that are not numbers. They will be deleted');
                    
                    event(~isnumbers,:) = [];
                    event(isemptys,:) = [];
                    
                end
                
                ts = [event{:,2}];
                
                if isnumeric(ts)
                   
                    [~, ind] = sort(ts); 
                    
                    self.event = event(ind,:);
                    
                end
                
            elseif isempty(event)
                
                self.event = {'Session start', self.b_ts(1); 'Session stop', self.b_ts(end) }; %#ok<MCSUP>
                
            end
            
        end
        
        function self = set.path_lfp(self, path_lfp)
            
            if isempty(path_lfp)
                self.path_lfp = [];
                return;
            end
            
            if ~iscell(path_lfp), path_lfp = {path_lfp}; end
            self.path_lfp = path_lfp;
        end
        
       
        function self = set.b_lfp(self, b_lfp)
           
            issue_warning = 0;
            
            if isa(b_lfp, 'CMBHOME.LFP')
            
                for i = 1:length(b_lfp)
                    
                    %if any(unique(diff(b_lfp(i).ts)).^-1 < .9*b_lfp(i).fs), issue_warning = 1; end
                    %!IMPLEMENTME: This line (1469) is broken in recent
                    %MATLAB upgrade. Need to fix definition of empty

                end
                
            end
            
            if issue_warning, warning('CMBH:notify', '%s', 'There are timestamps with uneven sampling intervals. Make sure this is expected for this session, and to use root.epoch properly.'); end
        
            self.b_lfp = b_lfp; 
            
        end
        
        
        %% METHODS
        
        function pos=spos(self,x,y)
            if ~exist('x','var')
                x = self.x;
                y = self.y;
            end
            
            if iscell(x)
                % Switched from cellfun to self call for each cell  wchapman 20131003
                pos = cell(size(x));
                for i = 1:numel(x)
                    pos{i} = self.spos(x{i},y{i});
                end
            else
                pos = self.spatial_scale*[x y]*[cosd(self.rotate) -sind(self.rotate);sind(self.rotate) cosd(self.rotate)]+repmat(self.shift,[numel(x) 1]);
                tolerance = 0.01;
                
                if size(pos,2)>1
                    count = 0;
                    while abs((sum(sqrt(sum(diff([x y]).^2,2)))*self.spatial_scale-sum(sqrt(sum(diff(pos).^2,2))))/(sum(sqrt(sum(diff([x y]).^2,2)))*self.spatial_scale))>tolerance
                        count = count+1;
                        pos = self.spatial_scale*[x-self.shift(1)/self.spatial_scale y-self.shift(2)/self.spatial_scale]*[cosd(self.rotate) -sind(self.rotate);sind(self.rotate) cosd(self.rotate)];
                        if count > size(pos,2)
                            warning('Session:spos fail','Max iterations reached in spos')
                            break
                        end
                    end
                end
            end
        end
       
        function self = AppendKalmanVel(self)
            
            if ~isempty(self.b_vel)
                
                warning('CMBH:error', 'User defined speed vector already created. Clear root.b_vel first.');
                
                return
                
            else
            
                [t,~,~,vx,vy] = self.KalmanVel(self.b_x,self.b_y,self.b_ts, 2); % quadratic

                [~, inds] = ismember(t, self.b_ts); % all indices in self.b_ts for which we have a kalman estimate

                b_vel = zeros(size(self.b_ts));

                b_vel(inds) = sqrt(vx.^2 + vy.^2); % pixels/sec
                
                self.b_vel = b_vel;
                
                %warning('CMBH:notify', '%s', 'Added Kalman Velocity to Session Object.');
            
            end
        end
        
        function ts = Label2Time(self, label)
        % ts = root.Label2Time(str)
        %
        % Returns ts, a vector of timestamp(s) from root.event that line up
        % with event label 'str'
        
        
        if iscell(label)
            
            ts = zeros(size(label, 1), size(label,2));
            
            for i = 1:size(label, 1)
                for j = 1:size(label, 2)
                    tmp = cell2mat(self.event(ismember(self.event(:,1), label{i,j}),2));
                    
                    if numel(tmp)==1, ts(i,j) = tmp;
                    elseif isempty(tmp), ts(i, j) = NaN;
                    else error('Multiple identical labels exist, cannot use cell array functionality');
                    end
                    
                end
            end
        else
        
            ts = cell2mat(self.event(ismember(self.event(:,1), label),2));
        
        end

        end
        
        function self = AddLFP(self, fname)
        % (1) root = root.AddLFP
        % (2) root = root.AddLFP(fname)
        %
        % Appends root.path_lfp with valid LFP file name via user input
        % prompts or filename, and adds it to root.path_lfp.
        %
        % (1) Prompts the user to select a file using the file browser
        % (2) Checks if string 'fname' is a valid LFP file (.ncs, .mat,
        % .plx files) and adds it to root.path_lfp.
        %
        % andrew 3 april 2010
                        
            if ~exist('fname', 'var')
                if isempty(self.path_lfp) 
                    i = 1;
                else
                    i = length(self.path_lfp)+1;
                end

                tmppath = self.name;

                tmppath = tmppath(1:find(tmppath==filesep, 1, 'last')-1); % directory of datafile

                [load_file,base_path] = uigetfile('*.eeg*;*.ncs;*.mat;*.plx','Select the LFP (EEG) file relevent to this project. Cancel to exit',tmppath, 'MultiSelect', 'on');

                if iscellstr(load_file)
                    
                    for j = 1:length(load_file)
                        self.path_lfp{i} = fullfile(base_path, load_file{j});
                        i = i+1;
                    end

                    self = self.AddLFP; % recursively call this function to load files
                    
                    warning('CMBH:notify', '%s', 'Added'), warning('CMBH:notify', '%s', load_file(:));
                    
                elseif load_file
                                       
                    self.path_lfp{i} = fullfile(base_path, load_file);

                    self = self.AddLFP; % recursively call this function to load files
                    
                    warning('CMBH:notify', '%s', 'Added'), warning('CMBH:notify', '%s', load_file(:));
                    
                end
            else
                
                if ~iscell(fname), fname = {fname}; end % convert to cell if char
                
                for i = 1:length(fname)
                    
                    if exist(fname{i}, 'file')
                        self.path_lfp{end+1} = fname{i};
                    else
                        warning('CMBH:error', [fname{i} ' does not exist. Not added']);
                        fname{i} = '';
                    end
                    
                end
                
                warning('CMBH:notify', '%s', ['Added ' int2str(length(fname)) ' lfp filenames.']);
                
            end
                                       
        end  
        
        function self = ClearLFP(self)
        % root = root.ClearLFP;
        %
        % Clears the LFP signals from root.b_lfp. This may be useful if the
        % LFP signals are very large, or highly sampled.
        %
        % andrew 3 april 2010
        
            self.b_lfp = [];
            
            self.active_lfp = [];
                         
        end 
        
        function Status(self)
        % root.Status
        % prints to screen relevant information about this Session class
        % object
        %
        % andrew 3 april 2010

            switch self.cell_thresh(2)
                case 1
                    str_cell = ' Hz continuously throughout the entire session';
                case 0
                    str_cell = ' Hz on average during the session';
            end

            disp(['This CMB Session is ' self.name '.' char(10) 'It has ' int2str(size(self.cells,1)) ...
                ' cells which fire above ' num2str(self.cell_thresh(1)) ...
                str_cell '.' char(10) char(10) 'There are ' int2str(size(self.event,1)) ' event flags' ...
                '.' char(10) 'The recording duration is ' num2str(self.b_ts(end)-self.b_ts(1)) ' seconds.']);
            
            disp([char(10) 'It was created on ' datestr(self.date_created, 1)]);
            
            disp([char(10) 'The tracking sampling frequency is ' num2str(self.fs_video) ' Hz.']);
            
            if isempty(self.b_vel), disp('root.vel returns the simple derivative of position. Edit root.b_vel to specify a more accurate vector of running speed,');
            else disp('User defined speed vector'); end
            
            switch self.raw_headdir
                case 1
                    str_headdir = 'The head direction signal has not been corrected or smoothed';
                case 0
                    str_headdir = 'The head direction signal has been corrected and smoothed';
            end
            
            switch self.raw_pos
                case 1
                    str_pos = 'The position data has not been corrected or smoothed.';
                case 0
                    str_pos = 'The position data has been corrected and smoothed.';
            end
            
            switch isempty(self.spatial_scale)
                case 1
                    str_spat_scale = [char(10) 'You have not indicated resolution (cm/pixel). ' ...
                        'Do so by assigning object.spatial_scale a value.'];
                case 0
                    str_spat_scale = [char(10) 'The resolution is set to ' num2str(self.spatial_scale) ...
                        ' cm/pixel.'];
            end
            
            disp(str_spat_scale);            
            disp(str_headdir);
            disp(str_pos);        
            
        end
        
        function self = SetEpoch(self)
        % root = root.SetEpoch;
        %
        % GUI for setting epoch start and stop times for plotting, etc
           
            width = 425;
            height = 600;

            %  Create and then hide the GUI as it is being constructed.
            f = figure('Visible','off','Position', [50, -150, width, height],'Color', 'w');

            sizes = [15, 180, 260];
            widths = [150, 75, 150];
            lineheight = 30;
            
            check_ind = 1;

            %  Construct the components.
            
            %labels
            hlabel1 = uicontrol('Style','text','String','Anchoring Event',...
                  'Position',[sizes(1),height-50,widths(1),lineheight], 'BackgroundColor', 'w',...
                  'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
            hlabel2 = uicontrol('Style','text','String','Lag (+/- seconds)',...
              'Position',[sizes(2),height-50,widths(2),lineheight], 'BackgroundColor', 'w',...
              'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
            hlabel3 = uicontrol('Style','text','String','Secondary Event',...
              'Position',[sizes(3),height-50,widths(3),lineheight], 'BackgroundColor', 'w',...
              'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
              
            str_help = ['If you specify start and stop above, we will ignore the selected labels.' ...
                ' Start and stop must lie between 0 and 1, and indicate percent session time.'];
          
            helplabel = uicontrol('Style','text','String',str_help,...
              'Position',[5,65,width-10,60], 'BackgroundColor', 'w',...
              'FontSize', 10, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
          
            h_objarr(check_ind, 1) = uicontrol('Style','popupmenu',... % one popup
                  'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                  'Position',[sizes(1),height-50-check_ind*lineheight,widths(1),lineheight],...
                  'HandleVisibility', 'off');
              
            h_objarr(check_ind, 2) = uicontrol('Style','edit',... % one textbox
                  'String','',...
                  'Position',[sizes(2),height-50-check_ind*lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');

            h_objarr(check_ind, 3) = uicontrol('Style','popupmenu',... % one popup
                  'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                  'Position',[sizes(3)+15,height-50-check_ind*lineheight,widths(3), lineheight],...
                  'HandleVisibility', 'off');
         
            h_objarr(check_ind, 4) = uicontrol('Style','edit',... % one textbox
                  'String','Start',...
                  'Position',[sizes(1),height-50-check_ind*lineheight*1.5,widths(2),22],...
                  'HandleVisibility', 'off');  
              
            h_objarr(check_ind, 5) = uicontrol('Style','edit',... % one textbox
                  'String','Stop',...
                  'Position',[sizes(2)-50,height-50-check_ind*lineheight*1.5,widths(2),22],...
                  'HandleVisibility', 'off');

              
            % buttons for add epoch or finish
            
            h_add = uicontrol('Style','pushbutton','String','Add',...
                  'Position',[width-100, height-50-(check_ind+1)*lineheight*1.5, 50, 25],...
                  'Callback',{@add_Callback}, 'HandleVisibility', 'off'); 
            h_finish = uicontrol('Style','pushbutton','String','Done',...
                  'Position',[width-50, height-50-(check_ind+1)*lineheight*1.5, 50, 25],...
                  'Callback',{@finish_Callback}, 'HandleVisibility', 'off');

            % Initialize the GUI.
            % Change units to normalized so components resize 
            % automatically.

            % set([f; h_objarr; h_add; h_finish; hlabel1; hlabel2; hlabel3],'Units','Pixels');

            %ha = axes('Position', [.1 .1 .6 .6]);

            % Assign the GUI a name to appear in the window title.
            set(f,'Name','SetEpoch')
            % Move the GUI to the center of the screen.
            movegui(f,'center')
            % Make the GUI visible.
            set(f,'Visible','on');
            
            uiwait(f); % wait to return until figure is closed
            
            function popup_lag_Callback(source, eventdata) 
          
                % Determine the selected data set.
                % clear old plots, and update header

                str = get(source, 'String');
                val = get(source,'Value');
                
            end
            
            function add_Callback(source, eventdata) % add another line for epoch
                
                check_ind = check_ind+1;
                
                h_objarr(check_ind, 1) = uicontrol('Style','popupmenu',... % one popup
                      'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)), char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                      'Position',[sizes(1),height-15-check_ind*lineheight*2,widths(1), lineheight],...
                      'HandleVisibility', 'off');
              
                h_objarr(check_ind, 2) = uicontrol('Style','edit',... % one textbox
                      'String','',...
                      'Position',[sizes(2),height-15-check_ind*lineheight*2+8,widths(2),22],...
                      'HandleVisibility', 'off');

                h_objarr(check_ind, 3) = uicontrol('Style','popupmenu',... % one popup
                      'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                      'Position',[sizes(3),height-15-check_ind*lineheight*2,widths(3), lineheight],...
                      'HandleVisibility', 'off');
                  
                h_objarr(check_ind, 4) = uicontrol('Style','edit',... % one textbox
                  'String','Start',...
                  'Position',[sizes(1),height-15-check_ind*lineheight*2-lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');  
              
                h_objarr(check_ind, 5) = uicontrol('Style','edit',... % one textbox
                  'String','Stop',...
                  'Position',[sizes(2)-50,height-15-check_ind*lineheight*2-lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');
                  
                set(h_add, 'Position', [width-100, height-15-(check_ind)*lineheight*2 - lineheight+8, 50, 25]);
                set(h_finish, 'Position', [width-50, height-15-(check_ind)*lineheight*2 - lineheight+8, 50, 25]);   
                
            end
            
            function finish_Callback(source, eventdata) % set self.epoch and close window
                
                epoch_ind = 1;
                
                dur = self.b_ts(end)-self.b_ts(1); % duration of recording

                for i = 1:size(h_objarr,1)
                    
                    user_start = get(h_objarr(i,4), 'String');
                    user_stop = get(h_objarr(i,5), 'String'); 
                    
                    if ~strcmp(user_start, 'Start') && ~strcmp(user_stop, 'Stop') %start and stop are declared and valid, 
                        
                        user_start = str2num(user_start);
                        user_stop = str2num(user_stop);
                        
                        if user_start<user_stop && user_start<=1 && user_start>=0 && user_stop<=1 && user_stop>=0
                            epoch(epoch_ind,:) = [dur*user_start+self.b_ts(1), dur*user_stop+self.b_ts(1)];
    
                            epoch_ind = epoch_ind+1;
                        else
                            disp('Start and Stop must lie between 0 and 1. Defaulting to 0 and 1 (start and end of session)');
                        end
                        
                    else
                    
                           
                        t1ind = get(h_objarr(i,1), 'Value');
                        t1 = self.event{t1ind,2};

                        if ~isempty(get(h_objarr(i,2), 'String'))
                            t2 = t1+str2num(get(h_objarr(i,2), 'String'));
                            if t2<0
                                disp('One of the timestamps was less than zero. Setting to zero');
                                t2=0;
                            end
                        end

                        if ~exist('t2', 'var')                        
                            t2ind = get(h_objarr(i,3), 'Value');
                            t2 = self.event{t2ind,2};
                            clear t2ind
                        end

                        if t1~=t2
                            epoch(epoch_ind,1) = min([t1, t2]);
                            epoch(epoch_ind,2) = max([t1, t2]);
                            epoch_ind = epoch_ind+1;
                        end
                        clear t1 t2 ind 
                        
                    end
                end
                
                self.epoch = epoch;
                
                close(f)
                
            end
            
        end
        
        
        function Save(root, params) 
        % (1) root.Save
        % (2) root.Save(params)
        %
        % Saves Session class object to root.name_formatted
        %
        % (1) Saves root to root.name, asks if overwrite is necessary. If
        % LFP exists, asks if it should be cleared.
        %
        % (2) Saves root to root.name, but will not prompt user as per
        % params:
        %
        % params = [ <clear_lfp> , <noprompt> ] for batch processing ex.
        % root.Save([1 0]) will clear the LFP, but not overwrite root.name,
        % if it exists
        %
        % If noprompt == 1, then root.Save will attempt to save to
        % root.name. If Directory does not exist, it will create it. If
        % file exists, it will overwrite.
        %
        % If you want to save to a new filename, change root.name before
        % calling root.save
        
            if ~exist('params', 'var')
                params = [NaN NaN];
            elseif numel(params)==1
                params(2) = NaN;
            end
        
        
            if ~isempty(root.b_lfp)
                
                if params(1)==1
                    root = root.ClearLFP;
                    warning('CMBH:notify', '%s', 'Cleared LFP from object. If you want to save it in the future, use root.Save(0)');
                elseif params(2)~=1
                    user_entry = input('Would you like to save the LFP data with the object? This could make for a large filesize. [y/n]: ', 's');
                
                    while( ~strcmpi(user_entry, 'y') && ~strcmpi(user_entry, 'n') )
                        user_entry = input('Would you like to save the LFP data with the object? This could make for a large filesize. [y/n]: ', 's');
                    end
                    
                    if strcmpi(user_entry, 'n')
                        root = root.ClearLFP;
                        warning('CMBH:notify', '%s', 'Cleared LFP from object.');
                    end
                end
                
            end
                        
           
            if exist(root.name_formatted, 'file') && (isnan(params(2)) || params(2)==0) % if exists and user didnt specify what to do
                reply = input(['You are about to overwrite ' root.name_formatted '. Are you sure? [y/n]: '], 's');
                if strcmpi(reply, 'y')            
                    save(root.name_formatted, 'root', '-v7.3');
                    warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
                else
                    warning('CMBH:notify', '%s', 'Save aborted.');
                end               
            elseif exist(root.name_formatted, 'file') && params(2)==1 % if exists and user specified overwrite
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            elseif exist(root.name_formatted(1:find(root.name_formatted==filesep, 1, 'last')-1), 'dir') % if directory exists, but file doesnt, save it
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            elseif params(2)~=1 % if directory doesnt exist and noprompt is false, ask user
                disp('Folder in root.name not found. Please select where you would like to save the file.');
                [file, path] = uiputfile('*.mat', 'Save Object As', pwd);
                root.name = fullfile(path, file);
        
                save(root.name_formatted, 'root','-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            else % directory doesnt exist, but noprompt is true, so mkdir and save file
                
                mkdir(root.name_formatted(1:find(root.name_formatted==filesep, 1, 'last')-1));
        
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
                
            end
            
        end
        
        function Session_array = SplitSession(self, varargin)
        % (1) Session_array = root.SplitSession;
        % (2) root.SplitSession(params);
        %
        % This method will split the present Session instance into more
        % than one Session Objects.
        %
        % (1) Splits the session object as per root.epochs and optionally
        % returns 'Session_array' of session objects of length(size(root.epoch,1))
        %
        % (2) Takes optional parameters
        % 'filenames' -  cellarray of strings, filenames to make each root.name
        % 'epoch' - Nx2 array of times to split. N = length off 'filenames'
        % 'save' - 1 or 0 (default 0). Save each session object
        %
        % andrew 24 may 2010
               
        import CMBHOME.*
        
        success = 0; %#ok<NASGU> % initialize variables
        Session_array = Session;
        
        p = inputParser;

        p.addRequired('self');
        p.addParamValue('filenames', [], @(x) isccellstr(x)); 
        p.addParamValue('epoch',  [], @(x) size(x,2)==2);
        p.addParamValue('save', 0, @ischar);

        p.parse(self, varargin{:});

        filenames = p.Results.filenames;
        epoch = p.Results.epoch;
        save = p.Results.save;
        
        if isempty(epoch)
            epoch = self.epoch;
        end
        
        if size(epoch,1)<2 % if there arent any epochs passed
            help CMBHOME.Session.SplitSession
            return;
        end
        
        if isempty(filenames)
            
            append = 1:size(epoch,1);
            
            filenames = cat(2, repmat(strrep(self.name, '.mat', ''), size(epoch,1), 1), repmat('_Split', size(epoch,1), 1), int2str(append'), repmat('.mat', size(epoch,1), 1));
            
            filenames = mat2cell(filenames, ones(size(filenames,1),1), size(filenames,2));

        end
        
        if length(filenames)~=size(epoch,1)
            help CMBHOME.Session.SplitSession
            error('length of filenames does not match number of epochs');
        end
        
        for i = 1:length(filenames) % parse through all Session Objects
            
            self.epoch = epoch(i,:);
            
            self.cell_thresh = [ 0 0 ];
            
            event = self.event( [self.event{:,2}]>= epoch(i,1) & [self.event{:,2}]<=epoch(i,2), : );
        
            Session_array(i) = Session('name', filenames{i},...
                                    'b_ts', self.ts,'b_x',self.x, 'b_y', self.y, 'b_headdir',...
                                    self.headdir, 'event', event, 'raw_pos', self.raw_pos, ...
                                    'raw_headdir', self.raw_headdir, 'date_created', now, ...
                                    'epoch', epoch(i,:), 'fs_video', self.fs_video,...
                                    'path_lfp', self.path_lfp, 'path_raw_data', self.path_raw_data); % put it all together

            Session_array(i).user_def = self.user_def; % copy over all user_def fields                    
                                
            Session_array(i).user_def.SplitHistory = ['This Session object was split from ' self.name ' using SplitSession.'];                   
                                
            for j = 1:size(self.cells,1) % import cells
                
                spike = self.spike(self.cells(j,1), self.cells(j,2));
                
                spike.ts = spike.ts(spike.ts>=epoch(i,1) & spike.ts<=epoch(i,2));
                
                spike.i = Utils.SpkInds(self.ts, spike.ts);
               
                Session_array(i).spike(self.cells(j,1), self.cells(j,2)) = spike;
                
            end
            
            if save, Session_array(i).Save; end % save, if specified                                
            
        end           
            
    end
    end
    
    properties (Hidden) % old or outdated properties, around for backwards compatibility
    
            b_kalmanvel
            
    end
end

function [i_ts, i_i] = IsolateEpoch(ts, i, epoch)
% returns indeces from vector i which correspond to timestamps in ts within
% epoch start and stop times t_start and t_stop
%
% example: we have indeces i that refer to vector ts of timestamps, give me just the i which
% lie betweeh t_start and t_stop?
%
% also: i_i are indeces in i for which ts is good (used in theta vector

% andrew 18 march 2010

if nargin<3
    error('IsolateEpoch: not enough input arguments');
end

if max(i)>length(ts)
    
    warning('CMBH:error', 'Your indeces do not align with the timestamp vector');
    return
    
end

t = ts(i);

if numel(epoch)>2 && numel(t)>1
    epoch = mat2cell(epoch, ones(1, size(epoch,1)), 2); 
		
		% if a stable sampling frequency can be identified, we can be way fast
		dts = round((unique(diff(t)))*1000); % ms per sample
        
    if ~any(floor(abs(dts-median(dts)))) && length(dts)>3 % all dts are within 1 ms of the median  
			
        fs = round(1000 / median(dts));
        endtime= cellfun(@(a) ceil(fs*a(2))>i(end),epoch);
        sttime = cellfun(@(a) ceil(fs*a(1))<i(1),epoch);
        epoch(endtime)=cellfun(@(a) [a(1) i(end)/fs],epoch(endtime),'uni',0);
        epoch(sttime)=cellfun(@(a) [i(1)/fs a(2)],epoch(sttime),'uni',0);

         i_i=cellfun(@(a) [ceil((a(1)-t(1))*fs+1):ceil((a(2)-t(1))*fs)]',epoch,'uni',0); %so now we can generate our own indices
      
    else
            i_i = cellfun(@(a) find(t>=a(1) & t<=a(2)), epoch, 'Unif', 0); %otherwise do the brute force (100x slower)
    end
		
    tmp_i_i = cell2mat(i_i); tmp_i_i(tmp_i_i==0) = [];
    i_ts=i(tmp_i_i); %these steps avoid doing the really slow thing twice
    i_ts=  mat2cell(i_ts(:),cellfun(@length,i_i),1);

else
    t_start = epoch(1);
    t_stop = epoch(2);

    i_ts = i(t>=t_start & t<=t_stop);
    
    i_i = find(t>=t_start & t<=t_stop);
end

end
function self = SetLFPInds(self)

    if ~iscell(self.p_ind)
        self.p_lfp_ind = cellfun(@(c) nan(size(c)), {self.p_ind}, 'unif', 0); % will eventually be the aligned indices in the LFP
    else
        self.p_lfp_ind = cellfun(@(c) nan(size(c)), self.p_ind, 'unif', 0); % will eventually be the aligned indices in the LFP
    end
                    
end

function self = SetCelInds(self)

    cel = self.cel;
    
    ind = cell(size(self.epoch,1), size(cel,1));
    
    spkind = cell(size(self.epoch,1), size(cel,1));
    
    ind_lfp = cell(size(self.epoch,1), size(cel,1));

    for j = 1:size(cel,1) % now set root.p_cel_ind

        [i, i2] = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);
        %[i, i2] = IsolateEpoch(self.spike(self.cel(1),self.cel(2)).ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);
        
        if iscell(i)

            if ~isempty(self.active_lfp)
                if isempty(self.spike(cel(j,1), cel(j,2)).i_lfp)
                  self = self.AlignSpike2LFP;
                end
                ind_lfp(:, j) = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{self.active_lfp}(c), i2, 'unif', 0); 
            else
                ind_lfp(:, j) = cell(length(i),1);
            end

            ind(:, j) = i; % video indices
            
            spkind(:,j) = i2;

        else

            if ~isempty(self.active_lfp)
                if isempty(self.spike(cel(j,1), cel(j,2)).i_lfp)
                  self = self.AlignSpike2LFP;
                end
                ind_lfp{1, j} = self.spike(cel(j,1), cel(j,2)).i_lfp{self.active_lfp}(i2); 
            else
                ind_lfp{1, j} = [];
            end

            ind{1, j} = i; % video indices
            
            spkind{1,j} = i2;

        end

    end % now we have cell arrays like {epochs, cells} of indices in video tracking and lfp

    self.p_cel_ind = ind; % cell array like {epochs, cells} of tracking alignment

    self.p_cel_lfp_ind = ind_lfp; % cell array like {epochs, cells} of LFP alignment
            
    self.p_cel_spkind = spkind;
    
end

function l = CheckBaseVarLength(self)

    l(1) = length(self.p_b_x);
    l(2) = length(self.p_b_y);
    l(3) = length(self.p_b_ts);
    l(4) = length(self.p_b_vel);
    l(5) = length(self.p_b_headdir);
    l(6) = length(self.p_b_myvar);

    l = max(l);

end
