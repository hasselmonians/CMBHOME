classdef Import
    
    methods (Static)
        
        function varargout = NL(varargin)

            import CMBHOME.NL2Mat.* % imports NL import utilities for windows only
            import CMBHOME.Utils.*
            import CMBHOME.Session
            import CMBHOME.Spike

            global noprompt

            if nargout>=1, varargout{1} = 0; end % returns whether successful save, if requested
            if nargout>1, varargout{2} = cell(0,0); end

            % Imports data from Neuralynx systems to instances of classes
            % Session and Spike,

            % ARGUMENTS

            % varargin uses input parser, and understands the following parameters:
            % 'n_tetrodes' int, default=12
            % 'base_path', str, default = current directory
            % 'save_path', str, default = current directory
            
            % RETURNS (optional)
            % [success, failures] where success is 1 or 0 for each import,
            % and failures are folder names in batch mode which failed
            
            p = inputParser;

            p.addParamValue('base_path', '', @(c) ischar(c));
            p.addParamValue('save_path', pwd, @(c) ischar(c));
            p.addParamValue('fix_pos', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('fix_headdir', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('batch', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('n_tetrodes', 12, @(c) isnumeric(c) & mod(c,1)==0 & numel(c)==1);
            p.addParamValue('add_lfp', 1, @(c) numel(c)==1 && (c==1 || c==0));

            p.parse(varargin{:});

            base_path = p.Results.base_path;
            save_path = p.Results.save_path;
            fix_pos = p.Results.fix_pos;
            fix_headdir = p.Results.fix_headdir;
            batch = p.Results.batch;
            n_tetrodes = p.Results.n_tetrodes;
            add_lfp = p.Results.add_lfp;

            if batch    % if batch is selected, process each folder in base_path
                         % if directory cluster_files exists in each folder,
                         % no user input will be necessary
                if isempty(base_path), error('For batch mode you must pass the base_path argument'); end
                
                noprompt=1;
                folders = dir(base_path);

                for folderi = 1:length(folders)
                    if folders(folderi).isdir && ~strcmp(folders(folderi).name,'.') && ~strcmp(folders(folderi).name,'..') % for every folder, call Import.NL

                        tmp_base_path = fullfile(base_path, folders(folderi).name);
                        
                        import CMBHOME.*

                        success = Import.NL('base_path', tmp_base_path, ...
                            'save_path', tmp_base_path, 'n_tetrodes', n_tetrodes, ...
                            'fix_headdir', fix_headdir, 'fix_pos', fix_pos);

                        if ~success
                            if nargin>1
                                varargout{2}{end+1} = tmp_base_path;
                                disp([tmp_base_path ' skipped.']);
                            end
                        else
                            disp([tmp_base_path ' imported successfully. ' folders(folderi+1).name ' is next.']);
                        end
                        
                    end
                end
                
                if nargout>1
                    if ~isempty(varargout{2})
                        disp('There were folders skipped. See second return variable for details');
                    end
                end

                return;
             end
            
            %Prompt the user to specify the name and location of the .mat
            %file to be created
            [save_file, save_path] = SaveObject(save_path);
            
            if exist(fullfile(save_path, save_file), 'file')
                if nargout==1
                    varargout{1} = 1;
                end
                return;
            end
            
            if isempty(base_path), base_path = save_path; end
            
            %% Import behavioral data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Prompt the user to specify the Neuralynx raw video file if one cannot be found automatically
            [video_file, base_path] = LoadVideo(base_path);
            if isempty(video_file), return, end

            %Prompt the user to specify the Neuralynx events file if one cannot be found automatically
            [events_file, base_path] = LoadEvent(base_path);
            if isempty(events_file), return, end

            %% Load Raw Video File
            % Note that these vectors may need error checking. These can be performed
            % later.
            
            if isempty(strfind(video_file,'mat'))
                FieldSelectionArray = [1 1 1 1 0 0];
                ExtractHeader = 0;
                ExtractionMode = 1;
                [t, x, y, headdir] = Nlx2MatVT(video_file, FieldSelectionArray,ExtractHeader,ExtractionMode);
                t = t/10^6; %convert from microseconds to seconds
            
            else
                t=load(video_file);
                x=t.x; y=t.y;headdir=t.hd;t=t.ts;
            end
                
            %% Load Event Files
            FieldSelectionArray = [1 0 0 0 1];
            ExtractHeader = 0;
            ExtractionMode = 1;
            [event_ts, event_str] = Nlx2MatEV(events_file, FieldSelectionArray,ExtractHeader,ExtractionMode);
            event_ts = event_ts';

            % Convert times from microseconds to seconds
            event_ts = event_ts/10^6;

            event = cell(length(event_ts), 2);
            event(:,1) = event_str;
            event(:,2) = num2cell(event_ts);

            root = Session('name', save_file, ...
                                'b_ts', t,'b_x',x, 'b_y', y, 'b_headdir',...
                                headdir, 'event', event, 'raw_pos', 1, ...
                                'raw_headdir', 1, 'date_created', now, ...
                                'epoch', [t(1) t(end)], 'fs_video', mean(diff(t)).^-1); % put it all together

            if fix_pos
                root = root.FixPos();
            end

            if fix_headdir
                root = root.FixDir();
            end
            
            %% Import spike data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ~exist('root', 'var'), root = Session(); end
                
            root.spike = Spike();

            [clusters, fname] = getClusterFiles(base_path, n_tetrodes);
            if isempty(clusters), return, end            
            
            for tetrode = 1:n_tetrodes;

                spikes = clusters{tetrode}; % [Nx2] array like [cell number, spiketime;...]

                if ~isempty(spikes)

                    n_cells = max(spikes(:,1)); % 0 indicates usorted spikes

                    for j = 1:n_cells;
                        % Add each cell's spikes to the struct
                        status_message = ['tetrode ' num2str(tetrode) ', cell ' num2str(j)];
                        disp(status_message)

                        spk_ts = spikes(spikes(:,1)==j,2);

                        spk_i = SpkInds(t, spk_ts);
                        
                        if ~isempty(spk_ts)

                            root.spike(tetrode, j) = Spike('i', spk_i(:), 'ts', spk_ts, 'tet_name', fname{tetrode});
                        
                        else
                            
                            disp('Skipped because there were no spikes.');
                        
                        end

                     end

                end

            end
            
            %% Offer to add LFP data files to the path_lfp property &&&&&&
            
            if noprompt | add_lfp
                
                lfp_files = dir(fullfile(base_path, '*.ncs'));
                
                lfp_files = {lfp_files(:).name};
                
                pathname = base_path;
                
            else
            
                [lfp_files, pathname] = uigetfile({'*.mat'; '*.ncs'; '*,plx'}, 'Would you like to make this object aware of LFP (EEG) files?', base_path,  'MultiSelect', 'on');

            end
            
            if iscell(lfp_files)
                root.path_lfp = cellfun(@(c) cat(2, pathname, c), lfp_files, 'UniformOutput', false);
            else
                root.path_lfp{1} = cat(2, pathname, lfp_files);
            end
            
            %% Add waveform to user_def %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plx = dir(fullfile(base_path, '*.plx'));
            disp('importing waveforms, please wait ...')
            waveform.mean = []; waveform.std = [];
            waveform = repmat(waveform,size(root.cells,1),4);
            for c = 1:size(root.cells,1)
                ttfile = fullfile(base_path, ['TT' num2str(root.cells(c,1)) '.plx']);
                cutNumber = root.cells(c,2);
                for k = 1:4
                    [~,~,~,wv] = CMBHOME.PLX2Mat.plx_waves_v(ttfile,k,cutNumber);
                    waveform(c,k).mean = mean(wv);
                    waveform(c,k).std = std(wv);
                end
            end
            
            root.user_def.waveform = waveform;

            %% Add files used and Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            datafiles = cat(2, video_file, events_file, shiftdim(fname,1), lfp_files);
            
            root.path_raw_data = datafiles;

            save(fullfile(save_path, save_file), 'root');
            
            if nargout==1
                varargout{1} = 1;
            end
%             
            function [clusters, fname_ordered] = getClusterFiles(base_dir, n_tetrodes)

            % looks through path for Sc* files that are either mat or ntt, and returns
            % a cell array of [Nx2] array like [cell number, spiketime;...]

                clusters = [];
                fname_ordered = [];
            
                import CMBHOME.*
                import CMBHOME.NL2Mat.*

                cluster_file_path = fullfile(base_dir, 'ClusterCut');

                cluster_files = cell(1,1);

                pathname = base_dir;
                
                if exist(cluster_file_path, 'dir')

                    cluster_files = dir(cluster_file_path);

                    cluster_files(strcmp('..', {cluster_files(:).name})) = []; % remove up dir character
                    cluster_files(strcmp('.', {cluster_files(:).name})) = []; % remove top dir character                

                    cluster_files = strcat(['ClusterCut' filesep], {cluster_files(:).name});
                else
                    
                    % assume cut files look like Sc*.mat
                    
                    a = dir(fullfile(base_dir, 'Sc*.mat'));
                    
                    cluster_files = {a(:).name};
                    
                    if isempty(cluster_files)
                        cluster_files{1} = [];
                    end
                    
                end

                if isempty(cluster_files{1}) || length(cluster_files)>n_tetrodes

                    if noprompt, disp('Batch mode must find a directory ''ClusterCut'' in every session folder with no more than n_tetrodes files, or follow naming convention Sc*.mat.'); return; end

                    [cluster_files, pathname] = uigetfile({'*.mat'; '*.ntt'}, 'Cluster File Selector', base_dir,  'MultiSelect', 'on');
                end

                clusters = cell(n_tetrodes,1);

                if length(cluster_files)>n_tetrodes
                    disp('Import.NL found more cluster files than n_tetrodes. Make sure you set n_tetrodes to the proper number to avoid this message. If it persists, something is wrong.');
                    n_tetrodes = length(cluster_files);
                end

                fname_ordered = TetInds(cluster_files, n_tetrodes); % collects and determines which clusters files belong to which tetrode
                ext = cell(size(fname_ordered));
                for i = 1:length(fname_ordered)
                    if length(fname_ordered{i})>3
                        ext{i} = fname_ordered{i}(end-3:end);
                    end
                end

                %if ~isnumeric(indeces), error('Error in getClusterFiles in Import'); end

                for i = 1:length(fname_ordered)
                    if ~isempty(fname_ordered{i})

                    disp([fname_ordered{i} ' is index ' num2str(i)]);

                        switch ext{i}
                            case '.mat'

                                tmpvar = whos('-file', fullfile(base_path, fname_ordered{i}));

                                if length(tmpvar)>1, disp('Badly formatted cluster file');  return; end

                                tmp = load(fullfile(pathname, fname_ordered{i}));

                                clusters{i} = tmp.(tmpvar.name);
                                
                            case '.ntt'

                                FieldSelectionArray = [1 0 1 0 0];
                                ExtractHeader = 0;
                                ExtractionMode = 1;

                                [timestamps, cell_ind] = Nlx2MatSpike(fullfile(pathname, fname_ordered{i}), FieldSelectionArray,ExtractHeader,ExtractionMode);

                                clusters{i} = [shiftdim(cell_ind), shiftdim(timestamps)/10^6;]; % convert time from microseconds to seconds
                        end

                    end

                end

                if all(cellfun(@(c) isempty(c), clusters)), disp('No cut cluster files with cells were found, but there were cells expected'); end

            end

            function fnames_ordered = TetInds(fnames, n_tetrodes)

            % returns n_tetrodes by 1 cell array of filenames corresponding to tetrode
            % of index they reside
                fnames_ordered = cell(n_tetrodes,1);
                matches = zeros(length(fnames), n_tetrodes);
                
                if ischar(fnames), fnames = {fnames}; end % if only one tetrode file
                
                for i=1:n_tetrodes
                    matches(:,i) = cellfun(@(c) sum(ismember(num2str(i), c))/length(num2str(i))==1 && ...
                        sum(ismember(c, num2str(i)))/length(num2str(i))==1, fnames);
                end
                matches = matches==1;

                for  i=1:length(fnames)

                    ind = find(matches(i,:)==1, 1, 'last');

                    if ~isempty(ind)
                        fnames_ordered{ind} = fnames{i};
                    end
                end
            end    

            function [save_file, save_path] = SaveObject(save_path)

                [pathstr, filename, ext] = fileparts(save_path);
                
                if noprompt & isempty(ext)
                    ind = strfind(save_path, filesep);

                    session_name = save_path(ind(end)+1:end);

                    save_file = ['CMBH_' session_name '.mat'];
                    
                elseif noprompt
                    
                    save_path = pathstr;
                    
                    save_file = [filename ext];
                    
                else
                    [save_file,save_path] = uiputfile('*.mat','Where would you like to save the session data?',save_path);
                end
            end

            function [video_file, base_path] = LoadVideo(base_path)
                
                video_file = [];
                
                video_files = [dir(fullfile(base_path, '*.nvt'));dir(fullfile(base_path, '*.mat'))];
                
                if noprompt & length(video_files)>1
                    
                    video_files = dir(fullfile(base_path, '*OrdCor.nvt'));
                    
                    if isempty(video_files)
                        video_files = dir(fullfile(base_path, '*.nvt'));
                        for i =1:length(video_files)
                            name_size(i) = length(video_files(i).name);
                        end
                        [trash, tmp_ind] = min(name_size);
                        video_files = video_files(tmp_ind).name;
                    end
                end        

                if length(video_files)~=1 

                    if noprompt, disp('Batch process could not find proper video file'); return; end

                    [load_file, base_path] = uigetfile('*.Nvt;*.nvt;*.mat','Please select Neuralynx video (.Nvt) file to load', base_path);

                    if load_file
                        video_file = fullfile(base_path, load_file);
                    else
                        error('No video file selected.');
                    end
                else
                    video_file = fullfile(base_path, video_files.name);
                end

            end

            function [events_file, base_path] = LoadEvent(base_path)

                events_files = dir(fullfile(base_path, 'Events.nev'));

                if length(events_files)~=1 

                    if noprompt, disp('Batch process could not find proper event file');  events_file = []; return; end

                    [load_file,base_path] = uigetfile('*.Nev;*.nev','Please select Neuralynx events (.Nev) file to load',base_path);

                    if load_file
                        events_file = fullfile(base_path, load_file);
                    else
                        error('No event file selected.');
                    end
                else
                    events_file = fullfile(base_path, events_files.name);
                end

            end
        
        end
        
        function varargout = PLX(varargin)
            
            % to do: head direction vector?
            
            import CMBHOME.PLX2Mat.* % imports NL import utilities for windows only
            import CMBHOME.Utils.*
            import CMBHOME.Session
            import CMBHOME.Spike

            p = inputParser;

            p.addParamValue('base_path', '', @(c) ischar(c));
            p.addParamValue('fix_pos', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('fix_headdir', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('batch', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('save_file', 1, @(c) numel(c)==1 && (c==1 || c==0));
           
            p.parse(varargin{:});

            base_path = p.Results.base_path;
            fix_pos = p.Results.fix_pos;
            fix_headdir = p.Results.fix_headdir;
            batch = p.Results.batch;
            save_file = p.Results.save_file;          
            
            [fopen, fsave] = CollectPLXFiles(batch,base_path); % cell array of directories/filenames to be imported
            
            if save_file, ca_roots = cell(size(fopen,1), 1); end
            
            for ii = 1 : size(fopen,1) % make objects, one by one
                                
                % this script reads all the spike timestamp and a/d info from a plx file into matlab
                % variables.

                [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information([fopen{ii,:}]);

                disp(['Opened File Name: ' OpenedFileName]);
                disp(['Version: ' num2str(Version)]);
                disp(['Frequency : ' num2str(Freq)]);
                disp(['Comment : ' Comment]);
                disp(['Date/Time : ' DateTime]);
                disp(['Duration : ' num2str(Duration)]);
                disp(['Num Pts Per Wave : ' num2str(NPW)]);
                disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
                % some of the information is only filled if the plx file version is >102
                if ( Version > 102 )
                    if ( Trodalness < 2 )
                        disp('Data type : Single Electrode');
                    elseif ( Trodalness == 2 )
                        disp('Data type : Stereotrode');
                    elseif ( Trodalness == 4 )
                        disp('Data type : Tetrode');
                    else
                        disp('Data type : Unknown');
                    end

                    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
                    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
                    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
                    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
                end   
                
                % get some counts
                [tscounts, wfcounts, evcounts] = plx_info(OpenedFileName,1);

                % tscounts, wfcounts are indexed by (unit+1,channel+1)
                % tscounts(:,ch+1) is the per-unit counts for channel ch
                % sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
                % [nunits, nchannels] = size( tscounts )
                % To get number of nonzero units/channels, use nnz() function

                % get events
                
                % header information, including event strings and
                % corresponding channel numbers
                headerinf = readPLXHeaders(OpenedFileName);                
                
                event_channels = [headerinf.evchans(:).channel]; % channel numbers
                evnames = {headerinf.evchans(:).name}; % channel names
                
                evnames(event_channels>256) = []; % only want those less than 257
                event_channels(event_channels>256) = [];
                
                event_channels(strcmp(evnames, 'Frame Marker')) = []; % get rid of frame markers as well
                evnames (strcmp(evnames, 'Frame Marker')) = []; 
                
                for iev = 1:299
                    if ( evcounts(iev) > 0 )
                        if ( iev == 257 )
                            % treat strobed channel seperately, just to avoid setting up a
                            % cell array to hold strobed values when only one channel will
                            % have them.
                            [nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(OpenedFileName, iev); 
                        else
                            [nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(OpenedFileName, iev);
                        end
                    end                    
                end
                                            
                x = [];
                y = [];
                t = [];
                headdir = [];
                
                % build event cell array
                event = cell(sum([nevs{event_channels}]), 2);

                if ~isempty(event_channels)
                    tmp_ind = 1;
                    for i_ch = 1:length(event_channels)
                        if nevs{event_channels(i_ch)}>0
                            event(tmp_ind:tmp_ind+nevs{event_channels(i_ch)}-1,1) = evnames(i_ch); % set labels
                            
                            event(tmp_ind:tmp_ind+nevs{event_channels(i_ch)}-1,2) = num2cell(tsevs{event_channels(i_ch)}); % set times
                                                        
                            tmp_ind = tmp_ind+nevs{event_channels(i_ch)}; % update index
                        end
                    end
                    
                    [~, tmpind] = sort([event{:,2}]); % sort by time

                    event = event(tmpind, :);
                    
                else
                    event = cell(1,2);
                end
                
                % get some other info about the spike channels
                [nspk,spk_filters] = plx_chan_filters(OpenedFileName);
                  [nspk,spk_gains] = plx_chan_gains(OpenedFileName);
                [nspk,spk_threshs] = plx_chan_thresholds(OpenedFileName);
                [nspk,spk_names] = plx_chan_names(OpenedFileName);

                spk_names = cellstr(spk_names);
                
                % get the a/d data into a cell array also.
                % This is complicated by channel numbering.
                % The presence/absence of slow analog data can be seen by looking at the
                % evcounts array at indexes 300-363. E.g. the number of samples for
                % analog channel 0 is stored at evcounts(300).
                % Note that analog ch numbering starts at 0, not 1 in the data, but the
                % 'allad' cell array is indexed by ich+1
                numads = 0;
                for ich = 0:63
                    if ( evcounts(300+ich) > 0 )
                        [adfreq, nad, tsad, fnad, allad{ich+1}] = plx_ad(OpenedFileName, ich);
                        numads = numads + 1;
                    end
                end

                if ( numads > 0 )
                    [nad,adfreqs] = plx_adchan_freqs(OpenedFileName);
                    [nad,adgains] = plx_adchan_gains(OpenedFileName);
                    [nad,adnames] = plx_adchan_names(OpenedFileName);

%                     % just for fun, plot the channels with a/d data
% 
%                     [adrows,nActiveADs] = size(allad);
%                     for ich = 1:nActiveADs
%                         if ( size(allad{ich}) > 0 )
%                             subplot(nActiveADs,1,ich); plot(allad{ich});
%                         end
%                     end
                end
                
                if nevs{257}>0 % check that there is position data
                    
                    [n_coords, trash, nVTMode, c] = plx_vt_interpret(tsevs{257}, svStrobed);    % added by andrew to spit out coordinates, etc
                    
                    % check plx_vt_interpret to learn about nVTMode
                    
                    switch nVTMode
                        case {1,2,3,4,5}
                    
                            y = c(:, 3);
                            x = c(:, 2);
                            t = c(:, 1);
                            
                        case {6,7,8} % there are two sets of coordinates, presumably one for each LED
                            % I take the mean of both sets of coordinates,
                            % so long as one of them isnt 0. if it is, then
                            % i just assume the other is right
                            
                            % if they are both zero, it stays zero and then
                            % hopefully the FixPos can fix it        
                                                        
                            c((c(:,2)==0 & c(:,3)==0),2:3) = NaN; % set zeros to NaN so we ignore them in mean
                            c((c(:,4)==0 & c(:,5)==0),4:5) = NaN;
                            
                            x = nanmean([c(:,2), c(:,4)],2);
                            y = nanmean([c(:,3), c(:,5)],2);
                            
                            x(isnan(x)) = 0;
                            y(isnan(y)) = 0;
                            
                            t = c(:,1);
                                                    
                        otherwise % nothing programmed for 3 leds yet
                            disp('Three LED signals detected, only using the first');
                            y = c(:, 3);
                            x = c(:, 2);
                            t = c(:, 1);

                    end
                    
                else
                    disp('Import.PLX: no tracking data found. This import script assumes you ran CinePlex before import.');
                end

                root = Session('name', [fsave{ii,:}], ...
                                'b_ts', t,'b_x',x, 'b_y', y, 'b_headdir',...
                                nan(length(y),1), 'event', event, 'raw_pos', 1, ...
                                'raw_headdir', 1, 'date_created', now, ...
                                'epoch', [t(1) t(end)], 'fs_video', mean(diff(t)).^-1, ...
                                'path_raw_data', fopen{ii}, 'path_lfp', {[fopen{ii,:}]}); % put all behavioral stuff together

                if fix_pos
                    root = root.FixPos();
                end
                
                if fix_headdir
                    root = root.FixDir();
                end                
                                                 
                clear x y t headdir                
                                
                root.spike = Spike(); % initialize spiking data
                
                % gives actual number of units (including unsorted) and actual number of
                % channels plus 1
                [nunits1, nchannels1] = size( tscounts );

                root.spike(1:(size(tscounts,2)-1)/Trodalness, 1:nunits1-1) = Spike(); % initialize the spike array

                % we will read in the timestamps of all units,channels into a two-dim cell
                % array named allts, with each cell containing the timestamps for a unit,channel.
                % Note that allts second dim is indexed by the 1-based channel number.
                for iunit = 1:nunits1-1   % starting with unit 1 (sorted). 0 is unsorted
                    ich = 1;
                    while ich<=nchannels1-1
                        if ( tscounts( iunit+1 , ich+1 ) > 0 )
                            % get the timestamps for this channel and unit
                            %[nts, allts{iunit+1,ich}] = plx_ts(OpenedFileName, ich , iunit );
                            
                            if lower(spk_names{ich}(1))=='t' % this is a tetrode channel
                                tet_ind = sscanf(spk_names{ich}, '%1s%f%1s%f');
                                tet_ind = tet_ind(2);
                                
                                [trash, spk_ts] = plx_ts(OpenedFileName, ich , iunit );
                                
                                spk_i = shiftdim(SpkInds(root.b_ts, spk_ts));
                                                               
                                if isempty(root.spike(tet_ind, iunit).ts) && ~isempty(spk_i) % if this cell hasnt been loaded (because this cells spikes are repeated four times over in tetrode event land)

                                    disp([spk_names{ich} ' is tetrode ' int2str(tet_ind) ' and cell ' int2str(iunit)]);
                                    root.spike(tet_ind, iunit) = Spike('i', spk_i, 'ts', spk_ts, 'tet_name', spk_names{ich});
                                    
                                end
                                
                                if isempty(spk_i), disp(['Skipping tetrode ' int2str(tet_ind) ' and cell ' int2str(iunit) ' because there were no spikes.']); end
                                
                                ich = ich+1; % skip next three channels, because they are just the identical spikes from other electrodes in the tetrode
                            else
                                ich = ich+1;
                            end
                        else
                            ich = ich+1;
                        end
                    end
                end
                
                if save_file
                    save([fsave{ii,:}], 'root'); % save
                end
                
                if nargout==1
                    ca_roots{ii} = root;
                end
                
            end           
        
            if nargout==1
                varargout{1} = ca_roots;
            end
            
            function [fopen, fsave] = CollectPLXFiles(batch, base_path)

                fopen = cell(1,2); % reasonably large cell array of str fnames
                fsave = cell(1,2);

                if batch        % if batch is selected, process each folder in base_path
                                % if directory cluster_files exists in each folder,
                                % no user input will be necessary
                    if isempty(base_path), error('For batch mode you must pass the base_path argument'); end

                    folders = dir(base_path);
                    
                    folders(end+1).name = base_path; % check base_path, too

                    f_ind = 1;

                    for folderi = 1:length(folders)
                        if folders(folderi).isdir & ~strcmp(folders(folderi).name,'.') & ~strcmp(folders(folderi).name,'..') % for every folder, call Import.NL

                            tmp_base_path = fullfile(base_path, folders(folderi).name);

                            tmp_file = dir([tmp_base_path '*.plx']);
                            
                            tmp_file = {tmp_file(:).name};
                            
                            for i=1:length(tmp_file)

                                fopen(i, :) = {tmp_base_path, tmp_file{i}};
                                fsave(i,:) = {tmp_base_path, ['CMBH_', strrep(tmp_file{i}, '.plx', '.mat')]};
                                
                            end
                       
                        end
                    end
                                        
                    return;
                end % end batch 

                % otherwise, open base_path (pwd if empty) and keep
                % returning files until user says done

                if isempty(base_path), base_path = pwd; end

                f_ind = 1;

                while 1 % loop until user cancels, and we get no files

                    [load_files, base_path] = uigetfile('*.plx','Please select PLX files to load. Exit to finish.',base_path, 'MultiSelect', 'on');

                    if iscell(load_files)
                            for i = 1:length(load_files)
                                fopen(f_ind,:) = {base_path, load_files{i}};
                                fsave(f_ind,:) = {base_path, ['CMBObj_' strrep(load_files{i}, '.plx', '.mat')]};
                                f_ind = f_ind+1;
                            end
                    elseif ischar(load_files)
                            fopen(f_ind,:) = {base_path, load_files};
                            fsave(f_ind,:) = {base_path, ['CMBObj_' strrep(load_files, '.plx', '.mat')]};
                            f_ind = f_ind+1;
                    else
                        break; % exit loop collecting files
                    end

                end

                if isempty(fopen)
                    disp('No .plx files found/selected. No CMBObjects will be saved.');
                end

            end   
        end
    end    
end

