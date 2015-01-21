function [control_f, f] = VisualizeLFP(self)
    
    % load all LFP files into object with waitbar if not already loaded
    
%     if load_lfp && length(self.b_lfp)~=length(self.path_lfp)
%         hwaitbar = waitbar(0,'Loading LFP data...');
%         for i = 1:length(self.path_lfp)
%             self.active_lfp = i;
%             self = self.LoadLFP; 
%             waitbar(i/length(self.path_lfp));
%         end
%         
%         close(hwaitbar);
%         
%     end
  
    width = 500;
    height = 650;
    
    f_low = 0; % Hz
    f_high = 20; % Hz
    bandwidth = .25; %Hz
    windowsize = 10; % seconds
    windowinc = 2; % seconds
    TBP = 10; 
    tapers = 9;
    
    lfp_file = self.path_lfp{1};
    lfp_ind = 1;
        
    loaded_lfp = CheckLoadedLFP(self);
    
    if isempty(self.active_lfp) && any(loaded_lfp), self.active_lfp = find(loaded_lfp==1, 1, 'first'); end
    
    % Create and then hide the Figure Screen as it is being constructed
    f = figure('Visible','off','Position', [50, -150, width, height],'Color', 'w');
    
    control_height = 600;
    control_width = 250;
    header_height = 215;
    
    %  Create and then hide the GUI as it is being constructed.
    control_f = figure('Visible','off','Position', [50, -150, control_width, control_height],'Color', 'w', 'Menubar', 'none');
    
    check_size = [225, 25];
    check_posx = repmat(15, 18, 1);
    
    check_ind = 1;
    
    %  Construct the components.
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Chronux Spectrogram',...
          'Position',[check_posx(check_ind), control_height - header_height - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@cspectrogram_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;  

%     check_objarr(check_ind) = uicontrol('Style','checkbox','String','MATLAB Spectrogram',...
%           'Position',[check_posx(check_ind), control_height - header_height - 25*check_ind, check_size],...
%           'BackgroundColor', 'w','Callback',{@mspectrogram_Callback}, 'HandleVisibility', 'off');
     check_ind = check_ind+1;   
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Chronux Spectrum',...
          'Position',[check_posx(check_ind), control_height - header_height-100 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@cspectrum_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;  
    
%     check_objarr(check_ind) = uicontrol('Style','checkbox','String','MATLAB Spectrum',...
%           'Position',[check_posx(check_ind), control_height - header_height-100 - 25*check_ind, check_size],...
%           'BackgroundColor', 'w','Callback',{@mspectrum_Callback}, 'HandleVisibility', 'off');
%     check_ind = check_ind+1;   
 
    plot_vec = zeros(check_ind-1, 1);
    
%% END OF CHECKBOXES

    hf_low = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-135, 50, 25],...
          'HandleVisibility', 'off', 'String', num2str(f_low));  
    
    uicontrol('Style', 'text', 'String', 'Frequency Range', 'Position',...
        [30, control_height-135, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);

    hf_high = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[180,control_height-135, 50, 25],...
          'HandleVisibility', 'off', 'String', num2str(f_high));

    hwindowsize = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-160, 100, 25],...
          'HandleVisibility', 'off', 'String', num2str(windowsize));
      
      uicontrol('Style', 'text', 'String', 'Time Window Width', 'Position',...
        [30, control_height-160, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);

    hwindowinc = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-185, 100, 25],...
          'HandleVisibility', 'off', 'String', num2str(windowinc));
      
      uicontrol('Style', 'text', 'String', 'Time Window Increment', 'Position',...
        [30, control_height-187, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);

    hbandwidth = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-210, 100, 25],...
          'HandleVisibility', 'off', 'String', num2str(bandwidth));
      
      uicontrol('Style', 'text', 'String', 'Bandwidth', 'Position',...
        [30, control_height-220, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);
    
    % Spectrum controlbox
    
    hTBP = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-320, 100, 25],...
          'HandleVisibility', 'off', 'String', num2str(TBP));
      
      uicontrol('Style', 'text', 'String', 'Time-Bandwidth Product', 'Position',...
        [30, control_height-320, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);
    
    htapers = uicontrol('Style','edit',... % one textbox
          'String','','BackgroundColor', 'w', ...
          'Position',[130,control_height-355, 100, 25],...
          'HandleVisibility', 'off', 'String', num2str(tapers));
      
      uicontrol('Style', 'text', 'String', 'Tapers', 'Position',...
        [30, control_height-355, 75, 25], 'BackgroundColor', 'w', 'FontSize', 8);

    htitle = uicontrol('Style','text','String','',...
          'Position',[15, control_height-25,300,23], 'BackgroundColor', 'w',...
          'FontSize', 16, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off'); 

    hspecs = uicontrol('Style','text','String','',...
          'Position',[15,control_height-85,300,30], 'BackgroundColor', 'w',...
          'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');       

    hpopup = uicontrol('Style','popupmenu',...
          'String',self.path_lfp(logical(loaded_lfp)),...
          'Position',[15, control_height-55, 225, 25],...
          'Callback',{@popup_menu_Callback}, 'HandleVisibility', 'off'); 

    hupdate = uicontrol('Style','pushbutton','String','Update',...
          'Position',[150, control_height-105, 70, 25],...
          'Callback',{@update_Callback}, 'HandleVisibility', 'off');
    
    hepoch = uicontrol('Style','pushbutton','String','Change Epochs',...
          'Position', [15, control_height-105, 125, 25],...
          'Callback', {@epoch_Callback}, 'HandleVisibility', 'off');
    
      

    if ~all(logical(loaded_lfp))
        
        hsave = uicontrol('Style','pushbutton','String','Load LFP',...
            'Position',[100, 5, 150, 25],...
            'Callback',{@load_Callback}, 'HandleVisibility', 'off');
    else
        
        hsave = [];
        
    end
        
   % Initialize the GUI.
   % Change units to normalized so components resize 
   % automatically.
   
    %set([control_f; hpopup; hupdate; check_objarr(:); htitle; hspecs; hepoch]);

    set(f, 'Units','normalized');

    set([control_f; hpopup; hupdate; check_objarr(:); htitle; hspecs; hepoch; hsave],...
    'Units','normalized');

    plot_vec = zeros(25, 1);
    current_cell = self.cells(1,:);

    str_title = self.path_lfp{1};
    set(htitle, 'String', str_title);

    % Update header specs
    
    figure(control_f);
    
    if size(self.epoch,1)>1
        str1 = [int2str(size(self.epoch,1)) ' Epochs Selected'];
        str2 = '';
    else
        str1 = ['Event Start: ' self.active_event{1}];
        str2 = ['Event Stop: ' self.active_event{2}];
    end

    str_specs = str2mat(str1, str2);
    set(hspecs, 'String', str_specs);
    
    %ha = axes('Position', [.1 .1 .6 .6]);
   
    movegui(control_f,'center')

    set(control_f,'Name','VisualizeLFP')

    set(control_f,'Visible','on'); 
    
    function load_Callback(source, eventdata)
       
        hwaitbar = waitbar(0,'Loading LFP data...');
        
        for i = 1:length(self.path_lfp)
            
            if loaded_lfp(i)~=1
            
                self.active_lfp = i;
                self = self.LoadLFP; 
                waitbar(i/length(self.path_lfp));
                
            end
            
        end
        
        if loaded_lfp ~= CheckLoadedLFP(self) % if we were successful loading new LFP
            set(hsave, 'String','Save Session with LFP',...
            'Callback',{@save_Callback}, 'HandleVisibility', 'on');
        else  % if we werent
            delete(source);
        end

        close(hwaitbar);

    end
    
    function save_Callback(source, eventdata)
       
        self.Save([0, 1]); % save by overwriting old Session and do not clear lfp
        
        delete(source); % clear save button
        
    end
    
    function popup_menu_Callback(source, eventdata) 

        % Determine the selected data set.
        % clear old plots, and update header

        str = get(source, 'String');
        val = get(source,'Value');

        % Set current data to the selected data set.

        lfp_file = str;
        lfpinds = find(loaded_lfp);
        lfp_ind = lfpinds(val);

        % Update header title

        h = get(f,'Children');

        delete(h)

        updateScreen;

    end

    function epoch_Callback(source, eventdata)
        self = self.SetEpoch;
        
        h = get(f,'Children');
        
        delete(h)
        
        updateScreen;
    end 

    function update_Callback(source, eventdata)  
        
        hand = get(f,'Children');
        
        delete(hand)
        
        f_low = str2double(get(hf_low, 'String'));
        f_high = str2double(get(hf_high, 'String'));
        bandwidth = str2double(get(hbandwidth, 'String'));
        windowsize = str2double(get(hwindowsize, 'String'));
        windowinc = str2double(get(hwindowinc, 'String'));
        TBP = str2double(get(hTBP, 'String'));
        tapers = str2double(get(htapers, 'String'));

        updateScreen;
               
    end

    function cspectrogram_Callback(source,eventdata)

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(1) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(1) = 0;
        end

    end

    function mspectrogram_Callback(source,eventdata)

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(2) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(2) = 0;
        end

    end

    function cspectrum_Callback(source,eventdata)

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(3) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(3) = 0;
        end

    end

    function mspectrum_Callback(source,eventdata)

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(4) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(4) = 0;
        end

    end

    function updateScreen 

        import CMBHOME.*

        figure(control_f);

        self.active_lfp = lfp_ind;

        str_title = self.path_lfp{lfp_ind};
        set(htitle, 'String', str_title);

        % Update header specs

        if size(self.epoch,1)>1
            str1 = [int2str(size(self.epoch,1)) ' Epochs Selected'];
            str2 = '';
        else
            str1 = ['Event Start: ' self.active_event{1}];
            str2 = ['Event Stop: ' self.active_event{2}];
        end

        str_specs = str2mat(str1, str2);
        set(hspecs, 'String', str_specs);

        % set the 'Position' property data for each axis

        % we need to orient all of our subplot axes 250 pixels down from
        % the top, and we need to normalize them since thats all subplot
        % takes

        figure(f);

        h_max = .95;
        h_min = .05; 
        buffer = .05;

        % make as many subplots as necessary, display them

        if sum(plot_vec)<3

            a_h = (h_max-h_min-(sum(plot_vec)-1)*buffer)/sum(plot_vec);
            a_hs = linspace(h_max-a_h,h_min, sum(plot_vec));

            a_w = .8;
            a_ws = [.1 .1];

        else

            a_h = (h_max-h_min-(ceil(sum(plot_vec)/2-1)*buffer))/ceil(sum(plot_vec)/2);
            a_hs = linspace(h_max-a_h,h_min, ceil(sum(plot_vec)/2));
            a_hs = reshape([a_hs;a_hs], 2*numel(a_hs), 1);

            a_w = .35;
            a_ws = repmat([.1, .5450], 1, ceil(sum(plot_vec)/2));

        end

        fig_ind = 1;

        if plot_vec(14)
            figure(f);
            subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

            axis off   

            str_title = self.path_lfp{lfp_ind};


            % Update header specs

            if size(self.epoch,1)>1
                str1 = [int2str(size(self.epoch,1)) ' Epochs Selected'];
                str2 = '';
            else
                str1 = ['Event Start: ' self.active_event{1}];
                str2 = ['Event Stop: ' self.active_event{2}];
            end

            str3 = self.path_lfp{1};

            str_specs = str2mat(str1, str2, str3);

            text(.1, .8, [str_title ' '], 'FontSize', 14);

            text(.1, .6, str_specs, 'FontSize', 8);

            text(.1, .4, strrep(self.name, '_', '\_'), 'FontSize', 8);

            fig_ind = fig_ind+1;
        end

        if plot_vec(1)
            figure(f);
            subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
            self.plot_lfp_spectrogram(windowsize, windowinc, bandwidth, [f_low f_high]);
            fig_ind = fig_ind+1;
        end

        if plot_vec(2)
            figure(f);
%             subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
%             %self.plot_polar_rate_map(current_cell);
%             fig_ind = fig_ind+1;
        end

        if plot_vec(3)
            figure(f);
            subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
            self.plot_lfp_spectrum(TBP, tapers,[f_low f_high], 0);
            fig_ind = fig_ind+1;
        end

        if plot_vec(4)
            figure(f);
%             subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
%             %self.plot_isi_dist(current_cell);
%             fig_ind = fig_ind+1;
        end

    end

end

function loaded_lfp = CheckLoadedLFP(self)

loaded_lfp = zeros(length(self.b_lfp), 1);

for i = 1:length(self.b_lfp)
    
    loaded_lfp(i) = ~isempty(self.b_lfp(i).signal);
    
end

end
    