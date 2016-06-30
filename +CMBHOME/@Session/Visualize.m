function Visualize(self)
% Displays a GUI for browsing through typical first-pass plots
%
% Allows user to set active epoch from within GUI, but starts with epochs
% set in root.epoch.
%
% Only cells that meet root.cell_thresh are available at each instance.
%
% root.Visualize


    if isempty(self.spatial_scale)
        self.spatial_scale = .3;
        disp('root.spatial_scale was not set. Defaulting to .3 cm/pixel');
    end
    
    % set plotting params, which will display off to the side to be shown
    % when screen is resized by checkbox

    xcorr_binsize = 5; % ms
    ratemap_binside = 3; % cm
    velocity_thresh = [-1 -1]; % cm/sec
    plot_vec = zeros(25, 1); % which plots are active
    plot_vec(14) = 1; % header information is active
    theta_skip = 0; % run theta skipping analysis with spike time autocorr
    current_cell = self.cells(1,:);
        
    width = 650;
    height = 875;
        
    % Create and then hide the Figure Screen as it is being constructed
    f = figure('Visible','off','Position', [50, -150, width, height],'Color', 'w', 'DeleteFcn', {@deletef});
    
    
    control_height = 600;
    control_width = 250;
    
    %  Create and then hide the GUI as it is being constructed.
    control_f = figure('Visible','off','Position', [50, -150, control_width, control_height],'Color', 'w', 'Menubar', 'none', 'DeleteFcn', {@deletecontrolf});
    
    check_size = [225, 25];
    check_posx = repmat(15, 18, 1);
    
    check_ind = 1;
    
    %%  Construct the components.
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Title Information',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@title_Callback}, 'HandleVisibility', 'off', 'Value', 1);
    check_ind = check_ind+1;  
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Trajectory',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@trajectory_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    
    
    if length(self.b_headdir)==length(self.b_x)
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Polar Rate Map',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@polarratemap_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    
    end
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Rate Map',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@ratemap_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','ISI Histogram',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@isihist_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Spike Time Auto Correlation',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@autocorr_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    if length(self.b_headdir)==length(self.b_x)
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Head Direction',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@headdir_Callback}, 'HandleVisibility', 'off');  
    check_ind = check_ind+1;  
    end
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Spike Raster',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@raster_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;
    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Rate Map Autocorrelogram',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@rate_map_ac_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Gridness Autocorrelogram',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@gridness_ac_Callback},'HandleVisibility', 'off'); 
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Gridness3 Score',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@gridness_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Gridness2 Score',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@gridness2_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1; 

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Speed vs. Firing Frequency',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@velocityrate_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;
       
      
%% THIRD COLUMN

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Intrinsic Frequency',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@intrinsicf_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Precession (grid)',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility','off','Callback',{@precgrid_callback});
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Precession (place)',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off','Callback',{@precplace_callback});  
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','thetaModGamma',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off','Callback',{@thetaModGamma_callback});
    check_ind = check_ind+1;
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','waveform',...
          'Position',[check_posx(check_ind), control_height - 125 - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off','Callback',{@waveform_callback});
    check_ind = check_ind+1;
 
    
%% END OF CHECKBOXES

    htitle = uicontrol('Style','text','String','',...
          'Position',[15, control_height-25,300,25], 'BackgroundColor', 'w',...
          'FontSize', 16, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off'); 

    hspecs = uicontrol('Style','text','String','',...
          'Position',[15,control_height-115,300,60], 'BackgroundColor', 'w',...
          'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');       

    hpopup = uicontrol('Style','popupmenu',...
          'String',cellstr(num2str(self.cells)),...
          'Position',[15, control_height-55, 100, 25],...
          'Callback',{@popup_menu_Callback}, 'HandleVisibility', 'off'); 
      
    hdisplayparams = uicontrol('Style','checkbox','String','show parameters', 'fontsize', 7,...
                  'Position',[115, control_height-50, 100, 25],...
                  'BackgroundColor', 'w','Callback',{@displayparams_Callback}, 'HandleVisibility', 'off');

    hupdate = uicontrol('Style','pushbutton','String','Update',...
          'Position',[150, control_height-120, 70, 25],...
          'Callback',{@update_Callback}, 'HandleVisibility', 'off');
    
    hepoch = uicontrol('Style','pushbutton','String','Change Epochs',...
          'Position', [15, control_height-120, 125, 25],...
          'Callback', {@epoch_Callback}, 'HandleVisibility', 'off');
      
%% PARAMETER Controls

paramprops.BackgroundColor = 'w';
paramprops.FontSize = 8;

thetaskip_obj = uicontrol('Style','checkbox','String','Calculate Theta Skipping with Spike Time AutoCorr', paramprops,...
                  'Position',[250, 575,200, 25],...
                  'Callback',{@thetaskip_Callback});
              

uicontrol(paramprops, 'Style','text','Position',[240, 523, 100, 25],'String', 'Rate Map Binsize');
ratemap_label = uicontrol(paramprops, 'Style','text','Position',[370, 523, 50, 25],'String', ['x ' num2str(ratemap_binside) ' cm']);
          
ratemap_binside_obj = uicontrol('Style','edit',... % one textbox
                  'BackgroundColor', 'w', 'Callback',{@ratemap_binside_Callback},...
                  'Position',[330, 530, 50, 25],'String', num2str(ratemap_binside));

uicontrol(paramprops, 'Style','text','Position',[240, 495, 100, 25],'String', 'Time Correlation Binsize');
uicontrol(paramprops, 'Style','text','Position',[370, 490, 50, 25],'String', 'ms');
                        
xcorr_binsize_obj = uicontrol('Style','edit',... % one textbox
                  'BackgroundColor', 'w','Callback', {@xcorr_binsize_Callback}, ...
                  'Position',[330, 497, 50, 25],'String', num2str(xcorr_binsize));

uicontrol(paramprops, 'Style','text','Position',[250, 435, 225, 50], 'HorizontalAlign', 'left','String', 'The velocity bounds below apply to the intrinsic frequency and spike time autocorrelation (enter -1 for no thresholding).');

uicontrol(paramprops, 'Style','text','Position',[247, 413, 40, 25],'String', 'Min Speed');
uicontrol('Style','edit',... % one textbox
                  'BackgroundColor', 'w', 'Callback',{@low_vel_Callback},...
                  'Position',[295, 410, 50, 25],'String', num2str(velocity_thresh(1)));
              
uicontrol(paramprops, 'Style','text','Position',[360, 413, 40, 25],'String', 'Max Speed');
uicontrol('Style','edit',... % one textbox
                  'BackgroundColor', 'w', 'Callback',{@high_vel_Callback},...
                  'Position',[400, 410, 50, 25],'String', num2str(velocity_thresh(2))); 
   
   % Initialize the GUI.
   % Change units to normalized so components resize 
   % automatically.
   
    %set([control_f; hpopup; hupdate; check_objarr(:); htitle; hspecs; hepoch]);

    set(f, 'Units','normalized');

    str_title = ['Tetrode ' int2str(current_cell(1)) ' Cell ' int2str(current_cell(2))];
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

    if size(self.epoch, 1)>1
        N_spikes = self.cel_ts;
        N_spikes = length(vertcat(N_spikes{:}));
    else
        N_spikes = length(self.cel_ts);
    end

    F = N_spikes / sum(self.epoch(:,2)-self.epoch(:,1));

    str3 = ['Spikes/Sec during this event ' num2str(F) ' F'];

    str_specs = str2mat(str1, str2, str3);
    set(hspecs, 'String', str_specs);
    
    %ha = axes('Position', [.1 .1 .6 .6]);
   
%    %Create a plot in the axes.
    % Move the GUI to the center of the screen.
    movegui(control_f,'center')
    % Assign the GUI a name to appear in the window title.
    set(control_f,'Name','Visualize2')
    % Make the GUI visible.
    set(control_f,'Visible','on'); 
    
 
   %  Callbacks for simple_gui. These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 
   %  Pop-up menu callback. Read the pop-up menu Value property
   %  to determine which item is currently displayed and make it
   %  the current data.
   
    function deletecontrolf(source, event)
    % clears all figures (in case any are hidden)
        delete(f)
    end
   
    function deletef(source, event)
    % when user closes figure window, make a new hidden one. that way
    % visualize2 isnt closed until the control window is deleted
        delete(control_f)     
    end

    function high_vel_Callback(source, event)
            str = get(source, 'String');
            if ~isempty(str)  
                set(source, 'String', str);
                velocity_thresh(2) = str2double(str);
            else
                set(source, 'String', num2str(velocity_thresh(2)));
            end
    end

    function low_vel_Callback(source, event)
            str = get(source, 'String');
            if ~isempty(str)
                set(source, 'String', str);
                velocity_thresh(1) = str2double(str);
            else
                set(source, 'String', num2str(velocity_thresh(1)));
            end
    end   

    function ratemap_binside_Callback(source, eventdata)
            str = get(source, 'String');
            if ~isempty(str) && str2double(str)>0, 
                set(source, 'String', str);
                set(ratemap_label, 'String', ['x ' str ' cm']);
                ratemap_binside = str2double(str);
            else               
                set(source, 'String', num2str(ratemap_binside));
            end
    end

    function xcorr_binsize_Callback(source, eventdata)       
            str = get(source, 'String');            
            if ~isempty(str) && str2double(str)>0,               
                set(source, 'String', str);
                xcorr_binsize = str2double(str);          
            else              
                set(source, 'String', num2str(xcorr_binsize));
            end
    end
     
    function thetaskip_Callback(source, eventdata)
        if (get(source,'Value') == get(source,'Max'))
           theta_skip = 1;
        else
           theta_skip = 0;
        end
    end
   
    function popup_menu_Callback(source, eventdata) 
          
        % Determine the selected data set.
        % clear old plots, and update header

        str = get(source, 'String');
        val = get(source,'Value');
        
        % Set current data to the selected data set.
       
        current_cell = str2num(str{val});
        
        % Update header title
        
        h = get(f,'Children');
        
        delete(h)
               
        updateScreen(self, current_cell, plot_vec, f, control_f);
        
        
        figure(control_f) % bring focus back to control panel

    end

%% check box callbacks. each one increments n_plots by 1 and plot the given
%% dealio

    function displayparams_Callback(source, eventdata)
       % resizes control panel and displays the follow parameter list:
       %
       % autocorrelation binsizes, ratemap binsize, velocity thresholds 

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           poss = get(control_f, 'Position');
           set(control_f, 'Position', poss.*[1 1 2 1]); % make control panel twice as wide
        else
           % Checkbox is not checked-take approriate action
           poss = get(control_f, 'Position');
           set(control_f, 'Position', poss.*[1 1 .5 1]); % make control panel twice as wide
        end
        
    end

    function epoch_Callback(source, eventdata)
        self = self.SetEpoch;
        
        h = get(f,'Children');
        
        delete(h)
        
        updateScreen(self, current_cell, plot_vec, f, control_f);
    end        

    function title_Callback(source,eventdata) 
   % Display surf plot of the currently selected data.

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(14) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(14) = 0;
        end
      
    end

    function trajectory_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(3) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(3) = 0;
        end
          
    end

    function polarratemap_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(2) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(2) = 0;
        end   
    end 

    function ratemap_Callback(source,eventdata) 
    % Display surf plot of the currently selected data.

        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(1) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(1) = 0;
        end
      
    end

    function isihist_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(4) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(4) = 0;
        end
          
    end

    function autocorr_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(5) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(5) = 0;
        end
       
    end

    function headdir_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(6) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(6) = 0;
        end

    end

    function raster_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(7) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(7) = 0;
        end

    end


    function rate_map_ac_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(9) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(9) = 0;
        end

    end

    function gridness_ac_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.

        if get(source+1, 'Value') ~= get(source, 'Value')

           set(source+1, 'Value', get(source, 'Value'));

        end
       
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(10) = 1;
           plot_vec(11) = 1;

        else
           % Checkbox is not checked-take approriate action
           plot_vec(10) = 0;
           plot_vec(11) = 0;
        end

    end

    function gridness_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        %if get(source-1, 'Value') ~= get(source, 'Value')

        %   set(source-1, 'Value', get(source, 'Value'));

        %end
        
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(10) = 1;
           plot_vec(11) = 1;

        else
           % Checkbox is not checked-take approriate action
           plot_vec(10) = 0;
           plot_vec(11) = 0;
        end

    end

    function velocityrate_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(12) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(12) = 0;
        end
    end

    function intrinsicf_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(15) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(15) = 0;
        end
    end

    function gridness2_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(16) = 1;
           plot_vec(17) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(16) = 0;
           plot_vec(17) = 0;
        end
    end

    function precgrid_callback(source,eventdata)
         if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(18) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(18) = 0;
        end
    end

    function precplace_callback(source,eventdata)
         if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(19) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(19) = 0;
        end
    end

    function thetaModGamma_callback(source,eventdata)
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(20) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(20) = 0;
        end
          
    end

    function waveform_callback(source,eventdata)
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(21) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(21) = 0;
        end
          
    end 

    function update_Callback(source, eventdata)  
        
        hand = get(f,'Children');
        
        delete(hand)     
        
        updateScreen(self, current_cell, plot_vec, f, control_f);
               
    end

    function updateScreen(self, current_cell, plot_vec, f, control_f)
            
            import CMBHOME.*
            self.cel = current_cell;
            figure(control_f);
            
            rate_map = []; % initialize reused variables
        
            str_title = ['Tetrode ' int2str(current_cell(1)) ' Cell ' int2str(current_cell(2))];
            set(htitle, 'String', str_title);

            % Update header specs

            if size(self.epoch,1)>1
                str1 = [int2str(size(self.epoch,1)) ' Epochs Selected'];
                str2 = '';
            else
                str1 = ['Event Start: ' self.active_event{1}];
                str2 = ['Event Stop: ' self.active_event{2}];
            end

            if size(self.epoch, 1)>1
                N_spikes = CMBHOME.Utils.ContinuizeEpochs(self.cel_ts);
                N_spikes = length(N_spikes);
            else
                N_spikes = length(self.cel_ts {1});
            end

            F = N_spikes / sum(self.epoch(:,2)-self.epoch(:,1));

            str3 = ['Spikes/Sec during this event ' num2str(F) ' F'];

            str_specs = str2mat(str1, str2, str3);
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
                
                str_title = ['Tetrode ' int2str(current_cell(1)) ' Cell ' int2str(current_cell(2))];
                
                
                % Update header specs

                if size(self.epoch,1)>1
                    str1 = [int2str(size(self.epoch,1)) ' Epochs Selected'];
                    str2 = '';
                else
                    str1 = ['Event Start: ' self.active_event{1}];
                    str2 = ['Event Stop: ' self.active_event{2}];
                end

                if size(self.epoch, 1)>1
                    N_spikes = CMBHOME.Utils.ContinuizeEpochs(self.cel_ts);
                    N_spikes = length(N_spikes);
                else
                    N_spikes = length(self.cel_ts{1});
                end

                F = N_spikes / sum(self.epoch(:,2)-self.epoch(:,1));

                str3 = ['Spikes/Sec during this event ' num2str(F) ' F'];

                str4 = ['Num Spikes:' num2str(N_spikes)];
                
                str5 = self.name_formatted;
                str5 = strrep(str5, '_', '\_');
                
                str_specs = str2mat(str1, str2, str3, str4, str5);
                
                text(.1, .65, ['\fontsize{14}' str_title ' ' , sprintf('\n') ...
                    '\fontsize{8}' str_specs(1,:) ' ' , sprintf('\n') ...
                    '\fontsize{8}' str_specs(2,:) ' ' , sprintf('\n') ...
                    '\fontsize{8}' str_specs(3,:) ' ' , sprintf('\n') ...
                    '\fontsize{8}' str_specs(4,:) ' ' , sprintf('\n') ...
                    '\fontsize{8}' str_specs(5,:) ' ' , sprintf('\n') ...
                ]);
                
                %{
                text(.1, .8, [str_title ' '], 'FontSize', 14);
                
                text(.1, .6, str_specs, 'FontSize', 8);
                
                text(.1, .5, ['Enclosure size: ' num2str((max(self.b_x) - min(self.b_x))*self.spatial_scale, 3) ' x ' ...
                    num2str((max(self.b_y) - min(self.b_y))*self.spatial_scale, 3) ' cm^2'], 'FontSize', 8)
                
                text(.1, .4, strrep(self.name, '_', '\_'), 'FontSize', 8, 'Interpreter', 'none');
                %}
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(1)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                rate_map = self.RateMap(current_cell, 'continuize_epochs', 1, 'supress_plot', 0, 'binside', ratemap_binside);
                fig_ind = fig_ind+1;
            end

            if plot_vec(2)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_polar_rate_map(current_cell);
                fig_ind = fig_ind+1;
            end

            if plot_vec(3)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_trajectory(current_cell);
                fig_ind = fig_ind+1;
            end

            if plot_vec(4)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_isi_dist(current_cell);
                fig_ind = fig_ind+1;
            end

            if plot_vec(5)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.AutoCorr(current_cell, 't_bin', xcorr_binsize*10^-3,...
                                'max_lag', .6, 'speed_thresh', velocity_thresh,...
                                'theta_skip', theta_skip, 'supress_plot', 0);
                fig_ind = fig_ind+1;
            end

            if plot_vec(6)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_head_dir;
                fig_ind = fig_ind+1;
            end

            if plot_vec(7)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_raster(current_cell);
                fig_ind = fig_ind+1;
            end


            if plot_vec(8)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_theta_phase_dist(current_cell);
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(9)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_rate_map_ac(current_cell, rate_map);
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(10)
                figure(f);
                [gridness_score, props] = self.Gridness(current_cell, 'binside', ratemap_binside, 'continuize_epochs', 1, 'grid3', 1);
                if ~isnan(gridness_score) % if analysis valid
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                    clims = [min(min(props.auto_corr)) max(max(props.auto_corr))];

                    pad = [-.03 .02]; % percent pad plot

                    [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');

                    ac = props.auto_corr;
                    
                    ac(~props.auto_corr_mask) = clims(1);

                    imagesc(ac, clims); hold on;

                    colormap(cbar);
                    axis equal
                    axis off

                    title('Affected Autocorrelogram'), 

                    %set(tmp, 'AlphaDataMapping', 'none');   %   set cut out to white
                    %set(gca, 'DrawMode', 'fast');
                    %set(tmp,'AlphaData', cut_ac(:,:,2));

                    [ymax, xmax] = size(props.auto_corr); % set axis limits

                    xlim([1, xmax]);
                    ylim([1, ymax]);

                    fig_ind = fig_ind+1;

                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                    line(props.periodicity(:,1), props.periodicity(:,2), 'Color', 'k', 'LineWidth', 1.5), hold on;

                    text(100, -.4, ['Gridness3: ' num2str(gridness_score)], 'FontSize', 8); 

                    title('Periodicity of Correlation of Rotated AutoCorr');

                    xlabel('Rotation Angle'); ylabel('Correlation');

                    xlim([0 180]);
                    ylim([-.5, 1.1]);

                    fig_ind = fig_ind+1;
                else
                    
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                    
                    axis off
                    
                    text(.1, .4, 'Invalid Gridness: No Peaks in Autocorrelogram');
                    
                    fig_ind = fig_ind+1;             
                    
                end
            end
            
            if plot_vec(12)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.plot_velocityrate(current_cell, 1);
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(15)
                figure(f);
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                self.IntrinsicFrequency(current_cell, 1, velocity_thresh); % self.IntrinsicFrequency(current_cell, screen_print, speed_thresh cm/s)
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(16)
                figure(f);
                [gridness_score, periodicity, cut_ac] = self.Gridness2(current_cell, 'binside', ratemap_binside, 'continuize_epochs', 1);
                if ~isnan(gridness_score) % if analysis valid
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                    clims = [min(min(cut_ac(:,:,1))) max(max(cut_ac(:,:,1)))];

                    pad = [-.03 .02]; % percent pad plot

                    [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');

                    ac = cut_ac(:,:,1);
                    
                    ac(~cut_ac(:,:,2)) = clims(1);

                    imagesc(ac, clims); hold on;

                    colormap(cbar);
                    axis equal
                    axis off

                    title('Affected Autocorrelogram'), 

                    %set(tmp, 'AlphaDataMapping', 'none');   %   set cut out to white
                    %set(gca, 'DrawMode', 'fast');
                    %set(tmp,'AlphaData', cut_ac(:,:,2));

                    [ymax, xmax] = size(cut_ac(:,:,1)); % set axis limits

                    xlim([1, xmax]);
                    ylim([1, ymax]);

                    fig_ind = fig_ind+1;

                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                    line(periodicity(:,1), periodicity(:,2), 'Color', 'k', 'LineWidth', 1.5), hold on;

                    text(100, -.4, ['Gridness2: ' num2str(gridness_score)], 'FontSize', 8); 

                    title('Periodicity of Correlation of Rotated AutoCorr');

                    xlabel('Rotation Angle'); ylabel('Correlation');

                    xlim([0 180]);
                    ylim([-.5, 1.1]);

                    fig_ind = fig_ind+1;
                else
                    
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                    
                    axis off
                    
                    text(.1, .4, 'Invalid Gridness2: No Peaks in Autocorrelogram');
                    
                    fig_ind = fig_ind+1;             
                    
                end
            end
            
            if plot_vec(18)
                figure(f);
                self.cel = current_cell;
                self.active_lfp = self.cel(1);
                if isempty(self.active_lfp)
                   [~,self.active_lfp] = max(cellfun(@median,{self.b_lfp.theta_amplitude})); 
                end
                results = self.pass_index('plots','Scatter plot','subplots',{'Position',[a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]});
                temp = {'NOT ',''};
                title([temp{results.is_precessing+1} 'Precessing']);
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(19)
                figure(f);
                self.cel = current_cell;
                self.active_lfp = self.cel(1);
                if isempty(self.active_lfp)
                   [~,self.active_lfp] = max(cellfun(@median,{self.b_lfp.theta_amplitude})); 
                end
                results = self.pass_index('method','place','plots','Scatter plot','subplots',{'Position',[a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]});
                temp = {'NOT ',''};
                title([temp{results.is_precessing+1} 'Precessing']);
                fig_ind = fig_ind+1;
            end
            
            if plot_vec(20)
                figure(f);
                self.cel = current_cell;
                self.active_lfp = self.cel(1);
                
                if isempty(self.active_lfp)
                   [~,self.active_lfp] = max(cellfun(@median,{self.b_lfp.theta_amplitude})); 
                end
                
                subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                ax = gca;
                disp('Takes a while, please wait ...')
                self.thetaModGamma('ifPlot',ax);
                
                fig_ind = fig_ind + 1;
            end
            
            if plot_vec(21)
                %try
                    figure(f);
                    self.cel = current_cell;

                    ind = find(sum(repmat(self.cel,size(self.cells,1),1) == self.cells,2) == 2);
                    
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                    m = {self.user_def.waveform(ind,:).mean};
                    s = {self.user_def.waveform(ind,:).std};
                    
                    m=cellfun(@(x) x(:)', m,'UniformOutput',0);
                    s=cellfun(@(x) x(:)', s,'UniformOutput',0);
                    
                    hold on
                    for i = 1:length(m)
                        t = [(i-1)*length(m{i})+11:(i)*length(m{i})+10] + (i-1)*5;
                        t = [t fliplr(t)];
                        patch(t,[m{i}+s{i} fliplr(m{i}-s{i})],[.8 .8 .8],'EdgeColor',[.8 .8 .8]);
                        t = [(i-1)*length(m{i})+11:(i)*length(m{i})+10] + (i-1)*5;
                        plot(t,m{i},'r','LineWidth',3)
                    end

                    title('Waveform')
                    fig_ind = fig_ind + 1;
                    
                    ylim([min(cellfun(@(x) min(x),m)) max(cellfun(@(x) max(x),m))])
                    xlim([0 t(end)+5])
                %end
            end
                
        end

end 
