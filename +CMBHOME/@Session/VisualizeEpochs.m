function VisualizeEpochs(self)
% root.Visualize2
%
% Displays a GUI for browsing through different cell analysis like rate
% maps, autocorrelations, and rasters.
%
% Allows user to set active epoch from within GUI, but starts with epochs
% set in root.epoch.
%
% Only cells that meet root.cell_thresh are available at each instance.
%
% andrew bogaard 3 april 2010

    epochs = self.epoch;

    if isempty(self.spatial_scale)
        self.spatial_scale = .3;
        disp('root.spatial_scale was not set. Defaulting to .3 cm/pixel');
    end
    
    width = 300;
    height = 875;
        
    % Create and then hide the Figure Screen as it is being constructed
    f = figure('Visible','off','Position', [50, -150, width, height],'Color', 'w');
    
    
    control_height = 600;
    control_width = 250;
    
    %  Create and then hide the GUI as it is being constructed.
    control_f = figure('Visible','off','Position', [50, -150, control_width, control_height],'Color', 'w', 'Menubar', 'none');
    
    check_size = [225, 25];
    check_posx = repmat(15, 18, 1);
    
    option_buffer = 160; % pixel buffer for checkboxes from top
    
    check_ind = 1;
    
    %  Construct the components.
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Title Information',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@title_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;  
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Rate Map',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@ratemap_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;  
    
    if length(self.b_headdir)==length(self.b_x)
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Polar Rate Map',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@polarratemap_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    
    end
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Trajectory',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@trajectory_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;  

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','ISI Histogram',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@isihist_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Spike Time Auto Correlation',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@autocorr_Callback}, 'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    if length(self.b_headdir)==length(self.b_x)
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Head Direction',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@headdir_Callback}, 'HandleVisibility', 'off');  
    check_ind = check_ind+1;  
    end
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Spike Raster',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@raster_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;
    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Rate Map Autocorrelogram',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@rate_map_ac_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Gridness Autocorrelogram',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@gridness_ac_Callback},'HandleVisibility', 'off'); 
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Gridness Score',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@gridness_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Speed vs. Firing Frequency',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@velocityrate_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1;
    
    check_objarr(check_ind) = uicontrol('Style','checkbox','String','Kalman Speed vs. Firing Frequency',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','Callback',{@velocityrate2_Callback},'HandleVisibility', 'off');
    check_ind = check_ind+1; 
   
      
%% THIRD COLUMN

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','empty',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','empty',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off');
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','empty',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off');  
    check_ind = check_ind+1;    

    check_objarr(check_ind) = uicontrol('Style','checkbox','String','empty',...
          'Position',[check_posx(check_ind), control_height - option_buffer - 25*check_ind, check_size],...
          'BackgroundColor', 'w','HandleVisibility', 'off');
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

    hupdate = uicontrol('Style','pushbutton','String','Update',...
          'Position',[150, control_height-120, 70, 25],...
          'Callback',{@update_Callback}, 'HandleVisibility', 'off');
    
    hepoch = uicontrol('Style','pushbutton','String','Change Epochs',...
          'Position', [15, control_height-120, 125, 25],...
          'Callback', {@epoch_Callback}, 'HandleVisibility', 'off');
      
    uicontrol('Style','text','String','Epoch Interval (min)',...
          'Position',[20,control_height-157,150,25], 'BackgroundColor', 'w',...
          'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
      
    hepochinterval = uicontrol('Style','edit',... % one textbox
          'BackgroundColor', 'w', ...
          'Position',[120,control_height-150, 100, 25],...
          'HandleVisibility', 'off', 'String', '');
      
%     htstart = uicontrol('Style', 'popupmenu', ...
%         'String', cat(1,{'Session Start'}, self.event{:,1}), ...
%         'Value', 1, 'Position', [350, height-115, 125, 20],...
%         'Callback', {@tstart_menu_Callback}, 'HandleVisibility', 'off');
%     
%     htstop = uicontrol('Style', 'popupmenu', ...
%         'String', cat(1, self.event{:,1}, {'Session Stop'}), ...
%         'Value', 1, 'Position', [500, height-115, 125, 20],...
%         'Callback', {@tstop_menu_Callback}, 'HandleVisibility', 'off');
                 
   %align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
   
%    % Create the data to plot.
%    peaks_data = peaks(35);
%    membrane_data = membrane;
%    [x,y] = meshgrid(-8:.5:8);
%    r = sqrt(x.^2+y.^2) + eps;
%    sinc_data = sin(r)./r;
   
   % Initialize the GUI.
   % Change units to normalized so components resize 
   % automatically.
   
    %set([control_f; hpopup; hupdate; check_objarr(:); htitle; hspecs; hepoch]);

    set(f, 'Units','normalized');

    set([control_f; hpopup; hupdate; check_objarr(:); htitle; hspecs; hepoch],...
    'Units','normalized');

    plot_vec = zeros(25, 1);
    current_cell = self.cells(1,:);

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
        N_spikes = self.spk_ts(current_cell);
        N_spikes = length(vertcat(N_spikes{:}));
    else
        N_spikes = length(self.spk_ts(current_cell));
    end

    F = N_spikes / sum(self.epoch(:,2)-self.epoch(:,1));

    str3 = ['Spikes/Sec during this event ' num2str(F) ' F'];

    str_specs = str2mat(str1, str2, str3);
    set(hspecs, 'String', str_specs);
    
    %ha = axes('Position', [.1 .1 .6 .6]);
   
%    %Create a plot in the axes.
%    current_data = peaks_data;
%    surf(current_data);
    % Move the GUI to the center of the screen.
    movegui(control_f,'center')
    % Assign the GUI a name to appear in the window title.
    set(control_f,'Name','VisualizeEpochs')
    % Make the GUI visible.
    set(control_f,'Visible','on'); 
    
 
   %  Callbacks for simple_gui. These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 
   %  Pop-up menu callback. Read the pop-up menu Value property
   %  to determine which item is currently displayed and make it
   %  the current data.
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

    end

%% check box callbacks. each one increments n_plots by 1 and plot the given
%% dealio

    function epoch_Callback(source, eventdata)
        self = self.SetEpoch;
        
        h = get(f,'Children');
        
        delete(h)
        
        updateScreen(self, current_cell, plot_vec);
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
       
        if get(source-1, 'Value') ~= get(source, 'Value')

           set(source-1, 'Value', get(source, 'Value'));

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

    function velocityrate2_Callback(source,eventdata) 
       % Display surf plot of the currently selected data.
       
        % Display surf plot of the currently selected data.
        if (get(source,'Value') == get(source,'Max'))
           % Checkbox is checked-take approriate action
           plot_vec(13) = 1;
        else
           % Checkbox is not checked-take approriate action
           plot_vec(13) = 0;
        end
    end

    function update_Callback(source, eventdata)  
        
        hand = get(f,'Children');
        
        delete(hand)     
        
        updateScreen(self, current_cell, plot_vec, f, control_f);
               
    end

    function updateScreen(self, current_cell, plot_vec, f, control_f)
            
            import CMBHOME.*
            
            epochs = get(hepochinterval, 'String');
            
            if ~isempty(epochs)
                
                epochs = str2double(epochs);
                
                epochs = self.b_ts(1):epochs*60:self.b_ts(end);
                
                epochs = epochs';
                
                if epochs(end)~=self.b_ts(end), epochs(end+1) = self.b_ts(end); end
                
                epochs = cat(2, epochs(1:end-1), epochs(2:end));
                
            end

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
                N_spikes = self.spk_ts(current_cell);
                N_spikes = length(vertcat(N_spikes{:}));
            else
                N_spikes = length(self.spk_ts(current_cell));
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
            
            set(f, 'Units', 'pixels');
            pos_tmp = get(f, 'Position');
            set(f, 'Position', [pos_tmp(1:2), width*size(epochs,1), height]);
            set(f, 'Units', 'normalized');
            
            h_max = .95;
            h_min = .05; 
            buffer = .05;
            
            w_max = .95;
            w_min = .05;            

            % make as many subplots as necessary, display them

            a_h = (h_max-h_min-(sum(plot_vec)-1)*buffer)/sum(plot_vec);
            a_hs = linspace(h_max-a_h,h_min, sum(plot_vec));
            a_hs = repmat(a_hs', 1, size(epochs, 1));

            a_w = (w_max-w_min-(size(epochs,1)-1)*buffer)/size(epochs,1);
            a_ws = linspace(w_min, w_max-a_w, size(epochs,1));
            a_ws = repmat(a_ws, sum(plot_vec), 1);

            fig_ind = 1;
            
            for j = 1:size(epochs,1)

                self.epoch = epochs(j,:);
                
                %uicontrol('Type', 'text', 'Position', [a_ws(fig_ind), h_max, a_w, .05],...
                    %'String', ['From ' num2str(epochs(j,1)) ' to ' num2str(epochs(j,2)) ' minutes']);

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
                        N_spikes = self.spk_ts(current_cell);
                        N_spikes = length(vertcat(N_spikes{:}));
                    else
                        N_spikes = length(self.spk_ts(current_cell));
                    end

                    F = N_spikes / sum(self.epoch(:,2)-self.epoch(:,1));

                    str3 = ['Spikes/Sec during this event ' num2str(F) ' F'];

                    str_specs = str2mat(str1, str2, str3);

                    text(.1, .8, [str_title ' '], 'FontSize', 14);

                    text(.1, .6, str_specs, 'FontSize', 8);

                    text(.1, .4, strrep(self.name, '_', '\_'), 'FontSize', 8);

                    fig_ind = fig_ind+1;
                end

                if plot_vec(1)
                    figure(f);
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                    rate_map = self.plot_rate_map(current_cell);
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
                    self.plot_auto_corr(current_cell, .6, .020);
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
                    [gridness_score, periodicity, cut_ac] = self.Gridness(current_cell);
                    if gridness_score % if analysis valid
                        subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                        tmp = imagesc(cut_ac(:,:,1));

                        title('Affected Autocorrelogram'), 

                        set(gca, 'Box', 'on');

                        axis equal

                        set(tmp, 'AlphaDataMapping', 'none');   %   set cut out to white
                        set(gca, 'DrawMode', 'fast');
                        set(tmp,'AlphaData', cut_ac(:,:,2));

                        [ymax, xmax] = size(cut_ac(:,:,1)); % set axis limits

                        xlim([1, xmax]);
                        ylim([1, ymax]);

                        fig_ind = fig_ind+1;

                        subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),

                        line(periodicity(:,1), periodicity(:,2), 'Color', 'k', 'LineWidth', 1.5), hold on;

                        text(100, -.4, ['Gridness: ' num2str(gridness_score)], 'FontSize', 8); 

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
                    self.plot_velocityrate(current_cell);
                    fig_ind = fig_ind+1;
                end

                if plot_vec(13)
                    figure(f);
                    subplot('Position', [a_ws(fig_ind), a_hs(fig_ind), a_w, a_h]),
                    self.plot_velocityrate(current_cell, [1 1]);
                    fig_ind = fig_ind+1;
                end
                
            end

        end

end 