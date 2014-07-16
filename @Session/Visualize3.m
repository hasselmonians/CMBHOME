function Visualize3(~)
% root.Visualize3
% 
% Used to visualize summary information about multiple cells/sessions next
% to each other.
%
% Asks for number of sessions, opens UI to select each session, prompts for
% each cell (in form of [1 1]) and then makes a new figure.
%
% To improve: 1.) Compeltely GUIize the system. 
            % 2.) Restrict Cell selections to the cells in the session
            % 3.) Allow user to select which parameters to plot
            % 4.) Improve Memory use (currently a huge memory hog)
            
%wchapman 2012.12.13


import CMBHOME.Utils.*

nsess = input('How many sessions to visualize from?: ');

% Load all of the sessions
for i  = 1:nsess
    [f2 f1] = uigetfile('*.mat',['select session ' num2str(i) ': ']);
    f = [f1 f2];
    
    c = load(f);
    roots(i) = c.root;  
end


roots2 = [];
for i = 1:nsess  
    nc = input(['How many cells for session' num2str(i) '?: ']);  
    
    for k = 1:nc      
        roots2 = [roots2;roots(i)];
        roots2(end).cel = input(['cell' num2str(k) ': ']);     
    end
end


numcells = length(roots2); %number of cells that we are plotting

% Ratemap, polarratemap, trajectory,isi,autocorr,r auto,grid auto grid
% score grid2 score, SPvsFR,Intrinsic.

holder = NaN(6,numcells);


figure
for i = 1:numcells
    % Ratemap
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[1 i]));
    roots2(i).plot_rate_map(roots2(i).cel);
    if iscell(roots2(i).name(end))
        title([roots2(i).name{end} num2str(roots2(i).cel)])
    else
        title([roots2(i).name(end) num2str(roots2(i).cel)])
    end
    
    % Polar ratemap
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[2 i]));
    roots2(i).plot_polar_rate_map(roots2(i).cel);
    title(roots2(i).HDWatsonsU2(roots2(i).cel))
    
    % Trajectory
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[3 i]));
    %roots2(i).plot_trajectory(roots2(i).cel)
    plot(h,roots2(i).x,roots2(i).y,'k')
    hold on
    plot(h,ContinuizeEpochs(roots2(i).cel_x),ContinuizeEpochs(roots2(i).cel_y{1}),'r.','MarkerSize',7);

    axis([min(roots2(i).x), max(roots2(i).x), min(roots2(i).y), max(roots2(i).y)])
    
    % Intrinsic Frequency2?
    %h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[6 i]));
    [it,t] = roots2(i).IntrinsicFrequency2(roots2(i).cel,'supress_plot',1);
    %title(['Intrinsic F2: ' num2str(it)])
    
    % Spike Autocorr
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[4 i]));
    roots2(i).AutoCorr(roots2(i).cel);
    title(['Thetaness: ' num2str(t)])
    
    % Spatial Autocorr
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[5 i]));
    rm = roots2(i).RateMap(roots2(i).cel);
    rm = moserac(rm);
    imagesc(rm);
    title(['Intrinsic Frequency: ' num2str(it) 'Hz'])
    
    % Speed Vs FR
    h = subplot(6,numcells,CMBHOME.Utils.SubplotSub2Ind(size(holder),[6 i]));
    roots2(i).plot_velocityrate(roots2(i).cel);
    

end

end