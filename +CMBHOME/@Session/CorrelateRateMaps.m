function R = CorrelateRateMaps(self, cel, obj2, cel2, supress_plot)
% (1) R = root.CorrelateRateMaps(cel);
% (2) R = root.CorrelateRateMaps(cel, obj2, cel2);
%
% (1) Calculates rate maps for cell in cel = [t, c]; in current root, prompts
% user to indicate another session object and cell to compare. If two cells
% are selected, compares both from same session
% (2) Calculates rate maps for cell cel in current session object, and cell
% "cel2" in session object "obj2", which must be loaded into the workspace.
%
% Finds the dimensions of the first rate map, and uses them to constrain
% the second. Also, removes values where occupancy was zero in either case
%
% To calculate rate map uses: root.plot_rate_map, and calculates the 
% pearsons correlation coefficient between the two
%
% andrew bogaard 9 sept 2010

if ~exist('supress_plot', 'var'), supress_plot = 0; end

if size(cel,1)==2
    obj2 = self; 
    cel2 = cel(2,:);
    cel = cel(1,:);
elseif ~exist('obj2', 'var')
    path = self.name(1:find(self.name=='/', 1, 'last'));
    [FileName,PathName] = uigetfile('*.mat','Select file containing second Session Object', path);
    obj2 = load(fullfile(PathName, FileName));
    obj2 = obj2.root;
    cel2(1) = input('Tetrode number: ');
    cel2(2) = input('Cell number: ');
elseif ~exist('cel2', 'var')
    cel2(1) = input('Tetrode number: ');
    cel2(2) = input('Cell number: ');
end

[rate_map1, ~, ~, xdim, ydim, occupancy1] = self.plot_rate_map(cel, 1); % rate map cell 1
[rate_map2, ~, ~, ~, ~, occupancy2] = obj2.plot_rate_map(cel2, 1, xdim, ydim); % rate map cell 2

occupancy = min(occupancy1, occupancy2);

R = corrcoef(rate_map1(occupancy), rate_map2(occupancy));

R = R(2);

if ~supress_plot % plot to compare alignment of rate maps
    subplot(1, 2, 1), imagesc(rate_map1), title('Cell 1');
    
    ys = ylim;
    xs = xlim;
    
    text(xs(2)+diff(xs)*.05, ys(2)+diff(ys)*.05, ['R: ' num2str(R)]);
    
    subplot(1, 2, 2), imagesc(rate_map2), title('Cell 2');
end   