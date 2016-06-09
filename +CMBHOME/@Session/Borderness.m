function [borderness, Cm, Dm] = Borderness(self,cel,plotIt)
% Calculates the borderness score, as in Solstad 2008

% Bill 2012.10.26

% METHOD:
% Dm = Mean Distance to wall:
% Firing Field: FR >= 0.3*Max(FR).
% Area of field > 200cm^2
% 
%
%
% Cm = Max coverage of a wall
%   (Number of pixels in wall also in field) / (Number Pixels in wall)

%
% If multiple recording sessions for the cell, should discard if spatial
% correlation < 0.5
%
% [borderness, Cm, Dm] = root.Borderness(cel, ifPlot)

%% Setup
self.cel = cel;
rm = self.RateMap(self.cel);

sp_thr = max(max(rm)) / 3;

b = bwboundaries(rm);
b = b{1};

rm_bord = zeros(size(rm,1),size(rm,2));
for i = 1:length(b)
    rm_bord(b(i,1),b(i,2)) = 1;
end

field = (rm>=sp_thr);
a = regionprops(field,{'Area' 'PixelList'});

b = [];
for i = 1:length(a)
    if ((a(i).Area / self.spatial_scale) > 200)
        b = [b;a(i).PixelList];
    end
end

field = zeros(size(rm,1),size(rm,2));
for i = 1:length(b)
   field(b(i,2),b(i,1)) = 1; 
end  

%% If the cell is a grid cell, then can not be a border cell
%{
gr = self.Gridness2(self.cel);
if gr > 0.34
    borderness = NaN;
    Cm = NaN;
    Dm = NaN;
end
%}
%% Check for fields that are too large (interneurons)
%%{
[l w] = size(rm);
ul = (l*w) / 3;
if ((sum(sum(field)) > ul) || isempty(b))
    borderness = NaN;
    Cm = NaN;
    Dm = NaN;
    return
end
%}
%% Major length of field along wall
rm_bord2 = rm_bord .* field;
bord_west = sum(sum(rm_bord2(:,(1:3))) ) / sum(sum(rm_bord(:,(1:3)))); 
bord_east = sum(sum(rm_bord2(:,(size(rm,2)-2:size(rm,2))))) / sum(sum(rm_bord(:,(size(rm,2)-2:size(rm,2)))));
bord_north = sum(sum(rm_bord2(1:3,:))) / sum(sum(rm_bord(1:3,:))) ;
bord_south = sum(sum(rm_bord2(size(rm_bord,1)-2:size(rm_bord,1),:))) / sum(sum(rm_bord(size(rm_bord,1)-2:size(rm_bord,1),:)));

Cm = max([bord_west,bord_east,bord_north,bord_south]);

%% Mean Distance from field to border
[D,~] = bwdist(rm_bord);
rm3 = rm / (nanmax(nanmax(field.*rm)));
rm3 = rm3.*D.*field;
shortL = min([size(rm,1) size(rm,2)]) /2;
Dm = nanmean(rm3(rm3~=0)) / shortL;
%% Borderness
borderness = (Cm - Dm) / (Cm + Dm);


%% Plot
if ~exist('plotIt','var');plotIt = 0;end
if plotIt ==1
    figure
    subplot(6,1,1)
    imagesc(rm)
    subplot(6,1,2)
    imagesc(rm_bord)
    subplot(6,1,3)
    imagesc(field)
    subplot(6,1,4)
    imagesc(D)
    subplot(6,1,5)
    imagesc(D.*field)
    subplot(6,1,6)
    imagesc(D.*field.*rm)

%keyboard
%}
end