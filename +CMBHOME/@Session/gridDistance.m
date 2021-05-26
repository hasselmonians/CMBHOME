function [d, ax] = gridDistance(self,Pd,thresh,fromCenter,ifplot)
% Calculates and returns the distance from the center Acorr peak to the
%
% Inputs: self -- CMBObject, function assumes self.cel is set
%         Pd   -- "PeakDistance", the minimum distance between peaks.
%         Defaults to 7, which works for most cells. Change if issues.
%         thresh -- If nonempty, don't use peaks smaller than this value.
%         fromCenter -- If false (default) compute from maximum peak; else compute from center peak.
%         ifplot -- For debugging mostly. Defaults to 0

%
% Outputs: d -- The 6 distances (pixels).

import CMBHOME.*

if ~exist('Pd','var');Pd = 7;end
if ~exist('thresh', 'var') || isempty(thresh), thresh = -inf; end
if ~exist('fromCenter', 'var') || isempty(fromCenter), fromCenter = false; end

rm = self.RateMap(self.cel);
rmA = CMBHOME.Utils.moserac(rm);
%%
[~,maxInds] = CMBHOME.Utils.extrema2(rmA); %linear index
[rowInd colInd] = ind2sub(size(rmA),maxInds);

%% Get rid of false peaks

peakHeight = NaN(length(rowInd),1);
for i = 1:length(rowInd)
    peakHeight(i) = rmA(rowInd(i),colInd(i));
end

[peakHeight inds] = sort(peakHeight,'descend'); % May be unnecessary, seems like they come out in order

idelete = peakHeight < thresh;
for i = find(~idelete, 1):length(peakHeight)
    ind1 = inds(i);
    for k = 1:i-1 %check to see if too close to any of the previous (larger) peaks
        ind2 = inds(k);
        dist = sqrt((rowInd(ind1) - rowInd(ind2)).^2 + (colInd(ind1) -colInd(ind2)).^2);
        if dist < Pd
            idelete(i) = true;
        end
    end
end
inds(idelete) = []; %indeces into rowInds and colInd 

%% Filter
rowInd = rowInd(inds);
colInd = colInd(inds); % now ordered from max to min peak that is still above thresh

% identify the field for computing distance - either maximum or closest to center of arena
if fromCenter
    center = size(rmA) / 2;
    dFromCenter = vecnorm([rowInd(:), colInd(:)] - center, 2, 2);
    [~, ifield] = min(dFromCenter);
    rowField = rowInd(ifield);
    colField = colInd(ifield);
else
    rowField = rowInd(1);
    colField = colInd(1);
end

d = sqrt((rowInd - rowField).^2 + (colInd - colField).^2);
[d inds] = sort(d,'ascend');

d = d(2:7);
%d = d./2;



%% Axes
y = rowInd(inds(2:7)) - rowField;
x = colInd(inds(2:7)) - colField;

theta = atan2(x,y);

%{
[ax(1)] = min(theta(theta>0));
nx = wrapToPi(ax(1)+(2*pi/3));
[~,ind1] = min(abs(theta - nx));
ax(2) = theta(ind1);
nx = wrapToPi(ax(2)+(2*pi/3));
[~,ind1] = min(abs(theta-nx));
ax(3) = theta(ind1);
%}
ax = theta;
%% ifplot
colInd = colInd(inds(2:7));
rowInd = rowInd(inds(2:7));

if ~exist('ifplot','var');ifplot=0;end

if ifplot==1
    figure
    imagesc(rmA)
    hold on
    plot(colInd,rowInd,'ko','MarkerFaceColor',[0 0 0])
end

end