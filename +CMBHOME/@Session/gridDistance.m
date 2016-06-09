function [d, ax] = gridDistance(self,Pd,ifplot)
% Calculates and returns the distance from the center Acorr peak to the
%
% Inputs: self -- CMBObject, function assumes self.cel is set
%         Pd   -- "PeakDistance", the minimum distance between peaks.
%         Defaults to 7, which works for most cells. Change if issues.
%         ifplot -- For debugging mostly. Defaults to 0
%
% Outputs: d -- The 6 distances (pixels).

import CMBHOME.*

if ~exist('Pd','var');Pd = 7;end

rm = self.RateMap(self.cel);
rmA = CMBHOME.Utils.moserac(rm);
%%
[~,maxInds] = CMBHOME.Utils.extrema2(rmA); %linear index
[rowInd colInd] = ind2sub(size(rmA),maxInds);

rr = rowInd - rowInd(1);
cr = colInd - colInd(1);

d = sqrt(rr.^2 + rr.^2);

%% Get rid of false peaks
locs = sqrt(rowInd.^2+colInd.^2);

peakHeight = NaN(length(rowInd),1);
for i = 1:length(rowInd)
    peakHeight(i) = rmA(rowInd(i),colInd(i));
end

[peakHeight inds] = sort(peakHeight,'descend'); % May be unnecessary, seems like they come out in order
rowInd = rowInd(inds);
colInd = colInd(inds);


idelete = zeros(size(peakHeight))<0; 
for i = 1:length(peakHeight)
    if ~idelete(i)
        if i >1
            for k = 1:i-1 %check to see if too close to any of the previous (larger) peaks
                d(i,k) = sqrt((rowInd(i) - rowInd(k)).^2 + (colInd(i) -colInd(k)).^2);
                if d(i,k) < Pd;
                    idelete(i) = idelete(i)+1;
                end
            end
        end
    end
end
inds(idelete) = []; %indeces into rowInds and colInd 

%% Filter
rowInd = rowInd(inds);
colInd = colInd(inds);
d = sqrt((rowInd - rowInd(1)).^2 + (colInd - colInd(1)).^2);
[d inds] = sort(d,'ascend');

d = d(2:7);
%d = d./2;



%% Axes
y = rowInd(inds(2:7)) - rowInd(inds(1));
x = colInd(inds(2:7)) - colInd(inds(1));

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