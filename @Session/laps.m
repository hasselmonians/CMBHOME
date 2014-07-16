function laps = laps(self)
% root = root.laps;
%
% This function parses data recorded on mark's circle track into laps.  It
% accepts the CMBHOME Session class, and finds times corresponding to
% circle start and circle stop
%
% andrew january 6 2010
% andrew june 1 2010 updated for CMBHOME

self.epoch = [self.Label2Time('circle start'), self.Label2Time('circle stop')];

x = self.x;
y = self.y;
t = self.ts;

center = self.user_def.circle_center;

ind_nan = unique(find(isnan(x) | isnan(y)));

while ~isempty(ind_nan) % replaces all NaNs with the most recent non NaN
    x(ind_nan) = x(ind_nan-1);
    y(ind_nan) = y(ind_nan-1);
    ind_nan = unique(find(isnan(x) | isnan(y)));
end

pos = [x,y];
clear x y

theta(1:size(pos,1)) = NaN;

ind_q1 = find(pos(:,1)-center(1)>0 & pos(:,2)-center(2)<=0); % quadrant I 

ind_q2 = find(pos(:,1)-center(1)>=0 & pos(:,2)-center(2)>0); % quadrant II

ind_q3 = find(pos(:,1)-center(1)<0 & pos(:,2)-center(2)>=0); % quadrant III

ind_q4 = find(pos(:,1)-center(1)<=0 & pos(:,2)-center(2)<0); % quadrant IV

theta(ind_q1) = atan((pos(ind_q1,1)-center(1))./(center(2)-pos(ind_q1,2)));

theta(ind_q2) = pi/2 + pi/2 - atan((pos(ind_q2,1)-center(1))./(pos(ind_q2,2)-center(2)));

theta(ind_q3) = pi + atan((center(1)-pos(ind_q3,1))./(pos(ind_q3,2)-center(2)));

theta(ind_q4) = 1.5*pi + pi/2 - atan((center(1)-pos(ind_q4,1))./(center(2)-pos(ind_q4,2)));

% Find the first lap by looking for the first time that theta passes 290
% degrees, or 5.06145 radians

% 250 degrees, or 4.36 radians marks the end of the lap

% finds where theta crosses 290, indexes until it crosses 250

% check where theta is greater than start, or 5.06145

start = 5.06145;
stop = 4.36;

current_lap = 1;
status = 'waiting_to_reach_290';
for n = 1:length(theta);
    switch status
        case 'waiting_to_reach_290';
            if theta(n) >= 5.06145;
                if abs(theta(n)-theta(n-1))<5.06145;
                    status = 'waiting_to_reach_250';
                    laps(current_lap,1) = t(n);
                end
            end

        case 'waiting_to_reach_250';
            if theta(n) >= 4.36 && theta(n)<5.06145;
                if abs(theta(n)-theta(n-1))<4.36;
                    status = 'waiting_to_reach_290';
                    laps(current_lap,2) = t(n);
                    current_lap = current_lap+1;
                end
            end
    end
end

% clear zeros
ind = [find(laps(:,1)==0) find(laps(:,2)==0)];

laps(ind,:)=[];

%clear very short or very long laps (laps shorter than 5 seconds, or longer
%than 11
laps(laps(:,2)-laps(:,1)<5 | laps(:,2)-laps(:,1)>11,:)=[];