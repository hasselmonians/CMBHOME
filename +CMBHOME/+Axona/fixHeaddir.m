function headdir = fixHeaddir(headdir)		

plotIt = 0;
shortEpochThresh = 50;
bigDeviationThresh = 100;
maxChunkSize2flip = 200; % larger starts flipping inappropriate data, smaller might be okay. . .

	if plotIt; hold off; plot(headdir(1:40000),'.b'); end

	%%% PHASE 1:
	% Drop all data points that deviate from the median HD of the surrounding
	% 2 seconds worth of data.
	devs = ones(size(headdir));
	for i = 1:length(headdir)
		inds = max([1 i-50]) : min([length(headdir) i+50]);
		devs(i) = min(abs(nanmedian(headdir(inds))-headdir(i) + [-360 0 360]));
	end
	bigDevs = abs(devs) > bigDeviationThresh;
	headdir(bigDevs) = nan;

	% fill in missing values as linear interpolation for gaps smaller than
	% defined threshold
	E = CMBHOME.Utils.OverThresholdDetect(isnan(headdir),0.5,1,1);
	shortEpochs = find(diff(E')<shortEpochThresh);
	flips2fix = [];
	for i = shortEpochs
		if E(i,1)==1 || E(i,2)==length(headdir), continue, end;
		span = circDiff([headdir(E(i,1)-1) headdir(E(i,2)+1)],2,'deg');
		if abs(span) > bigDeviationThresh, % if angle range to be corrected is too large, don't linearly interpolate
			if abs(span)>150,	flips2fix = [flips2fix i]; end % if it is very large, it probably needs to be flipped
			continue;
		end
		tmpVals = linspace(0,span,diff(E(i,:),[],2)+3);
		tmpVals = headdir(E(i,1)-1) + tmpVals;
		headdir(E(i,1):E(i,2)) = tmpVals(2:end-1);
	end

	INDS = [];% useful for debugging purposes

	% iteratively fix large flips of data
	while ~isempty(flips2fix)
		
		% find all possible locations of flips
		diffs = circDiff(headdir);
		bigDiffs = [find(abs(diffs)>130) E(flips2fix,1)'];
		nflips = length(flips2fix); % keep track of where we were at at the start
		while ~isempty(flips2fix)
			i = flips2fix(1);
			bigDiffs(bigDiffs==E(i,1)) = []; % drop current flip from list of possible flips
			if isempty(bigDiffs), break, end; % quit if out of possible flip points
			[flipOffset, bigDiffsInd] = minMag(bigDiffs-E(i,1),2); % find nearest possible flip point
			inds = sort([E(i,1)  E(i,1) + flipOffset]); % define indices to flip between
			INDS=[INDS; inds]; % save indices for later reference if needed
			if diff(inds)<maxChunkSize2flip, % Don't execute flips of huge chunks of data, they are less likely to be accurate
				headdir(inds(1):inds(2)) = mod(headdir(inds(1):inds(2)) + 180,360); % flip all data in window
				flips2fix(E(flips2fix,1)==bigDiffs(bigDiffsInd)) = []; % if used another flip listed to be fixed, drop it from list
				bigDiffs(bigDiffsInd) = []; % drop flip point that was used from list
			end
			flips2fix(E(flips2fix,1)==E(i,1)) = []; % drop this flip from list of flips to fix
		end

		% try to fill in nan gaps again now that flips have been corrected
		E = CMBHOME.Utils.OverThresholdDetect(isnan(headdir),0.5,1,1);
		shortEpochs = find(diff(E')<shortEpochThresh);
		flips2fix = [];
		for i = shortEpochs
			if E(i,1)==1 || E(i,2)==length(headdir), continue, end;
			span = circDiff([headdir(E(i,1)-1) headdir(E(i,2)+1)],2,'deg');
			if abs(span) > bigDeviationThresh,
				flips2fix = [flips2fix i];
				continue;
			end
			tmpVals = linspace(0,span,diff(E(i,:),[],2)+3);
			tmpVals = headdir(E(i,1)-1) + tmpVals;
			headdir(E(i,1):E(i,2)) = tmpVals(2:end-1);
		end
		headdir = mod(headdir,360);
		if(length(flips2fix)==nflips), break, end;

	end
	
	if plotIt; hold on; plot(headdir(1:40000),'.r');end;

end