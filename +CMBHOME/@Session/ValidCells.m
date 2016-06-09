function cells = ValidCells(self)
% returns matrix 'cells', an Nx2 matrix of [Tetrode Number, Cell Number] of viable cells 

% if f_min is set, the number of total spikes/total recording time (not given epoch set time) must be >f_min
% for that cell to be included

% if checkdrift~=1 && exist('checkdrift'), then the script will check to
% make sure that drift has not occurred during the run by sampling down to
% 1Hz and making sure that for a sliding window of 10 seconds we see at
% least f_min*10 spikes throughout the run

% andrew mar 2010

if length(self.cell_thresh)==2
    
    f_min = self.cell_thresh(1);
    checkdrift = self.cell_thresh(2);
    
end

if isempty(self.spike)
    disp('No cells in session object');
    cells = [];
    return;
end

cells=[0,0]; % determine cells to run

for i=1:size(self.spike,1)
    for j=1:size(self.spike,2)
        if ~isempty(self.spike(i,j).ts), cells(end+1,:) = [i, j]; end
    end
end
            
cells(1,:) = [];            
            
if exist('f_min', 'var')
    
    inds=[];
    
    for i = 1:size(cells,1)
        
        if length(self.spike(cells(i,1), cells(i,2)).ts)/(self.b_ts(end)-self.b_ts(1)) > f_min
            
            if checkdrift==1
                
                count = histc(self.spike(cells(i,1), cells(i,2)).ts, self.b_ts(1):10:self.b_ts(end));
                
                count = count(1:end-1);

                if all(count/10 >= f_min)
                    inds(end+1) = i;
                end
                
            else
                
                inds(end+1) = i;
                
            end                    
        end
    end
    
   cells = cells(inds,:);
   
end



