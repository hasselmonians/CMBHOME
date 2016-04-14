function epochs = SelectEpochs(epochs, validepochs, partials)
%
% CMBHOME.Utils.SelectEpochs(epochs, validepochs)
%
% Takes epochs that exist within valid epochs,
%
%   INPUTS
%       epochs - Nx2 array, epochs (time1 time2;...)
%       validepochs - Mx2 array, epochs
%       partials - (default 1) logical, 1==return parts of epochs within validepochs
%
%
% andrew 10/13/2011
% andrew 10/17/2011

if ~exist('partials', 'var'), partials = 1; end

if isempty(epochs), return; end

if isempty(validepochs), epochs = []; return; end

if any(epochs(:,2) - epochs(:,1)<0), error('epochs must be increasing from first column to second'); end

validepochs = CMBHOME.Utils.MergeEpochs2(validepochs, 0);

epochs = CMBHOME.Utils.MergeEpochs2(epochs, 0);

[~, inds] = sort(validepochs(:,1));

validepochs = validepochs(inds,:); % sort by starting times

S = repmat(validepochs(:,1), 1, size(epochs,1)); % start time matrix of valid epochs

E = repmat(validepochs(:,2), 1, size(epochs,1)); % start time matrix of valid epochs

eS = repmat(epochs(:,1)', size(validepochs,1), 1); % each column represents the i-th epochs

eE = repmat(epochs(:,2)', size(validepochs,1), 1);

if partials
    
    TF = (eS>=S & eS<=E) & ~(eE>=S & eE<=E); % epochs which start within a valid epochs, but dont end in a valid epoch (TAKE LAST ELEMENT in each column)
    
    [r, c] = find(TF); % where epochs start within a valid epoch, but the ending is nowhere
          
    [r, ind] = sort(r, 'ascend'); c = c(ind); % if we order everything so the largest rows are last, then the last value matlab will ascribe each column is the last row in that column

    epochsc = [epochs(c, 1), validepochs(r, 2)];
   
    TF = ~(eS>=S & eS<=E) & (eE>=S & eE<=E); % epochs which end within a valid epochs, but dont start on time (TAKE FIRST element in each column)

    [r, c] = find(TF);
    
    [r, ind] = sort(r, 'descend'); c = c(ind); % if we order everything so the largest rows are first, then the last value matlab will ascribe each column is the first row in that column

    epochsc2 = [validepochs(r, 1), epochs(c, 2)];
    
    epochssave = vertcat(epochsc, epochsc2);
 
    tmp = CMBHOME.Utils.SelectEpochs(validepochs, epochs, 0);
    
    tmp2 = CMBHOME.Utils.SelectEpochs(epochs, validepochs, 0);
    
    epochs = vertcat(tmp, tmp2, epochssave);

    epochs = unique(epochs, 'rows');
    
    return
    
end
    
TF = eS>=S & eE<=E; % matrix where epochs exist completely within a valid epochs (if 1 in any column, that epoch is OK with or without partial==1)

[r, c] = find(TF); % indices of valid epochs (entirely within validepochs)

epochs = epochs(c, :);

if any(sum(TF, 1)>=2), error('it appears that merge epochs (or something else) is not working in SelectEpochs.'); end % error checking

