function watsonsu2 = HDWatsonsU2(self, cel)
% watsonsu2 = root.HDWatsonsU2(cel)
%
% Returns U2 score for cel for all root.epoch. Vector if more than one
% epoch
%
% watsonsu2: Watson's U^2 test for uniformity in circular data. See
% redish,
% hippocampus 2005 to learn how it applies to head direction cells

% see
% http://www.ncbi.nlm.nih.gov/pmc/articles/mid/NIHMS71281/table/T2/
% for algorithm

% also see merriam and davis 2008 pg 10

% a common way to apply this to head direction/spiking data is to
% pass the array of head direction samples and spike head direction
% angles. the null hypothesis is thus that spikes head directions come 
% from the same distribution of all head direction angles. the null 
% hypothesis is rejected if U^2 is 'large'. >10 seems to be the 
% neuroscience HD cell literature standard. - andrew

% calculation by ehren newman 1 march 2010
% addition to CMBHOME by andrew abogaard 25 mat 2010
 
% remove any nan-s from x and y

import CMBHOME.*

epochs = self.epoch;

watsonsu2 = zeros(size(epochs,1), 1);

for i = 1:size(epochs, 1)
    
    self.epoch = epochs(i,:);
    
    watsonsu2(i) = Spike.WatsonsU2(self.spk_headdir(cel), self.headdir);
    
end