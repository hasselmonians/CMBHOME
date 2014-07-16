function arrout = DownSampleSum(arrin,binsize)

%% This function takes an array and sums its elements as per binsize. 
%% It trims off the cells at the end that do not constitute a full bin 

% andrew 8 oct 2009, inspired by neat trick somewhere in the www

arrin=arrin(1:end-rem(length(arrin),binsize));
numbins=length(arrin)/binsize;

if rem(numbins,1)~=0
    fprintf('somethings wrong in the binning function\n')
    return;
end

arrout=reshape(sum(reshape(arrin,binsize,numbins),1),numbins,1);
