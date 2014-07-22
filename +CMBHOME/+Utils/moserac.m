function CC = moserac(X, X2, threshold)
% AC = CMBHOME.Utils.moserac(M1, M2, thresh)
%
% Returns 2-d autocorrelation of matrix M with normalization standards as used by Mosers
%
% ARGUMENTS
%   M1           N x M matrix, like a ratemap
%   M2           (optional) O x P matrix, like a ratemap
%   threshold    (optional) value for pixel limit of edge effects (no lags
%                with less than 'threshold' pixels are eliminated
%                default=25.
%
% RETURNS
%   AC           2N-1 x 2M-1 matrix, Autocorrelation of M
%
% andrew august sometime 2010
% v2 - added cross correlation support 12/25/2010

if ~exist('threshold', 'var'), threshold = 25; end

if ~exist('X2', 'var'), X2 = X; end

n = conv2(ones(size(X)), ones(size(X2))); % number of bins at each lag

CC = (n.*conv2(X, rot90(conj(X2),2))) - (conv2(X, ones(size(X2))) .* conv2(fliplr(flipud(X2)), ones(size(X)))); % calculate 2d autocorr normalized like linear pearsons correlations coeff

CC = CC ./ ( sqrt( n.*conv2(X.^2, ones(size(X2))) - conv2(X, ones(size(X2))).^2 ) .* ...
             sqrt( n.*conv2(fliplr(flipud(X2)).^2, ones(size(X))) - conv2(fliplr(flipud(X2)), ones(size(X))).^2) );

CC(n<threshold) = NaN; % edge effects

