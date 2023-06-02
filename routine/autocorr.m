function [ autocorr, amdf ] = autocorr( noiseSignal, x, k )
%AUTOCORR Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 3,    k = mean(noiseSignal.^2); end
    n = length(noiseSignal);
    autocorr = noiseSignal(1:n-x)*noiseSignal(x+1:n)';
    autocorr = autocorr/(n-x);
     
    % Weighted Autocorrelation for Pitch Extraction of Noisy Speech
    amdf = mean(abs(noiseSignal(1:n-x) - noiseSignal(x+1:n)));
    amdf = autocorr / (amdf+k);
end

