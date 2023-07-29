function [ model] = FNLS( Signal, segmentLength, maxNoHarmonics, f0Bounds, bool )
% for simply using fast fundamental freqeuncy estimation
%   need folder: fastF0Nls
    model = struct();
    if nargin == 4
        f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds);
    elseif nargin == 5
        f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds, bool);
    end
    
    [f0Estimate]= f0Estimator.estimate(Signal);
    period = round(1/f0Estimate);
    model.period = period;
    
end


