function [ model] = FNLS( Signal, segmentLength, maxNoHarmonics, f0Bounds )
% for simply using fast fundamental freqeuncy estimation
%   need folder: fastF0Nls
    model = struct();
    f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds);
    
    [f0Estimate]= f0Estimator.estimate(Signal);
    period = round(1/f0Estimate);
    model.period = period;
    
end


