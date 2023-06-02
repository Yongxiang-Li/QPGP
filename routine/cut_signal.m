function [ shapeSignal ] = cut_signal( l, signal)
%CUT Summary of this function goes here
%   Detailed explanation goes here
    n = length(signal);
    if l >= n
        shapeSignal = signal;
    else
        shapeSignal = reshape(signal(1:l*floor(n/l)), l, floor(n/l));
    end
end